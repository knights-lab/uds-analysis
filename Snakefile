__author__ = "Benjamin Hillmann"
__license__ = "AGPL"

from snakemake.utils import min_version

min_version("3.11.2")

configfile: "config.yaml"

shogun_db_directories, shogun_db_files = glob_wildcards("data/references/{database}/{file}", database=config['database'])

rule all:
    input:
        expand("figs/fig{f}.pdf", f=[1,])

############# Analysis ##############
rule extract_uds:
    input:
        "data/hiseq4000/{uds_run}/{sample_name}.fastq.gz"
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}/{sample_name}.fastq")
    shell:
        "7z x {input} -so > {output}"

rule quality_control_uds:
    input:
        "/dev/shm/uds/{uds_run}.{sample_name}/{sample_name}.fastq"
    params:
        "/dev/shm/uds/{uds_run}.{sample_name}"
    priority: 1
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}/combined_seqs.fna"),
        temp("/dev/shm/uds/{uds_run}.{sample_name}/shi7.log"),
    shell:
        "shi7 -SE --combine_fasta True -i {params} -o {params} --adaptor Nextera -trim_q 32 -filter_q 36 --strip_underscore True -t 24"

rule place_on_ramdisk:
    input:
        "{dir}/{file}"
    priority: 2
    output:
        temp("/dev/shm/{dir}/{file}")
    shell:
        "cp {input} {output}"

rule shogun_filter:
    input:
        pass
    priority: 3

rule shogun_align:
    input:
        expand("/dev/shm/{dir}/{file}", zip, dir=shogun_db_directories, file=shogun_db_files),
        queries = "/dev/shm/uds/{uds_run}.{sample_name}/combined_seqs.fna"
    benchmark:
        "results/benchmarks/{sample_name}.shogun.log"
    priority: 4
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}.txt")
    shell:
        "shogun align --input {params} "

rule shogun_function:
    input:
        pass



################# Plots #################

########### Tables ##############

########### Figures #############
rule convert_svg:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.{fmt,(pdf|png)}"
    conda:
        "envs/cairosvg.yaml"
    shell:
        "cairosvg -f {wildcards.fmt} {input} -o {output}"