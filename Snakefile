__author__ = "Benjamin Hillmann"
__license__ = "AGPL"

from snakemake.utils import min_version

min_version("3.11.2")

configfile: "config.yaml"

shogun_db_files = glob_wildcards("data/references/{database}/".format(database=config['database']) + "{files}")

uds_runs, sample_names = glob_wildcards("data/hiseq4000/{uds_run}/{sample_name}.fastq.gz")

rule all:
    input:
        expand("figs/fig{f}.pdf", f=[1,]),
        expand("results/uds/{uds_run}.{sample_name}.txt", uds_run=uds_runs, sample_name=sample_names)

############# Analysis ##############

rule extract_uds:
    input:
        "data/hiseq4000/{uds_run}/{sample_name}.fastq.gz"
    conda:
        "envs/shogun.yaml"
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}/{sample_name}.fastq")
    shell:
        "7z x {input} -so > {output}"

rule quality_control_uds:
    input:
        "/dev/shm/uds/{uds_run}.{sample_name}/{sample_name}.fastq"
    conda:
        "envs/shogun.yaml"
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
        expand("/dev/shm/{dir}/{file}", zip, dir=shogun_db_directories, file=shogun_db_files),
        queries = "/dev/shm/uds/{uds_run}.{sample_name}/combined_seqs.fna"
    benchmark:
        "results/benchmarks/{sample_name}.shogun.log"
    priority: 4
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}.txt")
    shell:
        "shogun align --input {params} "

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
        expand("/dev/shm/{dir}/{file}", zip, dir=shogun_db_directories, file=shogun_db_files),
        queries = "/dev/shm/uds/{uds_run}.{sample_name}/combined_seqs.fna"
    benchmark:
        "results/benchmarks/{sample_name}.shogun.log"
    priority: 4
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}.txt")
    shell:
        "shogun align --input {params} "



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