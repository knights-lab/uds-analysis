__author__ = "Benjamin Hillmann"
__license__ = "AGPL"

from glob import glob
from snakemake.utils import min_version

min_version("3.11.2")

configfile: "config.yaml"

#shogun_db_files = glob(config['database'] + "/**/*.*", recursive=True)
shogun_db_dirs, shogun_db_files = glob_wildcards(config['database'] + "/{dir}/{file}")

uds_runs, sample_names = glob_wildcards("data/hiseq4000/{uds_run}/{sample_name}.fastq.gz")

rule all:
    input:
        #expand("figs/fig{f}.pdf", f=[1,]),
        expand("results/uds/{uds_run}.{sample_name}/alignment.burst.b6", zip, uds_run=uds_runs, sample_name=sample_names)

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

rule shogun_place_on_ramdisk:
    input:
        config['database'] + "/{dir}/{file}"
    priority:
        2
    output:
        temp("/dev/shm/{database_name}".format(database_name=config['database_name']) + "/{dir}/{file}")
    shell:
        "cp {input} {output}"

rule shogun_align:
    input:
        expand("/dev/shm/{database_name}/".format(database_name=config['database_name']) + "{dir}/{file}", zip, dir=shogun_db_dirs, file=shogun_db_files),
        queries = "/dev/shm/uds/{uds_run}.{sample_name}/combined_seqs.fna"
    conda:
        "envs/shogun.yaml"
    output:
        "results/uds/{uds_run}.{sample_name}/alignment.burst.b6",
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_name']),
        output="results/uds/{uds_run}.{sample_name}"
    benchmark:
        "results/benchmarks/{sample_name}.shogun.filter.log"
    priority:
        3
    output:
        temp("/dev/shm/uds/{uds_run}.{sample_name}.txt")
    shell:
        "shogun aligner --input {input.queries} --database {params.database} --aligner burst --function --capitalist --level strain --threads 32 --output {params.output}"

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