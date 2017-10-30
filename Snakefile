__author__ = "Benjamin Hillmann"
__license__ = "AGPL"

from glob import glob
from snakemake.utils import min_version
import multiprocessing
import re
import os

threads_max = multiprocessing.cpu_count()

min_version("3.11.2")

configfile: "config.yaml"

shogun_db_files = glob("{database}/{name}/**/*.*".format(database=config['database'], name=config['database_name']), recursive=True)
shogun_db_files = [file.replace(config['database'] + "/", '') for file in shogun_db_files]
shogun_genes_db_files = glob("{database}/{name}/**/*.*".format(database=config['database'], name=config['database_genes_name']), recursive=True)
shogun_genes_db_files = [file.replace(config['database'] + "/", '') for file in shogun_genes_db_files]
#shogun_db_dirs, shogun_db_files = glob_wildcards(config['database'] + "/{dir}/{file}")
#import ipdb; ipdb.set_trace()

uds_names = ['160729_K00180_0226_AH7WCCBBXX', '160729_K00180_0227_BHCT3LBBXX']
uds_runs, sample_names = glob_wildcards("data/hiseq4000/{uds_run}/{sample_name}.fastq.gz")

# Filter for all reads
#uds_runs, sample_names = zip(*((uds_run, sample_name) for uds_run, sample_name in zip(uds_runs, sample_names) if not sample_name.startswith("Undetermined") and uds_names[0] == uds_run))
uds_runs, sample_names = zip(*((uds_run, sample_name) for uds_run, sample_name in zip(uds_runs, sample_names) if not sample_name.startswith("Undetermined")))

def format_sample_name(sample_name):
    return re.sub(r"_R[1-2]", "", sample_name).replace("_", ".")

rule all:
    input:
        #expand("figs/fig{f}.pdf", f=[1,]),
        #expand("results/hiseq4000/{uds_run}/{sample_name_qc}/taxatable.burst.strain.txt", zip, uds_run=uds_runs, sample_name_qc=set(map(format_sample_name, sample_names)))
        expand("results/hiseq4000/{uds_run}/combined/taxatable.burst.strain.txt", uds_run=set(uds_runs))
        #expand("results/hiseq4000/{uds_run}/{sample_name}.fastq", zip, uds_run=uds_runs, sample_name=sample_names)



rule sample0_genes:
    input:
        expand("results/hiseq4000/{uds_run}/combined/taxatable.genes.kegg.txt", uds_run=uds_names[0])

rule sample1_genes:
    input:
        expand("results/hiseq4000/{uds_run}/combined/taxatable.genes.kegg.txt", uds_run=uds_names[1])

rule sample0_genomes:
    input:
        "results/hiseq4000/{uds_run}/combined/taxatable.burst.strain.kegg.txt".format(uds_run=uds_names[0]),
        "results/hiseq4000/{uds_run}/flat/taxatable.burst.strain.kegg.txt".format(uds_run=uds_names[0])


rule sample1_genomes:
    input:
        "results/hiseq4000/{uds_run}/combined/taxatable.burst.strain.kegg.txt".format(uds_run=uds_names[1]),
        "results/hiseq4000/{uds_run}/flat/taxatable.burst.strain.kegg.txt".format(uds_run=uds_names[1])
        #"results/hiseq4000/{uds_run}/combined/confidence.genomes.txt".format(uds_run=uds_names[1])

def debug(wildcards):
    import ipdb; ipdb.set_trace()
    return ""


############# Genomes ##############

rule extract_uds:
    input:
        # debug
        "data/hiseq4000/{uds_run}/{sample_name}.fastq.gz"
    output:
        "results/hiseq4000/{uds_run}/{sample_name}.fastq"
    shell:
        "7z x {input} -so > {output}"

rule quality_control_uds:
    input:
        lambda wildcards: glob("results/hiseq4000/{uds_run}".format(uds_run=wildcards.uds_run) + "/*.fastq")
    output:
        "results/hiseq4000/{uds_run}/quality_control/combined_seqs.fna",
        "results/hiseq4000/{uds_run}/quality_control/shi7.log"
    params:
        os.path.abspath("results/hiseq4000") + "/{uds_run}"
    priority: 1
    threads:
        int(threads_max/2)
    shell:
        "shi7 --input {params} --output {params} --adaptor Nextera --combine_fasta True --trim_qual 34 --filter_qual 36 --threads {threads} --allow_outies False --flash True"

rule split_quality_control:
    input:
        "results/hiseq4000/{uds_run}/quality_control/combined_seqs.fna"
    output:
        expand("results/hiseq4000/{{uds_run}}/quality_control/x{n}", n=range(10))
    shell:
        "split -a 1 -d -l 200000000 {input}"

rule shogun_place_on_ramdisk:
    input:
        config['database'] + "/{file}"
    priority:
        2
    output:
        temp("/dev/shm/{file}")
    shell:
        "cp {input} {output}"

rule shogun_filter_contaminates:
    input:
        expand("/dev/shm/{file}", file=shogun_db_files),
        queries="results/hiseq4000/{uds_run}/quality_control/x{n}"
    output:
        "results/hiseq4000/{uds_run}/filtered_{n}/combined_seqs.fna",
        "results/hiseq4000/{uds_run}/filtered_{n}/combined_seqs.filtered.fna",
        "results/hiseq4000/{uds_run}/filtered_{n}/alignment.burst.best.b6",
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_name']),
        output="results/hiseq4000/{uds_run}/filtered_{n}",
    benchmark:
        "results/benchmarks/filtered_{n}/shogun.filtered.log"
    priority:
        3
    threads:
        int(threads_max)
    shell:
        "shogun filter --input {input.queries} --database {params.database} --threads {threads} --output {params.output}"  

rule shogun_align:
    input:
        expand("/dev/shm/{file}", file=shogun_db_files),
        qc_log="results/hiseq4000/{uds_run}/shi7.log",
        queries="results/hiseq4000/{uds_run}/filtered_{n}/combined_seqs.fna"
    output:
        "results/hiseq4000/{uds_run}/alignment_{n}/alignment.burst.b6",
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_name']),
        output="results/hiseq4000/{uds_run}/alignment_{n}",
    benchmark:
        "results/benchmarks/alignment_{n}/shogun.align.log"
    priority:
        3
    threads:
        int(threads_max)
    shell:
        "shogun align --input {input.queries} --database {params.database} --aligner burst --threads {threads} --output {params.output}"

rule combine_alignments:
    input:
        lambda wildcards: expand("results/hiseq4000/{uds_run}/alignment_{n}/alignment.burst.b6", n=range(10), uds_run=wildcards.uds_run)
    output:
        "results/hiseq4000/{uds_run}/combined/alignment.burst.b6"
    shell:
        "cat {input} > {output}"

rule shogun_assign_taxonomy:
    input:
        expand("/dev/shm/{file}", file=shogun_db_files),
        alignment = "results/hiseq4000/{uds_run}/combined/alignment.burst.b6"
    output:
        "results/hiseq4000/{uds_run}/combined/taxatable.burst.strain.txt",
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_name'])
    benchmark:
        "results/benchmarks/{uds_run}/shogun.filter.log"
    priority:
        4
    shell:
        "shogun assign_taxonomy --input {input.alignment} --database {params.database} --output {output}"

rule shogun_functions:
    input:
        expand("/dev/shm/{file}", file=shogun_db_files),
        taxatable = "results/hiseq4000/{uds_run}/{style}/taxatable.burst.strain.txt"
    output:
        "results/hiseq4000/{uds_run}/{style}/taxatable.burst.strain.kegg.txt",
        "results/hiseq4000/{uds_run}/{style}/taxatable.burst.strain.normalized.txt",
        "results/hiseq4000/{uds_run}/{style}/taxatable.burst.strain.kegg.modules.txt",
        "results/hiseq4000/{uds_run}/{style}/taxatable.burst.strain.kegg.modules.coverage.txt"
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_name']),
        output="results/hiseq4000/{uds_run}/{style}"
    benchmark:
        "results/benchmarks/{style}/shogun.functional.log"
    shell:
        "shogun functional --input {input.taxatable} --database {params.database} --level strain --output {params.output}"

rule shogun_coverage:
    input:
        expand("/dev/shm/{file}", file=shogun_db_files),
        alignment = "results/hiseq4000/{uds_run}/combined/alignment.burst.b6"
    output:
        "results/hiseq4000/{uds_run}/{sample_name_qc}/confidence.genomes.txt"
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_name']),
    benchmark:
        "results/benchmarks/{sample_name_qc}/shogun.cover.log"
    shell:
        "shogun coverage --input {input.alignment} --database {params.database} --level strain --output {output}"

rule shogun_assign_taxonomy_flat:
    input:
        alignment = "results/hiseq4000/{uds_run}/combined/alignment.burst.b6"
    output:
        "results/hiseq4000/{uds_run}/flat/taxatable.burst.strain.txt",
    benchmark:
        "results/benchmarks/{uds_run}/shogun.flat.log"
    script:
        "scripts/locus2strain.py" 

################# Genes #################
rule shogun_genes_align:
    input:
        expand("/dev/shm/{file}", file=shogun_genes_db_files),
        qc_log="results/hiseq4000/{uds_run}/shi7.log",
        queries="results/hiseq4000/{uds_run}/filtered_{n}/combined_seqs.fna"
    output:
        alignment="results/hiseq4000/{uds_run}/genes_alignment_{n}/alignment.burst.b6",
        logfile="results/hiseq4000/{uds_run}/genes_alignment_{n}/shogun.log.txt"
    params:
        database="/dev/shm/{database_name}".format(database_name=config['database_genes_name']),
        output="results/hiseq4000/{uds_run}/genes_alignment_{n}"
    benchmark:
        "results/benchmarks/genes_alignment_{n}/shogun.align.log"
    priority:
        3
    threads:
        int(threads_max)
    shell:
        "shogun --log debug align --input {input.queries} --database {params.database} --aligner burst --threads {threads} --output {params.output} &> {output.logfile}"

rule combine_genes_alignments:
    input:
        lambda wildcards: expand("results/hiseq4000/{uds_run}/genes_alignment_{n}/alignment.burst.b6", n=range(9), uds_run=wildcards.uds_run)
    output:
        "results/hiseq4000/{uds_run}/combined/genes.alignment.burst.b6"
    shell:
        "cat {input} > {output}"

rule map_alignment2kos:
    input:
        alignment="results/hiseq4000/{uds_run}/combined/genes.alignment.burst.b6"
    output:
        "results/hiseq4000/{uds_run}/combined/taxatable.genes.kegg.txt"
    script:
        "scripts/locus2ko.py"


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