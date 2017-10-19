__author__ = "Benjamin Hillmann"
__license__ = "AGPL"

from snakemake.utils import min_version

min_version("3.11.2")

configfile: "config.yaml"

rule all:
    input:
        expand("figs/fig{f}.pdf", f=[1,])

############# Collect data ##############

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