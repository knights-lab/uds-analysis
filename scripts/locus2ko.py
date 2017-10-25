import pandas as pd
import csv
from collections import defaultdict

with open(snakemake.input.locus2ko) as inf:
    csv_inf = csv.reader(inf, delimiter="\t")
    locus2ko = dict(csv_inf)

i = 0
keggs_per_sample = defaultdict(lambda: defaultdict(int))
with open(snakemake.input.alignment) as inf:
    csv_inf = csv.reader(inf, delimiter="\t")
    for row in csv_inf:
        locus = row[1].split()[0]
        if locus in locus2ko:
            i += 1
            ko = locus2ko[locus].split(';')[1]
            sample = row[0].split("_")[0]
            keggs_per_sample[sample][ko] += 1

kegg_table_df = pd.DataFrame(keggs_per_sample)
kegg_table_df.to_csv(snakemake.output, sep='\t', float_format="%d", na_rep=0, index_label="#KEGG ID")
