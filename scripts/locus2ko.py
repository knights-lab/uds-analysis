import pandas as pd
import csv
from collections import defaultdict

keggs_per_sample = defaultdict(lambda: defaultdict(int))
with open(snakemake.input.alignment) as inf:
    csv_inf = csv.reader(inf, delimiter="\t")
    for row in csv_inf:
        ko = row[-1]
        if ko:
            ko = ko.split(';')
            if len(ko) == 2:
                sample = row[0].split("_")[0]
                keggs_per_sample[sample][ko[1]] += 1

kegg_table_df = pd.DataFrame(keggs_per_sample)
kegg_table_df.to_csv(snakemake.output[0], sep='\t', float_format="%d", na_rep=0, index_label="#KEGG ID")
