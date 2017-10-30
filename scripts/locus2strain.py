import pandas as pd
import csv
from collections import defaultdict

strains_per_sample = defaultdict(lambda: defaultdict(int))
with open(snakemake.input.alignment) as inf:
    csv_inf = csv.reader(inf, delimiter="\t")
    for row in csv_inf:
        strain = row[-1]
        if strain:
            if 't__' in strain:
                sample = row[0].split("_")[0]
                strains_per_sample[sample][strain] += 1

strain_table_df = pd.DataFrame(strains_per_sample)
strain_table_df.to_csv(snakemake.output[0], sep='\t', float_format="%d", na_rep=0, index_label="#OTU ID")
