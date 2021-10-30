import os
import pandas as pd
import config


seq_dict = {}

for strain in os.listdir(config.data):
    seq = pd.read_csv(os.path.join(config.data,
                                   strain,
                                   "pseudo_genes.csv"))["dna"]
    if not seq.empty:
        seq_dict[strain] = seq

seq_df = pd.DataFrame(seq_dict)
seq_df.to_csv(os.path.join(config.tables, "seq_pseudo_genes_dict.csv"), index=False)

print("script done")
