import pandas as pd
import os
import config

seq_list = []
strains = os.listdir(config.data)
strains_df = pd.DataFrame(columns=["strain"])
for index, strain in enumerate(strains):
    strains_df.loc[len(strains_df)] = [strain]
    genes = pd.read_csv(os.path.join(config.data, strain, "genes.csv"))
    for seq_index, row in genes.iterrows():
        header = str(index) + "|" + str(seq_index)
        protein = row["protein"]
        seq_list.append(">{}\n{}".format(header, protein))

strains_df.index.name = 'index'
strains_df.to_csv("strains_list.csv")

with open("protein_seq.fasta", "w") as file:
    file.write("\n".join(seq_list))
print("script done")
