import os
import pandas as pd
import config


strains_df = pd.read_csv(os.path.join(config.tables, "strains_list_names.csv"))
strains_df['# genes'] = 0
strains_df['# pseudogenes'] = 0
strains = os.listdir(config.data)
for index, strain in enumerate(strains):
    genes = pd.read_csv(os.path.join(config.data, strain, "genes.csv"))
    pseudo = pd.read_csv(os.path.join(config.data, strain, "pseudo_genes.csv"))
    genes_num = len(genes)
    pseudo_num = len(pseudo)
    strains_df.at[index, '# genes'] = genes_num
    strains_df.at[index, '# pseudogenes'] = pseudo_num

strains_df.to_csv(os.path.join(config.tables, "genes_pseudogenes_count.csv"), index=False)
print("script done")

