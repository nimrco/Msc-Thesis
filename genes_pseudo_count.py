import os
import pandas as pd
from matplotlib import pyplot as plt


root = "data"
strains = os.listdir(root)
genes_list = []
pseudo_list = []
strain_dict = {}
index_list = []
for index, strain in enumerate(strains):
    genes = pd.read_csv(os.path.join(root, strain, "genes.csv"))
    pseudo = pd.read_csv(os.path.join(root, strain, "pseudo_genes.csv"))
    genes_num = len(genes)
    pseudo_num = len(pseudo)
    strain_dict[index] = {"genes": genes_num, "pseudo": pseudo_num}
    genes_list.append(genes_num)
    pseudo_list.append(pseudo_num)
    index_list.append(index)

plt.plot(index_list, genes_list, color="blue", label="genes")
plt.plot(index_list, pseudo_list, color="red", label="pseudo")
plt.xlabel("strain number")
plt.ylabel("number of genes/pseudogenes")
plt.title("genes vs psuedogenes")
plt.legend()
plt.savefig("genes_vs_pseudo.png")
plt.show()
