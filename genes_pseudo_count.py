import os
import pandas as pd
from matplotlib import pyplot as plt


root = "data"
strains = os.listdir(root)
strain_dict = {}
print(os.getcwd())
for index, strain in enumerate(strains):
    genes = pd.read_csv(os.path.join(root, strain, "genes.csv"))
    pseudo = pd.read_csv(os.path.join(root, strain, "pseudo_genes.csv"))
    genes_num = len(genes)
    pseudo_num = len(pseudo)
    strain_dict[index] = {"genes": genes_num, "pseudo": pseudo_num}
    print(genes)
    print(pseudo)
