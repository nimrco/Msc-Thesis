import os
import pandas as pd
from matplotlib import pyplot as plt


root = "data"
print(os.getcwd())
for strain in os.listdir(root):
    genes = pd.read_csv(os.path.join(root, strain, "genes,csv"))
    print(genes)
