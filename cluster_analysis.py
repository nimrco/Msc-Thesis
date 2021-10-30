import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import config


def get_cluster_distribution():
    cluster_size = {}
    cluster_members = 0
    pseudo_flag = False
    gene_flag = False
    gene_counter = 0
    pseudo_counter = 0
    cluster = None
    with open(os.path.join(config.seq_files, "repr", "combined_clusters.clstr")) as cluster_file:
        for line in cluster_file:
            if line.startswith(">"):
                temp_cluster = line.split()[1]
                if cluster:
                    if pseudo_flag and gene_flag:
                        cluster_type = "mixed"
                    elif gene_flag:
                        cluster_type = "gene"
                    else:
                        cluster_type = "pseudo"
                    cluster_size[cluster] = (cluster_members, cluster_type, gene_counter, pseudo_counter)
                    cluster = temp_cluster
                    cluster_members = 0
                    pseudo_flag = False
                    gene_flag = False
                    gene_counter = 0
                    pseudo_counter = 0
                else:
                    cluster = temp_cluster
            else:
                cluster_members += 1
                member_type = line.split(">")[1].split("|")[0]
                if member_type == 'Gene':
                    gene_flag = True
                    gene_counter += 1
                else:
                    pseudo_flag = True
                    pseudo_counter += 1
    if pseudo_flag and gene_flag:
        cluster_type = "mixed"
    elif gene_flag:
        cluster_type = "gene"
    else:
        cluster_type = "pseudo"
    cluster_size[cluster] = (cluster_members, cluster_type, gene_counter, pseudo_counter)
    return cluster_size


cluster_dict = get_cluster_distribution()
cluster_df = pd.DataFrame.from_dict(cluster_dict, orient='index',
                                    columns=['cluster_size', 'cluster_type', 'gene_count', 'pseudo_count'])

sns.set_theme(style="darkgrid")
ax = sns.countplot(data=cluster_df, x="cluster_size", hue="cluster_type")
y = np.arange(0, 35000, 2000)
plt.legend(loc='upper right')
ax.set_yticks(y)
ax.set_yscale("log")
plt.title("distribution of genes/pseudogenes depending on cluster size")
plt.savefig("genes_pseudo_cluster_size.png")
plt.show()


