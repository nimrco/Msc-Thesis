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
cluster_df.to_csv("cluster_distribution.csv", index=False)


def change_type(type):
    if type == 'gene':
        return 'Gene-only'
    elif type == 'pseudo':
        return 'Pseudogene-only'
    return 'Mixed'


sns.set_theme(style="white")
cluster_df = pd.read_csv("cluster_distribution.csv")
cluster_df['cluster_type'] = cluster_df['cluster_type'].apply(change_type)
d = sns.color_palette(n_colors=3)
g = sns.countplot(data=cluster_df, x="cluster_size", hue="cluster_type", hue_order=['Mixed', 'Gene-only', 'Pseudogene-only'],
              palette={'Mixed': d[0], 'Gene-only': d[2], 'Pseudogene-only': d[1]})
plt.legend(loc='upper right')
y = np.arange(0, 35000, 2000)
plt.yticks(y)
plt.yscale("log")
plt.xticks(ha='right')
plt.xlabel("Cluster size")
plt.ylabel("Number of clusters")
plt.yticks(fontweight='bold')
plt.xticks(fontweight='bold')

plt.savefig("genes_pseudo_cluster_size.png", dpi=300)


