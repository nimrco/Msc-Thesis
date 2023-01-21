import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import config


def get_clusters_dict(mode=''):
    clusters_dict = {}
    cluster = None
    cluster_members = set()
    with open(os.path.join(config.seq_files, "cluster{m}_output.clstr".format(m=mode))) as clusters_file:
        for line in clusters_file:
            if line.startswith('>'):
                temp_cluster = line.split()[1]
                if cluster:
                    clusters_dict[cluster] = cluster_members
                    cluster = temp_cluster
                    cluster_members = set()
                else:
                    cluster = temp_cluster
            else:
                strain = line.split('>')[1].split('|')[0]
                cluster_members.add(strain)
    clusters_dict[cluster] = cluster_members
    return clusters_dict


def get_count_and_overlap():
    pseudo_overlap_dict = defaultdict(set)
    pseudo_count_dict = defaultdict(set)
    pseudo_clusters_dict = get_clusters_dict(mode='_pseudo')
    for cluster, strains in pseudo_clusters_dict.items():
        strains_list = list(strains)
        for i in strains_list:
            pseudo_count_dict[i].add(cluster)
            for j in strains_list:
                pseudo_count_dict[j].add(cluster)
                if i != j:
                    pseudo_overlap_dict[(i, j)].add(cluster)
                    pseudo_overlap_dict[(j, i)].add(cluster)

    pseudo_count_dict_counter = {strain: len(clusters) for strain, clusters in pseudo_count_dict.items()}
    strains_index = len(pseudo_count_dict_counter) + 1
    count_df = pd.DataFrame.from_dict(pseudo_count_dict_counter, orient='index', columns=['clusters count'])
    count_df.to_csv(os.path.join(config.tables, "pseudo_count_test.csv"))
    pseudo_overlap_dict_counter = {strains_pair: len(clusters) for strains_pair, clusters in pseudo_overlap_dict.items()}
    overlap_matrix = np.zeros((strains_index, strains_index), dtype=int)
    for strains_pair, count in pseudo_overlap_dict_counter.items():
        i, j = strains_pair
        i = int(i)
        j = int(j)
        overlap_matrix[i, j] = count
        overlap_matrix[j, i] = count
    pd.DataFrame(overlap_matrix).to_csv(os.path.join(config.tables, "pseudo_overlap_test.csv"))
    print("overlap")


overlap_df = pd.read_csv(os.path.join(config.tables, "pseudo_overlap.csv"))
counts = pd.read_csv(os.path.join(config.tables, "pseudo_count.csv"))
counts.rename(columns={"Unnamed: 0": "strain"}, inplace=True)
tips = pd.read_csv(os.path.join(config.tables, "tips.csv"))
index_list = []
for index, row in tips.iterrows():
    strain = int(row.values[0])
    index_list.append(strain)
overlap_matrix = overlap_df.to_numpy(dtype=int)[:, 1:]
overlap_matrix_jaccard = np.zeros((overlap_matrix.shape[0], overlap_matrix.shape[1]))
for row_index, row in enumerate(overlap_matrix):
    for column_index, column in enumerate(row):
        value_to_insert = overlap_matrix[row_index, column_index]
        if value_to_insert != 0:
            overlap_matrix[column_index, row_index] = value_to_insert
for row_index, row in enumerate(overlap_matrix_jaccard):
    for column_index, column in enumerate(row):
        intersec = overlap_matrix[index_list[row_index], index_list[column_index]]
        i = index_list[row_index]
        j = index_list[column_index]
        if i != 0 and j != 0:
            count_i = counts[counts["strain"] == i]["clusters count"].values[0]
            count_j = counts[counts["strain"] == j]["clusters count"].values[0]
        else:
            count_i = 0
            count_j = 0
        union = count_i + count_j - intersec
        value_to_insert = intersec / union
        overlap_matrix_jaccard[row_index, column_index] = value_to_insert
overlap_df_jaccard = pd.DataFrame(overlap_matrix_jaccard)
overlap_df_jaccard.columns = index_list
overlap_df_jaccard.index = index_list
overlap_df_jaccard.to_csv("pseudo_overlap_jaccard.csv", index=False)



overlap_df = pd.read_csv("pseudo_overlap.csv")
tips = pd.read_csv("tips.csv")
index_list = []
for index, row in tips.iterrows():
    strain = int(row.values[0])
    index_list.append(strain)
overlap_matrix = overlap_df.to_numpy(dtype=int)[:, 1:]
overlap_matrix_sym = np.zeros((overlap_matrix.shape[0], overlap_matrix.shape[1]), dtype=int)
for row_index, row in enumerate(overlap_matrix):
    for column_index, column in enumerate(row):
        value_to_insert = overlap_matrix[row_index, column_index]
        if value_to_insert != 0:
            overlap_matrix[column_index, row_index] = value_to_insert
for row_index, row in enumerate(overlap_matrix_sym):
    for column_index, column in enumerate(row):
        value_to_insert = overlap_matrix[index_list[row_index], index_list[column_index]]
        overlap_matrix_sym[row_index, column_index] = value_to_insert
        # overlap_matrix_sym[column_index, row_index] = value_to_insert
overlap_df_updated = pd.DataFrame(overlap_matrix_sym)
overlap_df_updated.columns = index_list
overlap_df_updated.index = index_list
overlap_df_updated.to_csv("pseudo_overlap_updated.csv", index=False)


get_count_and_overlap()
overlap_df = pd.read_csv("pseudo_overlap.csv")
overlap_matrix = overlap_df.to_numpy(dtype=int)[:, 1:]
tips = pd.read_csv("tips.csv")
overlap_matrix_tree = np.zeros((overlap_matrix.shape[0], overlap_matrix.shape[1]), dtype=int)
index_list = []
for index, row in tips.iterrows():
    strain = int(row.values[0])
    index_list.append(strain)
for row_index, row in enumerate(overlap_matrix_tree):
    for column_index, column in enumerate(row):
        overlap_matrix_tree[row_index, column_index] = overlap_matrix[index_list[row_index], index_list[column_index]]
overlap_df = pd.DataFrame(overlap_matrix_tree)
overlap_df.columns = index_list
overlap_df.index = index_list

overlap_df.to_csv("pseudo_overlap_tree.csv", index=False)

cmap = sns.color_palette("rocket_r", as_cmap=True)
# # plt.subplots(figsize=(30, 20))
# # # strains_index = np.arange(4699)
overlap_df = pd.read_csv("pseudo_overlap_updated.csv")
sns.heatmap(overlap_df, xticklabels=False, yticklabels=False, cmap=cmap)
plt.title("Count of shared pseudogenes between strains")
# # ax.set_xticklabels(ax.get_xticklabels(), fontsize=5, rotation=90)
# plt.savefig("heatmap_reversed_100.png", dpi=300)
# plt.show()

print("script done")
