import pandas as pd
from collections import Counter

strains_list = pd.read_csv("strains_list.csv")
num_strains = len(strains_list.index)
core_threshold = (num_strains * 90) / 100

clusters_dict = {}
with open("cluster_output.clstr") as cluster_file:
    for line in cluster_file:
        if line.startswith(">"):  # new cluster
            cluster = line.split()[1]
            clusters_dict.update({cluster: []})
        else:
            strain = line.split(">")[1].split("|")[0]
            strain = int(strain)
            clusters_dict[cluster].append(strain)

clusters_dict_count = {cluster: Counter(strain) for cluster, strain in clusters_dict.items()}
clusters_df_strains = pd.DataFrame(clusters_dict_count).transpose()
clusters_df = clusters_df_strains.copy(deep=True)
clusters_df["zero_copy"] = clusters_df_strains.isnull().sum(axis=1)
clusters_df = clusters_df.fillna(0)
one_copy = (clusters_df_strains == 1).sum(axis=1)
clusters_df["one_copy"] = one_copy
many_copies = (clusters_df_strains > 1).sum(axis=1)
clusters_df["many_copies"] = many_copies
clusters_df["core_cluster"] = clusters_df["one_copy"] >= core_threshold
clusters_df = clusters_df.append(pd.Series(name="num_of_core_clusters"))
clusters_df = clusters_df.append(pd.Series(name="single_copy_core_clusters"))

for strain in clusters_df_strains:
    strain_df = clusters_df.loc[:, [strain, "core_cluster"]]
    core_df = strain_df[(strain_df[strain] > 0) & (strain_df["core_cluster"] > 0)]
    clusters_df.loc["num_of_core_clusters", strain] = core_df[strain].count()

    single_copy_df = strain_df[(strain_df[strain] == 1) & (strain_df["core_cluster"] > 0)]
    clusters_df.loc["single_copy_core_clusters", strain] = single_copy_df[strain].count()

clusters_df.to_csv("core_cluster_matrix.csv")
print("script done")
