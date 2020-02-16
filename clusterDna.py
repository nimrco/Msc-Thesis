import pandas as pd
import os

root = "data"
core_df = pd.read_csv("core_cluster_matrix.csv")
core_df = core_df[core_df["core_cluster"] == 1].iloc[:, 1:]
strains_df = pd.read_csv("strains_list.csv")
clusters_df = pd.read_csv("cluster.csv").iloc[:, 1:]

for cluster, row in core_df.iterrows():
    seq_list = []
    for strain, count in row.iteritems():
        if count == 1:
            strain_name = strains_df.iloc[int(strain)]["strain"]
            seq = clusters_df.iloc[cluster][strain]
            seq = seq.strip(" '\'[]")
            dna = pd.read_csv(os.path.join(root, strain_name, "genes.csv")).iloc[int(seq)]["dna"]
            header = strain + "|" + seq
            seq_list.append(">{}\n{}".format(header, dna))

    with open(os.path.join("clusters", str(cluster)), "w") as cluster_file:
        cluster_file.write("\n".join(seq_list))

    print("cluster: {} done".format(cluster))
