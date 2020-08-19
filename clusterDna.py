import pandas as pd
import os

root = "data"
core_df = pd.read_csv("core_cluster_matrix.csv")
core_df = core_df[core_df["core_cluster"] == 1].iloc[:, 1:-4]
strains_df = pd.read_csv("strains_list.csv")
clusters_df = pd.read_csv("cluster.csv").iloc[:, 1:]
seq_dict = {}

for strain in os.listdir(root):
    seq_dict[strain] = pd.read_csv(os.path.join(root, strain, "genes.csv"))["dna"]

for cluster, row in core_df.iterrows():
    seq_list = []
    for strain, count in row.iteritems():
        if count == 1:
            strain_name = strains_df.iloc[int(strain)]["strain"]
            seq = clusters_df.iloc[cluster][strain]
            seq = seq.strip(" '\'[]")
            dna = seq_dict[strain_name].iloc[int(seq)]
            # dna = pd.read_csv(os.path.join(root, strain_name, "genes.csv")).iloc[int(seq)]["dna"]
            header = strain + "|" + seq
            seq_list.append(">{}\n{}".format(header, dna))

    os.mkdir(os.path.join("clusters", cluster))
    with open(os.path.join("clusters", cluster, str(cluster) + ".fasta"), "w") as cluster_file:
        cluster_file.write("\n".join(seq_list))

    print("cluster: {} done".format(cluster))

print("script done")
