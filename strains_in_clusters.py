import pandas as pd
from matplotlib import pyplot as plt
from collections import Counter


def clusters_generate():
    clusters_dict = {}
    with open("cluster_output.clstr") as cluster_file:
        for line in cluster_file:
            if line.startswith(">"):  # new cluster
                cluster = line.split()[1]
                clusters_dict.update({cluster: []})
            else:
                strain = line.split(">")[1].split("|")[0]
                clusters_dict[cluster].append(strain)

    clusters_dict_count = {cluster: Counter(strain) for cluster, strain in clusters_dict.items()}
    clusters_df = pd.DataFrame(clusters_dict_count).transpose()
    clusters_df["strain_num"] = clusters_df.count(axis=1)
    clusters_df.to_csv("clusters_plot.csv")


def plot_clusters():
    clusters_df = pd.read_csv("clusters_plot.csv")
    clusters_df = clusters_df[[clusters_df.columns[0], "strain_num"]]
    clusters_df["strain_num"].plot.hist(bins=30)
    plt.xlabel("number of strains")
    plt.ylabel("number of clusters")
    plt.title("clusters distribution by size")
    plt.savefig("clusters_strains.png")
    plt.show()


plot_clusters()
