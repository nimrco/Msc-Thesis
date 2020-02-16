import pandas as pd


clusters_dict = {}
with open("cluster_output.clstr") as cluster_file:
    for line in cluster_file:
        if line.startswith(">"):
            cluster = line.split()[1]
            clusters_dict.update({cluster: {}})
        else:
            data = line.split(">")[1].split("|")
            strain = data[0]
            seq = data[1].split(".")[0]
            if strain not in clusters_dict[cluster]:
                clusters_dict[cluster].update({strain: [seq]})
            else:
                clusters_dict[cluster][strain].append(seq)

clusters_df = pd.DataFrame({cluster: pd.Series(strains) for cluster, strains in clusters_dict.items()}).transpose()
clusters_df.to_csv("cluster.csv")
