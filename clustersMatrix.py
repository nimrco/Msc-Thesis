import pandas as pd

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
            clusters_dict[cluster].append(strain)

clusters_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in clusters_dict.items()]))
clusters_df.transpose()
print(clusters_df)