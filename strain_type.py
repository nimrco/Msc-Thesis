import os
import pandas as pd
import multiprocessing as mp
from cluster_analysis import get_cluster_distribution
from mlst import get_mlsts
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


def cluster_file_for_tree(cluster):
    strains_df = pd.read_csv(os.path.join(config.tables, "strains_list_names.csv"))
    strains_df['gene'] = 0
    strains_df['pseudo'] = 0
    strains_df['MLST'] = ''
    gene_cluster = clusters_map.get(cluster)[0]
    gene_members = results_dict[0].get(gene_cluster)
    pseudo_cluster = clusters_map.get(cluster)[2]
    pseudo_members = results_dict[1].get(pseudo_cluster)
    strains = list(gene_members.union(pseudo_members))
    for strain in strains:
        if strain in gene_members:
            strains_df.at[int(strain), 'gene'] = 1
        if strain in pseudo_members:
            strains_df.at[int(strain), 'pseudo'] = 1
    return strains_df


with mp.Pool(processes=2) as pool:
    results_dict = pool.map(get_clusters_dict, ['', '_pseudo'])


cluster_size = get_cluster_distribution()
mlsts = get_mlsts()
cluster_size = {cluster: val for cluster, val in cluster_size.items() if val[0] == 2
                and val[1] == 'mixed'}
clusters_map = {}
mixed_clusters = list(cluster_size.keys())

with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
    for line in combined_file:
        if line.startswith('>'):
            cluster = line.split()[1]
            if cluster in mixed_clusters:
                line = combined_file.readline().split('>')[1].split('|')
                if line[0] == 'Gene':
                    origin_gene_cluster = line[1]
                    origin_gene_cluster_size = int(line[2])
                    line = combined_file.readline().split('>')[1].split('|')
                    origin_pseudo_cluster = line[1]
                    origin_pseudo_cluster_size = int(line[2])
                    clusters_map[cluster] = (origin_gene_cluster,
                                             origin_gene_cluster_size,
                                             origin_pseudo_cluster,
                                             origin_pseudo_cluster_size)
                else:
                    origin_pseudo_cluster = line[1]
                    origin_pseudo_cluster_size = int(line[2])
                    line = combined_file.readline().split('>')[1].split('|')
                    origin_gene_cluster = line[1]
                    origin_gene_cluster_size = int(line[2])
                    clusters_map[cluster] = (origin_gene_cluster,
                                             origin_gene_cluster_size,
                                             origin_pseudo_cluster,
                                             origin_pseudo_cluster_size)

clusters_map = {cluster: val for cluster, val in clusters_map.items() if val[1] <= 1500 and val[3] >= 50}
clusters_map = {cluster: val for cluster, val in sorted(clusters_map.items(),
                                                        key=lambda item: item[1][1],
                                                        reverse=True)}

for cluster in list(clusters_map.keys())[:10]:
    strains_df = cluster_file_for_tree(cluster)
    for index, row in strains_df.iterrows():
        mlst_series = mlsts[mlsts['Strain'] == row['strain']]['MLST']
        if mlst_series.empty:
            continue
        strains_df.at[index, 'MLST'] = mlst_series.values[0]
    strains_df.to_csv(os.path.join(config.tables, "cluster_{}.csv".format(cluster)), index=False)

print("script done")
