import os
import pandas as pd
import multiprocessing as mp
import config


def get_mlsts(groupby_flag=True):
    mlst_data = pd.read_excel(os.path.join(config.tables, "strain_summary.xlsx"))
    mlst_data = mlst_data[['Species', 'Strain', 'Refseq assembly accession', 'MLST Sequence Type']]
    mlst_data = mlst_data[mlst_data['Species'] == 'Pseudomonas aeruginosa']
    mlst_data = mlst_data.dropna(subset=['MLST Sequence Type'])
    mlst_data['MLST'] = mlst_data['MLST Sequence Type'].apply(lambda r: r.split('|')[-1])
    mlst_data = mlst_data[mlst_data['MLST'] != '-']
    mlst_data['Refseq'] = mlst_data['Refseq assembly accession'].apply(lambda r: r.split('.')[0])
    if groupby_flag:
        mlst_types = mlst_data.groupby('MLST')
        mlsts = mlst_types.filter(lambda x: len(x) >= 5)
        return mlsts
    return mlst_data


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


def cluster_file_for_tree(cluster, ratio_flag=False):
    gene_cluster = clusters_map.get(cluster)[0]
    gene_members = results_dict[0].get(gene_cluster)
    pseudo_cluster = clusters_map.get(cluster)[2]
    pseudo_members = results_dict[1].get(pseudo_cluster)
    strains = list(gene_members.union(pseudo_members))
    strains_df = pd.read_csv(os.path.join(config.tables, "strains_list_names.csv"))
    strains_df['MLST'] = ''
    if ratio_flag:
        ratio_df = pd.read_csv(os.path.join(config.tables, "supplementary_table.csv"))
        strains_df['ratio'] = 0
        for strain in strains:
            strains_df.at[int(strain), 'ratio'] = ratio_df.at[int(strain), "ratio genes/pseudogenes"]
    else:
        count_df = pd.read_csv(os.path.join(config.tables, "genes_pseudogenes_count.csv"))
        strains_df['gene'] = 0
        strains_df['pseudo'] = 0
        for strain in strains:
            if strain in gene_members:
                strains_df.at[int(strain), 'gene'] = count_df.at[int(strain), "# genes"]
            if strain in pseudo_members:
                strains_df.at[int(strain), 'pseudo'] = count_df.at[int(strain), "# pseudogenes"]
    return strains_df


with mp.Pool(processes=2) as pool:
    results_dict = pool.map(get_clusters_dict, ['', '_pseudo'])

cluster_size = get_cluster_distribution()
mlsts = get_mlsts()
# mixed with one gene and one pseudogene
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

clusters_map = {cluster: val for cluster, val in clusters_map.items() if val[1] <= 2500 and val[3] >= 2}
clusters_map = {cluster: val for cluster, val in sorted(clusters_map.items(),
                                                        key=lambda item: item[1][1],
                                                        reverse=True)}

for cluster in list(clusters_map.keys()):
    strains_df = cluster_file_for_tree(cluster)
    for index, row in strains_df.iterrows():
        mlst_series = mlsts[mlsts['Strain'] == row['strain']]['MLST']
        if mlst_series.empty:
            continue
        strains_df.at[index, 'MLST'] = mlst_series.values[0]
    strains_df.to_csv(os.path.join(config.tables, "clusters", "cluster_{}.csv".format(cluster)), index=False)

print("script done")
