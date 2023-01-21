import os
import pandas as pd
from collections import Counter
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
                line_proccess = line.split(">")[1].split("|")
                member_type = line_proccess[0]
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
    cluster_members_set = set()
    cluster_members_list = list()
    with open(os.path.join(config.seq_files, "cluster{m}_output.clstr".format(m=mode))) as clusters_file:
        for line in clusters_file:
            if line.startswith('>'):
                temp_cluster = line.split()[1]
                if cluster:
                    clusters_dict[cluster] = (cluster_members_set, cluster_members_list)
                    cluster = temp_cluster
                    cluster_members_set = set()
                    cluster_members_list = list()
                else:
                    cluster = temp_cluster
            else:
                strain = line.split('>')[1].split('|')[0]
                cluster_members_set.add(strain)
                cluster_members_list.append(strain)
    clusters_dict[cluster] = (cluster_members_set, cluster_members_list)
    return clusters_dict


def get_pseudo_clusters():
    cluster_size_pseudo = {cluster: val for cluster, val in cluster_size.items() if val[1] == 'pseudo'}
    origin_members_dict = get_clusters_dict('_pseudo')
    origin_members_dict_sets = {cluster: members[0] for (cluster, members) in origin_members_dict.items()}
    origin_members_dict_lists = {cluster: members[1] for (cluster, members) in origin_members_dict.items()}

    clusters_map = {}
    pseudo_clusters = list(cluster_size_pseudo.keys())
    strain_counter_other = Counter()
    strain_counter_singles = Counter()

    singles = {cluster: members for (cluster, members) in origin_members_dict_sets.items() if len(members) == 1}

    with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
        for line in combined_file:
            if line.startswith('>'):
                cluster = line.split()[1]
                if cluster in pseudo_clusters:
                    line = combined_file.readline().split('>')[1].split('|')
                    origin_cluster = line[1]
                    origin_cluster_members = origin_members_dict_lists.get(origin_cluster)
                    if origin_cluster in singles:
                        strain_counter_singles.update(origin_cluster_members)
                    else:
                        strain_counter_other.update(origin_cluster_members)
                    # origin_size = len(origin_members_dict.get(origin_cluster))

                    # origin_size = int(line[2])
                    # clusters_map[cluster] = origin_size
    singles_df = pd.DataFrame.from_dict(strain_counter_singles, orient='index').reset_index()
    singles_df.rename(columns={'index': 'strain', 0: 'count'}, inplace=True)
    other_df = pd.DataFrame.from_dict(strain_counter_other, orient='index').reset_index()
    other_df.rename(columns={'index': 'strain', 0: 'count'}, inplace=True)
    singles_df.to_csv(os.path.join(config.tables, "pseudo_only_single.csv"), index=False)
    other_df.to_csv(os.path.join(config.tables, "pseudo_only_other.csv"), index=False)
    # pseudo_df = pd.DataFrame(clusters_map.items(), columns=['cluster', 'size'])
    # pseudo_df.to_csv(os.path.join(config.tables, "pseudo_clusters_size.csv"), index=False)


def get_gene_clusters():
    cluster_size_gene = {cluster: val for cluster, val in cluster_size.items() if val[1] == 'gene'}
    origin_members_dict = get_clusters_dict()

    clusters_map = {}
    gene_clusters = list(cluster_size_gene.keys())

    with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
        for line in iter(combined_file.readline, ''):
            if line.startswith('>'):
                cluster = line.split()[1]
                if cluster in gene_clusters:
                    line = combined_file.readline()
                    cluster_members = set()
                    while not line.startswith('>'):
                        processed_line = line.split('>')[1].split('|')
                        origin_cluster = processed_line[1]
                        origin_members = origin_members_dict.get(origin_cluster)
                        cluster_members.update(origin_members)
                        # origin_size = int(line.split('>')[1].split('|')[2])
                        prev_line = combined_file.tell()
                        line = combined_file.readline()
                        if line == "":
                            break
                    combined_file.seek(prev_line)
                    clusters_map[cluster] = len(cluster_members)
    gene_df = pd.DataFrame(clusters_map.items(), columns=['cluster', 'size'])
    gene_df.to_csv(os.path.join(config.tables, "gene_origin_clusters_size.csv"), index=False)


def get_mixed_clusters():
    cluster_size_mixed = {cluster: val for cluster, val in cluster_size.items() if val[1] == 'mixed'}
    # with mp.Pool(processes=2) as pool:
    #     results_dict = pool.map(get_clusters_dict, ['', '_pseudo'])
    results_dict = get_clusters_dict('_pseudo')
    results_dict_sets = {cluster: members[0] for (cluster, members) in results_dict.items()}
    results_dict_lists = {cluster: members[1] for (cluster, members) in results_dict.items()}
    strain_counter_other = Counter()
    strain_counter_singles = Counter()
    singles = {cluster: members for (cluster, members) in results_dict_sets.items() if len(members) == 1}

    clusters_map = {}
    mixed_clusters = list(cluster_size_mixed.keys())

    with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
        for line in iter(combined_file.readline, ''):
            if line.startswith('>'):
                cluster = line.split()[1]
                if cluster in mixed_clusters:
                    line = combined_file.readline()
                    # pseudogene_members = set()
                    while not line.startswith('>'):
                        processed_line = line.split('>')[1].split('|')
                        origin_type = processed_line[0]
                        if origin_type == 'PseudoGene':
                            origin_cluster = processed_line[1]
                            origin_cluster_members = results_dict_lists.get(origin_cluster)
                            if origin_cluster in singles:
                                strain_counter_singles.update(origin_cluster_members)
                            else:
                                strain_counter_other.update(origin_cluster_members)
                        prev_line = combined_file.tell()
                        line = combined_file.readline()
                        if line == "":
                            break
                    combined_file.seek(prev_line)
    singles_df = pd.DataFrame.from_dict(strain_counter_singles, orient='index').reset_index()
    singles_df.rename(columns={'index': 'strain', 0: 'count'}, inplace=True)
    other_df = pd.DataFrame.from_dict(strain_counter_other, orient='index').reset_index()
    other_df.rename(columns={'index': 'strain', 0: 'count'}, inplace=True)
    singles_df.to_csv(os.path.join(config.tables, "mixed_single.csv"), index=False)
    other_df.to_csv(os.path.join(config.tables, "mixed_other.csv"), index=False)

    with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
        for line in iter(combined_file.readline, ''):
            if line.startswith('>'):
                cluster = line.split()[1]
                if cluster in mixed_clusters:
                    line = combined_file.readline()
                    # cluster_members = set()
                    gene_members = set()
                    pseudogene_members = set()
                    while not line.startswith('>'):
                        processed_line = line.split('>')[1].split('|')
                        origin_type = processed_line[0]
                        origin_cluster = processed_line[1]
                        if origin_type == 'Gene':
                            origin_members = results_dict[0].get(origin_cluster)
                            gene_members.update(origin_members)
                        else:
                            origin_members = results_dict[1].get(origin_cluster)
                            pseudogene_members.update(origin_members)
                        # cluster_members.update(origin_members)
                        prev_line = combined_file.tell()
                        line = combined_file.readline()
                        if line == "":
                            break
                    combined_file.seek(prev_line)
                    clusters_map[cluster] = (len(gene_members), len(pseudogene_members))
    mixed_df = pd.DataFrame.from_dict(clusters_map,
                                      orient='index',
                                      columns=['gene count', 'pseudogene count'])
    mixed_df.to_csv(os.path.join(config.tables, "mixed_clusters_size.csv"), index=False)

cluster_size = get_cluster_distribution()
get_gene_clusters()
get_pseudo_clusters()
get_mixed_clusters()
print("script done")
# sns.set_theme(style="darkgrid")
# data_table = pd.read_csv("mixed_size_type_repr.csv")
# g = sns.JointGrid(data=data_table,
#                   x="Number of genes",
#                   y="Number of pseudogenes",
#                   hue="representative type"
#                   )
# g.plot_joint(sns.scatterplot, alpha=.7, s=50)
# g.plot_marginals(sns.kdeplot, fill=True)
# # sns.jointplot(data=data_table, x="gene count", y="pseudogene count", marginal_ticks=True, marginal_kws=dict(bins=100))
# # # sns.displot(data=data_table, x='size', bins=50)
# # plt.ylabel('# of pseudogenes')
# # plt.xlabel('# of genes')
# # # plt.title('Mixed representatives clusters sizes')
# # plt.xlim(-50, 5000)
# # # # plt.ylim(0, 800)
# plt.show()
# plt.savefig('mixed_clusters_repr_type.png', dpi=300)
