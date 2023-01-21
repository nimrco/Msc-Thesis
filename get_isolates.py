import os
import pandas as pd
import multiprocessing as mp
import config

window_size = 5

tips_origin = pd.read_csv("tips.csv")
tips_origin["gene state"] = 0
tips_origin["pseudo state"] = 0


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


def agg_stats(df):
    # df = df.replace({"gene state": 1,
    #                  "gene state": 4,
    #                  "pseudo state": 1,
    #                  "pseudo state": 4}, "With same type")
    # df['gene state'] = df['gene state'].replace([1, 4], "With same type")
    # df['pseudo state'] = df['pseudo state'].replace
    df = df.replace([1, 4], "With same type")
    df = df.replace([2, 5], "With opposite type")
    df = df.replace([3, 6], "Alone")

    gene_counts = df.groupby("gene state").agg(
        gene_count=pd.NamedAgg(column="x", aggfunc="count")
    )
    gene_counts.drop(index=0, inplace=True)
    total_count = gene_counts.sum().values[0]
    gene_counts["gene_ratio"] = gene_counts.apply(lambda x: (x / total_count) * 100)
    gene_counts = gene_counts.reindex(["Alone",
                                       "With opposite type",
                                       "With same type"], fill_value=0)
    pseudo_counts = df.groupby("pseudo state").agg(
        pseudo_count=pd.NamedAgg(column="x", aggfunc="count")
    )
    pseudo_counts.drop(0, inplace=True)
    total_count = pseudo_counts.sum().values[0]
    pseudo_counts["pseudo_ratio"] = pseudo_counts.apply(lambda x: (x / total_count) * 100)
    pseudo_counts = pseudo_counts.reindex(["Alone",
                                           "With opposite type",
                                           "With same type"], fill_value=0)
    joined = pseudo_counts.join(gene_counts)
    # joined.rename_axis('State', inplace=True)
    # joined_melt = joined.melt(ignore_index=False)
    # # joined_melt['variable'] = joined_melt['variable'].apply(lambda v: v.replace(["_count", "_ratio"], ""))
    # joined_melt.rename(columns={"variable": "Type"}, inplace=True)
    # joined = joined_melt.reset_index(level=0)
    return joined


def reshape_stats(df):
    df = df.reset_index()
    df = df.pivot(index='level_0', columns=['State', 'Type'], values='value')
    df.rename_axis('Cluster', inplace=True)
    return df

folder = os.path.join(config.tables, "clusters")
states_folder = os.path.join(config.tables, "states")
clusters_dfs = {}


def stats_for_graph():
    combined_stats = pd.DataFrame(columns=['Cluster',
                                           'State',
                                           'Gene Count',
                                           'Pseudo Count',
                                           'Gene Ratio',
                                           'Pseudo Ratio']
                                  )
    # combined_stats.drop(['gene state', 'pseudo state'], axis=1, inplace=True)
    # combined_stats.rename(columns={"x": "Cluster"}, inplace=True)
    # combined_stats['Type'] = 0
    # combined_stats['State'] = 0
    # combined_stats['Count'] = 0
    # combined_stats['Ratio'] = 0
    for cluster_stats_folder in os.listdir(states_folder):
        for file in os.listdir(os.path.join(states_folder, cluster_stats_folder)):
            if file.startswith("cluster"):
                cluster = pd.read_csv(os.path.join(states_folder, cluster_stats_folder, file))
                joined = agg_stats(cluster)
                # joined_list = [list(value[1].values) for value in joined.iterrows()]
                # joined_list = [value for sublist in joined_list for value in sublist]
                for row in joined.itertuples():
                    combined_stats.loc[len(combined_stats)] = [int(cluster_stats_folder),
                                                               row.Index,
                                                               row.gene_count,
                                                               row.pseudo_count,
                                                               row.gene_ratio,
                                                               row.pseudo_ratio]
    combined_stats.to_csv("stats_for_graph.csv", index=False)


def create_stats():
    for file in os.listdir(folder):
        cluster_df = pd.read_csv(os.path.join(folder, file))
        tips = tips_origin.copy()
        for index, row in tips.iterrows():
            if cluster_df.iloc[row["x"]]["pseudo"]:  # pseudo
                tips.iloc[index]["pseudo state"] = 6
                for i in range(1, window_size + 1):
                    forward = index + i
                    if forward < len(tips):
                        strain = tips.iloc[forward]["x"]
                        if cluster_df.iloc[strain]["gene"]:
                            tips.iloc[index]["pseudo state"] = 5
                        if cluster_df.iloc[strain]["pseudo"]:
                            tips.iloc[index]["pseudo state"] = 4
                    reverse = index - i
                    if reverse >= 0:
                        strain = tips.iloc[reverse]["x"]
                        if cluster_df.iloc[strain]["gene"]:
                            tips.iloc[index]["pseudo state"] = 5
                        if cluster_df.iloc[strain]["pseudo"]:
                            tips.iloc[index]["pseudo state"] = 4
            if cluster_df.iloc[row["x"]]["gene"]:  # gene
                tips.iloc[index]["gene state"] = 3
                for i in range(1, window_size + 1):
                    forward = index + i
                    if forward < len(tips):
                        strain = tips.iloc[forward]["x"]
                        if cluster_df.iloc[strain]["pseudo"]:
                            tips.iloc[index]["gene state"] = 2
                        if cluster_df.iloc[strain]["gene"]:
                            tips.iloc[index]["gene state"] = 1
                    reverse = index - i
                    if reverse >= 0:
                        strain = tips.iloc[reverse]["x"]
                        if cluster_df.iloc[strain]["pseudo"]:
                            tips.iloc[index]["gene state"] = 2
                        if cluster_df.iloc[strain]["gene"]:
                            tips.iloc[index]["gene state"] = 1
        cluster = file.split(".csv")[0].split("_")[1]
        cluster_stats_folder = os.path.join(states_folder, cluster)
        os.mkdir(cluster_stats_folder)
        tips.to_csv(os.path.join(cluster_stats_folder, file), index=False)
        stats_df = agg_stats(tips)
        stats_df.to_csv(os.path.join(cluster_stats_folder, "stats.csv"))
        clusters_dfs[cluster] = stats_df


cluster_size = get_cluster_distribution()
cluster_size = {cluster: val for cluster, val in cluster_size.items() if val[0] == 2
                and val[1] == 'mixed'}
clusters_map = {}
mixed_clusters = list(cluster_size.keys())


def get_clusters_dict(mode=''):
    clusters_dict = {}
    cluster = None
    cluster_members = set()
    members_counter = 0
    with open(os.path.join(config.seq_files, "cluster{m}_output.clstr".format(m=mode))) as clusters_file:
        for line in clusters_file:
            if line.startswith('>'):
                temp_cluster = line.split()[1]
                if cluster:
                    ratio = members_counter / len(cluster_members)
                    clusters_dict[cluster] = (cluster_members, members_counter, ratio)
                    cluster = temp_cluster
                    cluster_members = set()
                    members_counter = 0
                else:
                    cluster = temp_cluster
            else:
                strain = line.split('>')[1].split('|')[0]
                cluster_members.add(strain)
                members_counter += 1
    ratio = members_counter / len(cluster_members)
    clusters_dict[cluster] = (cluster_members, members_counter, ratio)
    return clusters_dict


def get_mixed_clusters():
    cluster_size = get_cluster_distribution()
    cluster_size_mixed = {cluster: val for cluster, val in cluster_size.items() if val[1] == 'mixed'}
    with mp.Pool(processes=2) as pool:
        results_dict = pool.map(get_clusters_dict, ['', '_pseudo'])

    clusters_map = {}
    mixed_clusters = list(cluster_size_mixed.keys())

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
                            origin_members = results_dict[0].get(origin_cluster)[0]
                            gene_members.update(origin_members)
                        else:
                            origin_members = results_dict[1].get(origin_cluster)[0]
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
    mixed_df.index.name = 'Cluster'
    mixed_df.to_csv(os.path.join(config.tables, "mixed_clusters_size_index.csv"))


def get_mixed_size():
    cluster_size = get_cluster_distribution()
    cluster_size_mixed = {cluster: val for cluster, val in cluster_size.items() if val[1] == 'mixed' and val[0] == 2}
    with mp.Pool(processes=2) as pool:
        results_dict = pool.map(get_clusters_dict, ['', '_pseudo'])

    clusters_map = {}
    mixed_clusters = list(cluster_size_mixed.keys())

    with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
        for line in iter(combined_file.readline, ''):
            if line.startswith('>'):
                cluster = line.split()[1]
                if cluster in mixed_clusters:
                    line = combined_file.readline()
                    gene_members = set()
                    pseudogene_members = set()
                    while not line.startswith('>'):
                        processed_line = line.split('>')[1].split('|')
                        origin_type = processed_line[0]
                        origin_cluster = processed_line[1]
                        if origin_type == 'Gene':
                            gene_size = processed_line[2]
                            gene_ratio = results_dict[0].get(origin_cluster)[2]
                            origin_members = results_dict[0].get(origin_cluster)[0]
                            gene_members.update(origin_members)
                        else:
                            pseudo_size = processed_line[2]
                            pseudo_ratio = results_dict[1].get(origin_cluster)[2]
                            origin_members = results_dict[1].get(origin_cluster)[0]
                            pseudogene_members.update(origin_members)
                        prev_line = combined_file.tell()
                        line = combined_file.readline()
                        if line == "":
                            break
                    combined_file.seek(prev_line)
                    clusters_map[cluster] = (len(gene_members), len(pseudogene_members), gene_size, gene_ratio, pseudo_size, pseudo_ratio)
    mixed_df = pd.DataFrame.from_dict(clusters_map,
                                      orient='index',
                                      columns=['gene count', 'pseudogene count',
                                               'gene cluster size (seq)', 'gene cluster ratio',
                                               'pseudo cluster size (seq)', 'pseudo cluster ratio'])
    mixed_df.index.name = 'Cluster'
    mixed_df.to_csv(os.path.join(config.tables, "mixed_clusters_ratio.csv"))


get_mixed_size()
get_mixed_clusters()
stats_for_graph()

mixed_clusters = pd.read_csv(os.path.join(config.tables, "mixed_clusters_ratio.csv"))
stats_df = pd.read_csv("stats_for_graph.csv")
gene_origin_clusters = get_clusters_dict()
stats_df['Gene Cluster Size'] = 0
stats_df['Gene Cluster Sequences'] = 0
stats_df['Gene Cluster Ratio'] = 0.0
stats_df['Pseudo Cluster Size'] = 0
stats_df['Pseudo Cluster Sequences'] = 0
stats_df['Pseudo Cluster Ratio'] = 0.0
for index, row in stats_df.iterrows():
    cluster = row['Cluster']
    gene_size = mixed_clusters[mixed_clusters['Cluster'] == cluster]['gene count'].values[0]
    gene_sequences = mixed_clusters[mixed_clusters['Cluster'] == cluster]['gene cluster size (seq)'].values[0]
    gene_ratio = mixed_clusters[mixed_clusters['Cluster'] == cluster]['gene cluster ratio'].values[0]
    pseudo_size = mixed_clusters[mixed_clusters['Cluster'] == cluster]['pseudogene count'].values[0]
    pseudo_sequences = mixed_clusters[mixed_clusters['Cluster'] == cluster]['pseudo cluster size (seq)'].values[0]
    pseudo_ratio = mixed_clusters[mixed_clusters['Cluster'] == cluster]['pseudo cluster ratio'].values[0]
    stats_df.at[index, 'Gene Cluster Size'] = gene_size
    stats_df.at[index, 'Gene Cluster Sequences'] = gene_sequences
    stats_df.at[index, 'Gene Cluster Ratio'] = gene_ratio
    stats_df.at[index, 'Pseudo Cluster Size'] = pseudo_size
    stats_df.at[index, 'Pseudo Cluster Sequences'] = pseudo_sequences
    stats_df.at[index, 'Pseudo Cluster Ratio'] = pseudo_ratio
#
stats_df.to_csv("stats_with_ratio.csv", index=False)




# mixed_clusters = pd.read_csv(os.path.join(config.tables, "mixed_clusters_ratio.csv"))
# stats_df = pd.read_csv("stats_for_graph.csv")
# gene_origin_clusters = get_clusters_dict()
# stats_df['Gene Cluster Size'] = 0
# stats_df['Gene Cluster Sequences'] = 0
# stats_df['Gene Cluster Ratio'] = 0
# stats_df['Pseudo Cluster Size'] = 0
# for index, row in stats_df.iterrows():
#     cluster = row['Cluster']
#     gene_size = mixed_clusters[mixed_clusters['Cluster'] == cluster]['gene count'].values[0]
#     # origin_sequences =
#     pseudo_size = mixed_clusters[mixed_clusters['Cluster'] == cluster]['pseudogene count'].values[0]
#     stats_df.at[index, 'Gene Cluster Size'] = gene_size
#     stats_df.at[index, 'Pseudo Cluster Size'] = pseudo_size
# #
# stats_df.to_csv("stats_with_size.csv", index=False)
# stats_for_graph()
# result = pd.concat(clusters_dfs, keys=clusters_dfs.keys())
# result = reshape_stats(result)
# result.to_csv("mixed_stats.csv")

# stats = pd.read_csv("stats_sorted.csv", header=[0, 1, 2])
# stats_drop = stats.drop("NA", axis=1, level=0)
# stats_drop.columns = stats_drop.columns.remove_unused_levels()
# stats_drop.columns = stats_drop.columns.set_names(['Cluster', 'State', 'Type'])
# states = ['Alone', 'With opposite type', 'With same type']
# types = ['gene_count', 'pseudo_count']
# for index, row in stats_drop.iterrows():
#     for type_count in types:
#         for state in states:
#             stats_drop.loc[index, (state, type_count.replace("count", "ratio"))] = 100 * (row[state][type_count].values[0] / sum([row[state_val][type_count].values[0] for state_val in states]))

# stats_drop.to_csv("stats_without_NA.csv", index=False)

# stats_names = stats_drop['Alone']
# state = ["State",
#          "Alone",
#          "With opposite type",
#          "With same type",
#          "NA"]
# types = ["Type",
#          "gene_count",
#          "pseudo_count",
#          "gene_ratio",
#          "pseudo_ratio"]
# stats = stats.reindex(state, axis=1, level=0)
# stats = stats.reindex(types, axis=1, level=1)
# # stats = stats.reindex("Cluster", axis=1, level=2)
# # stats.to_csv("stats_sorted.csv")

# For article
# with mp.Pool(processes=2) as pool:
#     results_dict = pool.map(get_clusters_dict, ['', '_pseudo'])
print("script done")



