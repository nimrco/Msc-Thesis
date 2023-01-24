import pandas as pd
import seaborn as sns
from scipy.stats import variation
from matplotlib import pyplot as plt
from collections import Counter
import config
import os

sns.set_theme(style="white")


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


def plot_pseudo_clusters():
    clusters_df = pd.read_csv("clusters_pseudo_plot.csv")
    clusters_df = clusters_df[[clusters_df.columns[0], "strain_num"]]
    # single = clusters_df[clusters_df["strain_num"] == 1]
    # single.rename_axis("cluster", inplace=True)
    # single.to_csv("singletons_list_pseudo.csv")
    sns.displot(data=clusters_df, x='strain_num', bins=100)
    plt.xlabel("Number of strains")
    plt.ylabel("Number of clusters")
    plt.yticks(fontweight='bold')
    plt.xticks(fontweight='bold')
    plt.xlim(-10, 400)  # 800 bins
    plt.ylim(0, 5000)
    plt.savefig("clusters_pseudo_dist_zoom.png", dpi=300)


def plot_clusters():
    clusters_df = pd.read_csv("clusters_plot.csv")
    clusters_df = clusters_df[[clusters_df.columns[0], "strain_num"]]
    sns.displot(data=clusters_df, x='strain_num', bins=50)
    plt.xlabel("Number of strains")
    plt.ylabel("Number of clusters")
    plt.yticks(fontweight='bold')
    plt.xticks(fontweight='bold')
    # plt.title("Gene clusters distribution by size")
    plt.xlim(-50, 400)  # 800 bins
    plt.ylim(0, 40000)
    plt.savefig("clusters_gene_dist.png", dpi=300)


def plot_coeff():
    clusters_df = pd.read_csv("pseudo_coefficient_of_variation_single.csv")
    sns.displot(data=clusters_df, x='Coefficient of variation')
    plt.ylabel("Number of pseudogene clusters")
    # plt.ylim(0, 2500)
    # plt.xlim(-0.1, 1)
    mean_value = clusters_df['Coefficient of variation'].mean()
    plt.axvline(x=mean_value, color="red", label="Mean = {:.2f}".format(mean_value))
    plt.yticks(fontweight='bold')
    plt.xticks(fontweight='bold')
    plt.legend()
    plt.savefig("pseudo_coeff.png", dpi=300)


def get_table():
    clusters_df = pd.read_csv("clusters_pseudo_plot.csv")
    clusters_df = clusters_df[[clusters_df.columns[0], "strain_num"]]
    counts_series = clusters_df["strain_num"].value_counts()
    pd.DataFrame({'Number of strains': counts_series.index, 'Number of clusters': counts_series.values}).to_csv("pseudogenes_count.csv", index=False)


def get_lengths():
    lengths_dict = {}
    cluster = None
    lengths_list = []
    single = pd.read_csv("singletons_list_pseudo.csv")["cluster"].tolist()
    with open("cluster_pseudo_output.clstr") as cluster_file:
        for line in cluster_file:
            if line.startswith('>'):
                temp_cluster = line.split()[1]
                if cluster:
                    if int(cluster) not in single:
                        lengths_dict[cluster] = variation(lengths_list, ddof=1)
                    cluster = temp_cluster
                    lengths_list = []
                else:
                    cluster = temp_cluster
            else:
                length = int(line.split()[1].split("nt")[0])
                lengths_list.append(length)
    if int(cluster) not in single:
        lengths_dict[cluster] = variation(lengths_list, ddof=1)
    lengths_df = pd.DataFrame(lengths_dict.items(), columns=['cluster', 'Coefficient of variation'])
    lengths_df.to_csv('pseudo_coefficient_of_variation_single.csv', index=False)


def get_repr():
    repr_dict = {}
    cluster_size_mixed = {cluster: val for cluster, val in cluster_size.items() if val[1] == 'mixed'}
    mixed_clusters = list(cluster_size_mixed.keys())
    with open(os.path.join(config.seq_files, 'repr', 'combined_clusters.clstr')) as combined_file:
        for line in combined_file:
            if line.startswith('>'):
                cluster = line.split()[1]
                if cluster in mixed_clusters:
                    repr_flag = False
                    while not repr_flag:
                        line = combined_file.readline()
                        if line.endswith('*\n'):
                            member_type = line.split('>')[1].split('|')[0]
                            gene_count = cluster_size_mixed.get(cluster)[2]
                            pseudo_count = cluster_size_mixed.get(cluster)[3]
                            repr_dict[cluster] = (member_type, gene_count, pseudo_count)
                            repr_flag = True
    repr_df = pd.DataFrame.from_dict(repr_dict,
                                     orient='index',
                                     columns=['representative type', 'gene count', 'pseudogene count'])
    repr_df.to_csv(os.path.join(config.tables, 'repr_type.csv'), index=False)


def merge_mixed_clusters():
    repr_df = pd.read_csv(os.path.join(config.tables, 'repr_type.csv')).reset_index()
    mixed_df = pd.read_csv(os.path.join(config.tables, "mixed_clusters_size.csv")).reset_index()
    merged_df = repr_df.merge(mixed_df, on='index')
    merged_df = merged_df[["representative type", "gene count_y", "pseudogene count_y"]]
    merged_df.rename(columns={"gene count_y": "Number of genes", "pseudogene count_y": "Number of pseudogenes"}, inplace=True)
    merged_df.to_csv(os.path.join(config.tables, "mixed_size_type_repr.csv"), index=False)


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
        mlsts_filter = mlst_types.filter(lambda x: len(x) >= 5)
        # mlsts = mlst_types['MLST'].count()
        # mlsts_dict = mlsts.to_dict()
        # mlsts_df = mlsts.to_frame('count')
        # mlsts_df.reset_index(inplace=True)
        # mlsts_df['MLST'] = mlsts_df['MLST'].astype(float)
        return mlst_data, mlst_types
    return mlst_data


# mlst, mlst_groups = get_mlsts()
# mlst_stats = mlst_groups.agg(
#     number_of_strains=pd.NamedAgg(column='MLST', aggfunc=pd.Series.count),
#     mean_of_strains=pd.NamedAgg(column='MLST', aggfunc=np.mean),
#     std_of_strains=pd.NamedAgg(column='MLST', aggfunc=np.std)
# )


# plot_pseudo_clusters()
# clusters_generate()
# plot_clusters()
# get_table()
# get_lengths()
# print("script done")
# plot_coeff()
# cluster_size = get_cluster_distribution()
# get_repr()
# merge_mixed_clusters()
print("script done")
