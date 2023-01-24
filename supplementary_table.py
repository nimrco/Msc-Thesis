import os
import pandas as pd
import numpy as np
from scipy.stats import variation
import config

window_size = 15

def get_mlsts(groupby_flag=True):
    # mlst_data = pd.read_excel(os.path.join(config.tables, "strain_summary.xlsx"))
    mlst_data = pd.read_excel("strain_summary.xlsx")
    mlst_data = mlst_data[['Species', 'Strain', 'Refseq assembly accession', 'MLST Sequence Type']]
    mlst_data = mlst_data[mlst_data['Species'] == 'Pseudomonas aeruginosa']
    mlst_data = mlst_data.dropna(subset=['MLST Sequence Type'])
    mlst_data['MLST'] = mlst_data['MLST Sequence Type'].apply(lambda r: r.split('|')[-1])
    mlst_data = mlst_data[mlst_data['MLST'] != '-']
    mlst_data['Refseq'] = mlst_data['Refseq assembly accession'].apply(lambda r: r.split('.')[0])
    if groupby_flag:
        mlst_types = mlst_data.groupby('MLST')
        # mlsts_filter = mlst_types.filter(lambda x: len(x) >= 5)
        # mlsts = mlst_types['MLST'].count()
        # mlsts_dict = mlsts.to_dict()
        # mlsts_df = mlsts.to_frame('count')
        # mlsts_df.reset_index(inplace=True)
        # mlsts_df['MLST'] = mlsts_df['MLST'].astype(float)
        return mlst_data, mlst_types
    return mlst_data

supp_df = pd.DataFrame(columns=['Strain',
                                'BioProject',
                                'Level',
                                'Size(Mb)',
                                'Scaffolds',
                                '# of genes',
                                '# of pseudogenes',
                                'ratio genes/pseudogenes',
                                'MLST',
                                'Clinical/Enviromental'])

assembly_summary = pd.read_table(os.path.join(config.tables, 'assembly_summary.txt'), skiprows=0, header=1)
prokaryotes = pd.read_csv(os.path.join(config.tables, 'prokaryotes.csv'))
genes_pseudo_count = pd.read_csv(os.path.join(config.tables, 'genes_pseudogenes_count.csv'))
isolation_df = pd.read_csv(os.path.join(config.tables, 'isolates.csv'))

isolation_df.dropna(subset=['Assembly'], inplace=True)
isolation_df['RefSeq'] = isolation_df['Assembly'].apply(lambda r: r.split('.')[0].strip('GCA_'))

assembly_summary['RefSeq'] = assembly_summary['# assembly_accession'].apply(lambda r: r.split('.')[0].strip('GCF_'))
bioprojects = assembly_summary.groupby('bioproject')
assembly_summary = bioprojects.filter(lambda x: len(x) >= 5)

prokaryotes.dropna(subset=['RefSeq FTP'], inplace=True)
prokaryotes['RefSeq'] = prokaryotes['RefSeq FTP'].apply(lambda r: r.split('/')[-1].split('.')[0].strip('GCF_'))

mlsts = get_mlsts(False)

for index, row in genes_pseudo_count.iterrows():
    refseq = assembly_summary.at[index, 'RefSeq']
    bioproject_series = prokaryotes[prokaryotes['RefSeq'] == refseq]['BioProject']
    level_series = prokaryotes[prokaryotes['RefSeq'] == refseq]['Level']
    size_series = prokaryotes[prokaryotes['RefSeq'] == refseq]['Size(Mb)']
    scaffolds_series = prokaryotes[prokaryotes['RefSeq'] == refseq]['Scaffolds']
    if row['# pseudogenes'] == 0:
        ratio = 0
    else:
        ratio = row['# genes']/row['# pseudogenes']
    mlst_series = mlsts[mlsts['Strain'] == row['strain']]['MLST']
    isolation_type_series = isolation_df[isolation_df['RefSeq'] == refseq]['Isolation type']

    supp_df.at[index, 'Strain'] = row['strain']
    if not bioproject_series.empty:
        supp_df.at[index, 'BioProject'] = bioproject_series.values[0]
    if not level_series.empty:
        supp_df.at[index, 'Level'] = level_series.values[0]
    if not size_series.empty:
        supp_df.at[index, 'Size(Mb)'] = size_series.values[0]
    if not scaffolds_series.empty:
        supp_df.at[index, 'Scaffolds'] = scaffolds_series.values[0]
    supp_df.at[index, '# of genes'] = row['# genes']
    supp_df.at[index, '# of pseudogenes'] = row['# pseudogenes']
    supp_df.at[index, 'ratio genes/pseudogenes'] = ratio
    if not mlst_series.empty:
        supp_df.at[index, 'MLST'] = mlst_series.values[0]
    if not isolation_type_series.empty:
        supp_df.at[index, 'Clinical/Enviromental'] = isolation_type_series.values[0]

supp_df.to_csv(os.path.join(config.tables, "supplementary_table.csv"), index=False)



# For article
mlsts, mlsts_groups = get_mlsts()
sizes = mlsts_groups.size()
filter_size = sizes[sizes >= 5]
supp_df = pd.read_csv("supplementary_table.csv")
stats = pd.Series({'count': supp_df['# of pseudogenes'].count(),
                   'mean': np.mean(supp_df['# of pseudogenes']),
                   'std': np.std(supp_df['# of pseudogenes'], ddof=1),
                   'cv': variation(supp_df['# of pseudogenes'], ddof=1)})
# For article
stats_genes = pd.Series({'count': supp_df['# of genes'].count(),
                   'mean': np.mean(supp_df['# of genes']),
                   'std': np.std(supp_df['# of genes'], ddof=1),
                   'cv': variation(supp_df['# of genes'], ddof=1)})

tips = pd.read_csv("tips.csv")
tree_df = pd.DataFrame(columns=['Strains',
                                '# strains',
                                'Mean',
                                'STD',
                                'CV'])

supp_df['MLST count'] = mlsts_groups['MLST'].transform('count')
mlst_filter = supp_df[supp_df['MLST count'] > 4]
mlst_filter = mlst_filter[['Strain', '# of genes', '# of pseudogenes', 'MLST']]
mlst_filter.rename(columns={"Strain": "strain", "# of genes": "genes", "# of pseudogenes": "pseudogenes"}, inplace=True)
mlst_filter.index.name = 'index'
mlst_filter.to_csv("genes_pseudo_mlst.csv")
mean = np.mean(supp_df['MLST count'])
strains_df = pd.read_csv("strains_list.csv")
pseudo_single = pd.read_csv("pseudo_only_single.csv")
pseudo_other = pd.read_csv("pseudo_only_other.csv")
mixed_single = pd.read_csv("mixed_single.csv")
mixed_other = pd.read_csv("mixed_other.csv")
supp_df_join = supp_df.join(strains_df['strain'])
gcf = supp_df_join.pop('strain')
supp_df_join.insert(0, 'GCF_id', gcf)

random_df = pd.DataFrame(columns=['Mean',
                                  'STD',
                                  'CV'])

for i in range(4684):
    strains = tips.sample(n=15)
    strains = strains['x'].to_list()
    pseudogenes = []
    for strain in strains:
        pseudogenes.append(int(supp_df_join.at[strain, '# of pseudogenes']))
    mean = np.mean(pseudogenes)
    std = np.std(pseudogenes, ddof=1)
    cv = variation(pseudogenes, ddof=1)
    random_df.loc[len(random_df)] = [mean, std, cv]

random_df.to_csv("random_groups_15.csv", index=False)

for index, _ in tips.iterrows():
    last_strain = index + window_size
    if last_strain >= 4698:
        last_strain = 4698
    strains = tips.iloc[index:last_strain]['x'].to_list()
    pseudogenes = []
    for strain in strains:
        pseudogenes.append(int(supp_df_join.at[strain, '# of pseudogenes']))
    mean = np.mean(pseudogenes)
    std = np.std(pseudogenes, ddof=1)
    cv = variation(pseudogenes, ddof=1)
    strains_range = supp_df_join.at[index, 'GCF_id'] + ' - ' + supp_df_join.at[last_strain, 'GCF_id']
    tree_df.loc[len(tree_df)] = [strains_range,
                                 len(pseudogenes),
                                 mean,
                                 std,
                                 cv]
    if last_strain == 4698:
        break
tree_df.to_csv("tree_window_15.csv", index=False)

# For article
for index, _ in tips.iterrows():
    last_strain = index + window_size
    if last_strain >= 4698:
        last_strain = 4698
    strains = tips.iloc[index:last_strain]['x'].to_list()
    genes = []
    for strain in strains:
        genes.append(int(supp_df_join.at[strain, '# of genes']))
    mean = np.mean(genes)
    std = np.std(genes, ddof=1)
    cv = variation(genes, ddof=1)
    strains_range = supp_df_join.at[index, 'GCF_id'] + ' - ' + supp_df_join.at[last_strain, 'GCF_id']
    tree_df.loc[len(tree_df)] = [strains_range,
                                 len(genes),
                                 mean,
                                 std,
                                 cv]
    if last_strain == 4698:
        break
tree_df.to_csv("tree_window_15_genes.csv", index=False)


supp_df_join['# pseudo_only_single'] = 0
supp_df_join['# pseudo_only_other'] = 0
supp_df_join['# mixed_single'] = 0
supp_df_join['# mixed_other'] = 0
for _, row in pseudo_single.iterrows():
    strain = int(row['strain'])
    supp_df_join.at[strain, '# pseudo_only_single'] = row['count']
for _, row in pseudo_other.iterrows():
    strain = int(row['strain'])
    supp_df_join.at[strain, '# pseudo_only_other'] = row['count']
for _, row in mixed_single.iterrows():
    strain = int(row['strain'])
    supp_df_join.at[strain, '# mixed_single'] = row['count']
for _, row in mixed_other.iterrows():
    strain = int(row['strain'])
    supp_df_join.at[strain, '# mixed_other'] = row['count']

mean_test = mlsts_groups.agg(
    mean_of=pd.NamedAgg(column="# of pseudogenes", aggfunc='mean')
)
bio_all = pd.Series({'count': supp_df['# of pseudogenes'].count(),
                     'mean': np.mean(supp_df['# of pseudogenes']),
                     'std': np.std(supp_df['# of pseudogenes'], ddof=1),
                     'cv': variation(supp_df['# of pseudogenes'], ddof=1)})
supp_df['Bio mean'] = bioprojects['# of pseudogenes'].transform(np.mean)
supp_df['Bio std'] = bioprojects['# of pseudogenes'].transform(
    lambda x: np.std(x, ddof=1) if len(x) > 1 else 0
)
supp_df['Bio cv'] = bioprojects['# of pseudogenes'].transform(
    lambda x: variation(x, ddof=1) if len(x) > 1 else 0
)
supp_df_join['MLST mean of pseudogenes'] = mlsts_groups['# of pseudogenes'].transform(np.mean)
mlsts_all_des = pd.Series({'count': mlsts['# of pseudogenes'].count(),
                           'mean': np.mean(mlsts['# of pseudogenes']),
                           'std': np.std(mlsts['# of pseudogenes'], ddof=1),
                           'cv': variation(mlsts['# of pseudogenes'], ddof=1)})
supp_df_join['MLST std of pseudogenes'] = mlsts_groups['# of pseudogenes'].transform(
    lambda x: np.std(x, ddof=1) if len(x) > 1 else 0
)
supp_df_join['MLST CV of pseudogenes'] = mlsts_groups['# of pseudogenes'].transform(
    lambda x: variation(x, ddof=1) if len(x) > 1 else 0
)
#
# # subset_cols = ['BioProject', 'BioProject count', 'Bio mean', 'Bio std', 'Bio cv']
# # bio_df = supp_df[subset_cols]
# # bio_df.drop_duplicates(inplace=True)
# # bio_df.to_csv("bioproject_stats.csv", index=False)
# # mlst_df = supp_df_join[subset_cols]
# # subset_cols = ['MLST mean of pseudogenes', 'MLST std of pseudogenes', 'MLST CV of pseudogenes']
# # supp_df_join.drop(columns=subset_cols, inplace=True)
# # mlst_df.drop_duplicates(inplace=True)
# # mlst_df.dropna(inplace=True)
# # supp_df_join.to_csv("supplementary_table_with_counters.csv", index=False)
# # mlst_df.to_csv("mlst_stats.csv", index=False)


# For article
# supp_df = pd.read_csv("supplementary_table_with_counters.csv")
# pseudo_only_high = supp_df.sort_values(by=['# pseudo_only_single'], ascending=False)
print("script done")
