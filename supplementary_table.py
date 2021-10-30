import os
import pandas as pd
from mlst import get_mlsts
import config


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
# bioprojects = assembly_summary.groupby('bioproject')
# assembly_summary = bioprojects.filter(lambda x: len(x) >= 5)

prokaryotes.dropna(subset=['RefSeq FTP'], inplace=True)
prokaryotes['RefSeq'] = prokaryotes['RefSeq FTP'].apply(lambda r: r.split('/')[-1].split('.')[0].strip('GCF_'))
bioprojects = prokaryotes.groupby('BioProject')
bioprojects_df = bioprojects.filter(lambda x: len(x) >= 5)

mlsts = get_mlsts()

for index, row in genes_pseudo_count.iterrows():
    refseq = assembly_summary.at[index, 'RefSeq']
    bioproject_series = bioprojects_df[bioprojects_df['RefSeq'] == refseq]['BioProject']
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

print("script done")
