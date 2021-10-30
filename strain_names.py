import os
import pandas as pd
import config

strains_df = pd.DataFrame(columns=['strain'])

strain_summary = pd.read_excel(os.path.join(config.tables, "strain_summary.xlsx"))
strain_summary = strain_summary[['Species', 'Strain', 'Refseq assembly accession']]
strain_summary = strain_summary[strain_summary['Species'] == 'Pseudomonas aeruginosa']
x = strain_summary['Refseq assembly accession']
strain_summary['Refseq'] = strain_summary['Refseq assembly accession'].apply(lambda r: r.split('.')[0])

strains = pd.read_csv(os.path.join(config.tables, "strains_list.csv"))
strains['strain_normal'] = strains['strain'].apply(lambda r: r.split('.')[0])
for index, row in strains.iterrows():
    strain_series = strain_summary[strain_summary['Refseq'] == row['strain_normal']]['Strain']
    if strain_series.empty:
        strains_df.loc[len(strains_df)] = [row['strain_normal']]
        continue
    strain_name = strain_series.values[0]
    strains_df.loc[len(strains_df)] = [strain_name]

strains_df.index.name = 'index'
strains_df.to_csv(os.path.join(config.tables, "strains_list_names.csv"))
print("script done")
