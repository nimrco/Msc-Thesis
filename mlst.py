import pandas as pd
import random
import os
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

strains = pd.read_csv(os.path.join(config.tables, "strains_list.csv"))
strains['strain_normal'] = strains['strain'].apply(lambda r: r.split('.')[0])

mlsts = get_mlsts()
types = mlsts['MLST'].unique().tolist()

rand_hex = lambda: random.randint(0, 255)
mlsts_dict = {mlst: '#{:02x}{:02x}{:02x}'.format(rand_hex(), rand_hex(), rand_hex()) for mlst in types}  # mlst: color

strain_color = []
for index, row in mlsts.iterrows():
    strain_index_series = strains[strains['strain_normal'] == row['Refseq']]['index']
    if strain_index_series.empty:
        continue
    strain_index = strain_index_series.values[0]
    color = mlsts_dict[row['MLST']]
    tree_color = '\t'.join([str(strain_index), 'clade', color, 'normal', '4'])
    strain_color.append(tree_color)


with open(os.path.join(config.seq_files, "TREE_COLORS.txt"), "a") as tree_file:
    tree_file.write('\n'.join(strain_color))
