from ftplib import FTP
import os
import parseData


strains_list_full_path = []
strains_list = []
with open("assembly_summary.txt") as parent_dir:
    for line in parent_dir:
        if not line.startswith("#"):
            line_tab = line.split('\t')
            strain = line_tab[19].split(":")
            strains_list_full_path.append(strain[1][2:])

host = strains_list_full_path[0].split("/")
os.chdir('data')

for strain in strains_list_full_path:
    ftp = FTP(host[0])
    ftp.login()
    strain = strain.split("/", 1)[1]
    ftp.cwd(strain)
    prefix = strain.split("/")[-1]
    strains_list.append(prefix)
    os.mkdir(prefix)
    os.chdir(prefix)
    feature_table = prefix + "_" + "feature_table.txt.gz"
    protein = prefix + "_" + "protein.faa.gz"
    cds = prefix + "_" + "cds_from_genomic.fna.gz"
    local_feature_table = open(feature_table, 'wb')
    local_protein = open(protein, 'wb')
    local_cds = open(cds, 'wb')
    ftp.retrbinary('RETR ' + feature_table, local_feature_table.write)
    local_feature_table.close()
    ftp.retrbinary('RETR ' + protein, local_protein.write)
    local_protein.close()
    ftp.retrbinary('RETR ' + cds, local_cds.write)
    local_cds.close()
    ftp.quit()
    parseData.merge_files(feature_table, cds, protein)
    os.chdir('..')

with open("strains_list", 'w') as strains_file:
    strains_file.write('\n'.join(strains_list))



