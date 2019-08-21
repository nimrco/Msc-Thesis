from ftplib import FTP
import os

strains_list = []
with open("assembly_summary.txt") as parent_dir:
    for line in parent_dir:
        if not line.startswith("#"):
            line_tab = line.split('\t')
            strain = line_tab[19].split(":")
            strains_list.append(strain[1][2:])

host = strains_list[0].split("/")
ftp = FTP(host[0])
ftp.login()
root = ftp.pwd()

for strain in strains_list:
    strain = strain.split("/", 1)[1]
    ftp.cwd(strain)
    prefix = strain.split("/")[-1]
    os.mkdir(prefix)
    os.chdir(prefix)
    feature_table = prefix + "_" + "feature_table.txt.gz"
    protein = prefix + "_" + "protein.faa.gz"
    cds = prefix + "_" + "cds_from_genomic.fna.gz"
    local_feature_table = open(feature_table, 'wb')
    local_protein = open(protein, 'wb')
    local_cds = open(cds, 'wb')
    ftp.retrbinary('RETR ' + feature_table, local_feature_table.write)
    ftp.retrbinary('RETR ' + protein, local_protein.write)
    ftp.retrbinary('RETR ' + cds, local_cds.write)
    local_feature_table.close()
    local_protein.close()
    local_cds.close()
    os.chdir('..')
    ftp.cwd(root)

ftp.quit()

