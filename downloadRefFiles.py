from ftplib import FTP
import os
import parseData


def download(file, strain_name):
    ftp = FTP(host[0])
    ftp.login()
    ftp.cwd(strain_name)
    ftp.retrbinary("RETR " + file, open(file, 'wb').write)
    ftp.quit()


strains_list_full_path = []
strains_list = []
ftp = FTP("ftp.ncbi.nlm.nih.gov")
ftp.login()
ftp.cwd("genomes/refseq/bacteria/Pseudomonas_aeruginosa")
with open("assembly_summary.txt", 'w') as parent_dir:
    ftp.retrbinary("RETR assembly_summary.txt", parent_dir.write)
    ftp.quit()
    for line in parent_dir:
        if not line.startswith("#"):
            line_tab = line.split('\t')
            strain = line_tab[19].split(":")
            strains_list_full_path.append(strain[1][2:])

host = strains_list_full_path[0].split("/")
os.chdir('data')

for strain in strains_list_full_path:
    strain = strain.split("/", 1)[1]
    prefix = strain.split("/")[-1]
    strains_list.append(prefix)
    os.mkdir(prefix)
    os.chdir(prefix)
    feature_table = prefix + "_" + "feature_table.txt.gz"
    protein = prefix + "_" + "protein.faa.gz"
    cds = prefix + "_" + "cds_from_genomic.fna.gz"
    download(feature_table, strain)
    download(cds, strain)
    download(protein, strain)
    parseData.merge_files(feature_table, cds, protein)
    os.chdir('..')

with open("strains_list", 'w') as strains_file:
    strains_file.write('\n'.join(strains_list))

