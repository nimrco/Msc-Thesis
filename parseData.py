import pandas as pd
import gzip
from Bio import SeqIO
import re


def extract_file(file):
    uncompressed = open(file[:-3], "wb")
    with gzip.open(file, "rb") as compressed:
        data = compressed.read()
    uncompressed.write(data)
    uncompressed.close()


def merge_files(features_file, cds_file_zip, protein_file_zip):
    print(features_file)
    with gzip.open(features_file, "rt") as features:
        features_data = pd.read_csv(features, sep="\t")

    # genes
    gene = features_data[((features_data["# feature"] == 'gene') & (features_data["class"] == 'protein_coding'))]
    cds = features_data[((features_data["# feature"] == 'CDS') & (features_data["class"] == 'with_protein'))]
    gene_join = gene.merge(cds, on='locus_tag')

    # pseudogenes
    pseudo_gene = features_data[((features_data["# feature"] == 'gene') & (features_data["class"] == 'pseudogene'))]
    pseudo_cds = features_data[((features_data["# feature"] == 'CDS') & (features_data["class"] == 'without_protein'))]
    pseudo_join = pseudo_gene.merge(pseudo_cds, on='locus_tag')

    # proteins
    protein_df = pd.DataFrame(columns=['id', 'protein'])
    extract_file(protein_file_zip)
    for record in SeqIO.parse(protein_file_zip[:-3], "fasta"):
        protein_df.loc[len(protein_df)] = [record.id, record.seq]

    # cds
    dna_df = pd.DataFrame(columns=['locus_tag', 'dna'])
    extract_file(cds_file_zip)
    for record in SeqIO.parse(cds_file_zip[:-3], "fasta"):
        try:
            locus_tag = re.findall('\[locus_tag=(.+?)\]', record.description)[0]
        except AttributeError:
            locus_tag = ''
        dna_df.loc[len(dna_df)] = [locus_tag, record.seq]

    # create genes merged file
    gene_join = gene_join.merge(dna_df, on='locus_tag')
    gene_join = gene_join.merge(protein_df, left_on='product_accession_y', right_on='id')
    gene_join.to_csv("genes.csv")

    # create pseudogenes merged file
    pseudo_join = pseudo_join.merge(dna_df, on='locus_tag')
    pseudo_join.to_csv("pseudo_genes.csv")
