import os
import sys
import subprocess
import pandas as pd
import config


def preprocess():
    seq_list = []
    strains_df = pd.read_csv("strains_list.csv")
    for strain_index, strain_row in strains_df.iterrows():
        strain = strain_row.values[1]
        pseudo_genes = pd.read_csv(os.path.join(config.data, strain, "pseudo_genes.csv"))
        if not pseudo_genes.empty:
            for seq_index, seq_row in pseudo_genes.iterrows():
                header = str(strain_index) + "|" + str(seq_index)
                seq_dna = seq_row["dna"]
                seq_list.append(">{}\n{}".format(header, seq_dna))

    with open("pseudo_seq.fasta", "w") as pseudo_file:
        pseudo_file.write("\n".join(seq_list))

    print("preprocess done")


def cluster(input_file, output_file):
    cd_hit_args = " ".join(["cd-hit-est", "-i", input_file, "-o", output_file, "-c 0.8", "-n 5", "-M 0",
                            "-g 1", "-p 1", "-d 0", "-T 0"])
    subprocess.run(cd_hit_args, shell=True)
    print("cd-hit done")


def main():
    option = sys.argv[1]
    if option == "preprocess":
        preprocess()
    else:
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        cluster(input_file, output_file)


if __name__ == "__main__":
    main()
