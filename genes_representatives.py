import os
import pandas as pd
import config

strains_df = pd.read_csv(os.path.join(config.tables, "strains_list.csv"))
seq_list = []
seq_dict = {}

for strain in os.listdir("data"):
    seq_dict[strain] = pd.read_csv(os.path.join("data",
                                                strain,
                                                "genes.csv"))["dna"]

with open(os.path.join("seq_files", "cluster_output.clstr")) as cluster_file:
    for line in cluster_file:
        if line.startswith(">"):
            cluster = line.split()[1]
        if line.endswith("*\n"):
            length = line.split()[1].split("a")[0]
            data = line.split(">")[1].split("|")
            strain_index = int(data[0])
            seq_index = int(data[1].split(".")[0])
            strain_name = strains_df.iloc[strain_index]["strain"]
            seq = seq_dict[strain_name].iloc[seq_index]

            # header: Gene|cluster|length|strain|seq
            header = "|".join(['Gene',
                               cluster,
                               length,
                               str(strain_index),
                               str(seq_index)])
            seq_list.append(">{}\n{}".format(header, seq))

with open(os.path.join("seq_files", "genes_repr.fasta"), "w") as genes_file:
    genes_file.write("\n".join(seq_list))

print("script done")
