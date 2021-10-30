import sys
import os
import pandas as pd
import config


def get_cluster_size(mode):
    cluster_size = {}
    cluster_members = 0
    cluster = None
    with open(os.path.join(config.seq_files, "cluster{}_output.clstr".format(mode))) as cluster_file:
        for line in cluster_file:
            if line.startswith(">"):
                temp_cluster = line.split()[1]
                if cluster:
                    cluster_size[cluster] = cluster_members
                    cluster = temp_cluster
                    cluster_members = 0
                else:
                    cluster = temp_cluster
            else:
                cluster_members += 1
    cluster_size[cluster] = cluster_members
    return cluster_size


def get_repr(mode):
    strains_df = pd.read_csv(os.path.join(config.tables, "strains_list.csv"))
    seq_list = []
    seq_dict = pd.read_csv(os.path.join(config.tables, "seq{}_genes_dict.csv".format(mode))).to_dict()
    cluster_size_dict = get_cluster_size(mode)

    with open(os.path.join(config.seq_files, "cluster{}_output.clstr".format(mode))) as cluster_file:
        for line in cluster_file:
            if line.startswith(">"):
                cluster = line.split()[1]
            if line.endswith("*\n"):
                size = str(cluster_size_dict[cluster])
                length = line.split()[1].split("nt" if mode == '_pseudo' else "a")[0]
                data = line.split(">")[1].split("|")
                strain_index = int(data[0])
                seq_index = int(data[1].split(".")[0])
                strain_name = strains_df.iloc[strain_index]["strain"]
                seq = seq_dict[strain_name][seq_index]

                # header: Gene/PseudoGene|cluster|cluster size|length|strain|seq
                header = "|".join(['{}Gene'.format(mode.strip('_').capitalize()),
                                   cluster,
                                   size,
                                   length,
                                   str(strain_index),
                                   str(seq_index)])
                seq_list.append(">{}\n{}".format(header, seq))

    file_path = 'pseudo_genes_repr.fasta' if mode == '_pseudo' else 'genes_repr.fasta'
    with open(os.path.join(config.seq_files, file_path), "w") as repr_file:
        repr_file.write("\n".join(seq_list))

    print("script done")


def main():
    if len(sys.argv) > 1:
        get_repr(sys.argv[1])
    else:
        get_repr('')


if __name__ == "__main__":
    main()
