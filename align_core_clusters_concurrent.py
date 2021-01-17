import subprocess
import os
import sys
import pandas as pd
import multiprocessing as mp
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import config


def get_alignment_dict():
    alignment_dict_local = {}
    for cluster in os.listdir(config.clusters):
        with open(os.path.join(config.clusters, cluster, cluster + "_alignment-gb")) as alignment_file:
            alignment_dict_local[cluster] = AlignIO.read(alignment_file, "fasta")
    print("alignment dict")
    return alignment_dict_local


def get_first_char(alignment):
    for index in range(len(alignment)):
        c = alignment[index]
        if c != '-':
            return c, index
    return alignment[0], -1


def check_variable_col(c, row_index, alignment):
    for char_row in alignment[row_index+1:]:
        if char_row != c and char_row != '-':
            return True
    return False


def filter_alignment(cluster):
    alignment_lengths_df = pd.DataFrame(columns=["cluster", "original", "edited"])
    lengths = [int(cluster)]
    alignment = alignment_dict[cluster]
    original_length = alignment.get_alignment_length()
    lengths.append(original_length)
    edited_alignment = None
    for col_idx in range(original_length):
        col = alignment[:, col_idx:col_idx + 1]
        col_str = alignment[:, col_idx]
        c, row_index = get_first_char(col_str)
        if row_index == -1:
            continue
        if check_variable_col(c, row_index, col_str):
            if not edited_alignment:
                edited_alignment = col
            else:
                edited_alignment += col
    if not edited_alignment:
        open(os.path.join(config.clusters, cluster, "alignment_filtered_all"), "w").close()
        lengths.append(0)
        alignment_lengths_df.loc[len(alignment_lengths_df)] = lengths
        alignment_lengths_df.to_csv(os.path.join(config.clusters, cluster, "stats_lengths.csv"), index=False)
    else:
        edited_length = edited_alignment.get_alignment_length()
        lengths.append(edited_length)
        alignment_lengths_df.loc[len(alignment_lengths_df)] = lengths
        strain_idx = 0
        strains_num = len(os.listdir(config.data))
        while strain_idx < strains_num:
            if len(edited_alignment) > strain_idx:
                seq = edited_alignment[strain_idx]
                seq_strain_idx = int(seq.id.split("|")[0])
                if strain_idx < seq_strain_idx:
                    for padding_index in range(seq_strain_idx - strain_idx):
                        edited_alignment._records.insert(strain_idx + padding_index, SeqRecord(Seq(edited_length * "-"),
                                                                                               id="{} padding".format(
                                                                                                   strain_idx + padding_index)))
                    strain_idx += (seq_strain_idx - strain_idx + 1)
                    continue
                strain_idx += 1
            else:
                edited_alignment._records.append(SeqRecord(Seq(edited_length * "-"),
                                                           id="{} padding".format(strain_idx)))
                strain_idx += 1
        with open(os.path.join(config.clusters, cluster, cluster + "_alignment_filtered_without_gaps"), "w") as filtered_file:
            AlignIO.write(edited_alignment, filtered_file, "fasta")
        alignment_lengths_df.to_csv(os.path.join(config.clusters, cluster, "stats_lengths.csv"), index=False)


def concat_stats():
    stats_df = pd.DataFrame(columns=["cluster", "original", "edited"])
    for i, cluster in enumerate(os.listdir(config.clusters)):
        print("cluster {} start".format(i))
        cluster_stats = pd.read_csv(os.path.join(config.clusters, cluster, "stats_lengths.csv"))
        if stats_df.empty:
            stats_df = cluster_stats
        else:
            stats_df = pd.concat([stats_df, cluster_stats])
        print("cluster {} done".format(i))
    stats_df.to_csv(os.path.join(config.tables, "stats_lengths_without_gaps.csv"), index=False)


def filter_pool(clusters):
    with mp.Pool() as pool:
        pool.map(filter_alignment, clusters)


def filter_start():
    global alignment_dict
    alignment_dict = get_alignment_dict()
    clusters = os.listdir(config.clusters)
    filter_pool(clusters)
    concat_stats()


def align(mode):
    clusters_num = len(os.listdir(config.clusters))
    i = 1
    for cluster in os.listdir(config.clusters):
        print("cluster: {} start".format(cluster))
        print("cluster num {} of {} start".format(i, clusters_num))
        if mode == "mafft":
            mafft_args = ["mafft", "--auto", str(os.path.join(config.clusters, cluster, cluster + ".fasta"))]
            with open(os.path.join(config.clusters, cluster, cluster + "_alignment"), "w") as aligned_file:
                subprocess.run(mafft_args, stdout=aligned_file)
        elif mode == "gblocks":
            gblocks_args = "~/Gblocks_0.91b/Gblocks " + str(os.path.join(config.clusters, cluster, cluster + "_alignment")) + " -t=d -b5=a"
            with open(os.path.join(config.clusters, cluster, cluster + "_alignment_gblocks"), "w") as gblocks_file:
                subprocess.run(gblocks_args, shell=True, stdout=gblocks_file)


def main():
    arg = str(sys.argv[1])
    if arg == 'filter':
        filter_start()
        print("filter done")
    else:
        align(arg)
        print("align done")


if __name__ == "__main__":
    main()
