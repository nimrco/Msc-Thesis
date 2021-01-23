import os
from Bio import AlignIO
import config


def get_alignment_filtered_dict():
    alignment_dict_local = {}
    for cluster_local in os.listdir(config.clusters):
        with open(os.path.join(config.clusters, cluster_local, cluster_local + "_alignment_filtered_without_gaps")) as alignment_file:
            alignment_local = AlignIO.read(alignment_file, "fasta")
            for record in alignment_local:
                record.id = record.id.split('|')[0]
                record.description = record.id
            alignment_dict_local[cluster_local] = alignment_local
    print("alignment_dict")
    return alignment_dict_local


alignment_dict = get_alignment_filtered_dict()
alignments_num = len(os.listdir(config.clusters))
concat_alignment = None
for i, cluster in enumerate(os.listdir(config.clusters)):
    print("alignment: {} start".format(cluster))
    print("alignmnet num {} of {} start".format(i, alignments_num))
    alignment = alignment_dict[cluster]
    if not concat_alignment:
        concat_alignment = alignment[:, :]
    else:
        concat_alignment += alignment[:, :]
    print("alignment num {} of {} done".format(i, alignments_num))
with open(os.path.join(config.seq_files, "concatenated_alignmnet_id"), "w") as concatenated_file:
    AlignIO.write(concat_alignment, concatenated_file, "fasta")
print("script done")
