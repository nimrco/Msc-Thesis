import os
from Bio import AlignIO


root = "clusters"
alignments_num = len(os.listdir(root))
i = 1
concat_alignment = None
for cluster in os.listdir(root):
    print("alignment: {} start".format(cluster))
    print("alignmnet num {} of {} start".format(i, alignments_num))
    alignment = AlignIO.read(open(os.path.join(root, cluster, cluster + "_alignment_filtered")), "fasta")
    if not concat_alignment:
        concat_alignment = alignment[:, :]
    else:
        concat_alignment += alignment[:, :]
    print("alignment num {} of {} done".format(i, alignments_num))
    i += 1
AlignIO.write(concat_alignment, open("concatenated_alignmnet", "w"), "fasta")

print("script done")
