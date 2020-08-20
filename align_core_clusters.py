import subprocess
import os


root = "clusters"
for cluster in os.listdir(root):
    print("cluster: {} start".format(cluster))
    mafft_args = ["mafft", "--auto", str(os.path.join(root, cluster, cluster + ".fasta"))]
    with open(os.path.join(root, cluster, cluster + "_alignment"), "w") as aligned_file:
        subprocess.run(mafft_args, stdout=aligned_file)
    print("cluster: {} done".format(cluster))

print("script done")
