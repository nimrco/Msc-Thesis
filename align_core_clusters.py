import subprocess
import os


root = "test_dir"
mafft_args = ["mafft", "--auto"]
for cluster in os.listdir(root):
    with open(os.path.join(root, cluster, cluster + ".fasta")) as cluster_file:
        with open(os.path.join(root, cluster, cluster + "_alignment"), "w") as aligned_file:
            subprocess.run(mafft_args, stdin=cluster_file, stdout=aligned_file)
