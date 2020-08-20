import subprocess
import os


root = "test_dir"
for cluster in os.listdir(root):
    mafft_args = ["mafft", "--auto", str(os.path.join(root, cluster, cluster + ".fasta"))]
    with open(os.path.join(root, cluster, cluster + "_alignment"), "w") as aligned_file:
        subprocess.run(mafft_args, stdout=aligned_file)
