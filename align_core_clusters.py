import subprocess
import os
import sys


def align(mode):
    root = "clusters"
    clusters_num = len(os.listdir(root))
    i = 1
    for cluster in os.listdir(root):
        print("cluster: {} start".format(cluster))
        print("cluster num {} of {} start".format(i, clusters_num))
        if mode == "mafft":
            mafft_args = ["mafft", "--auto", str(os.path.join(root, cluster, cluster + ".fasta"))]
            with open(os.path.join(root, cluster, cluster + "_alignment"), "w") as aligned_file:
                subprocess.run(mafft_args, stdout=aligned_file)
        else:
            gblocks_args = ["Gblocks", str(os.path.join(root, cluster, cluster + "_alignment")), "-t=d", "-b5=a"]
            with open(os.path.join(root, cluster, cluster + "_alignment_gblocks"), "w") as gblocks_file:
                subprocess.run(gblocks_args, stdout=gblocks_file)
        print("cluster: {} done".format(cluster))
        print("cluster num {} of {} done".format(i, clusters_num))
        i += 1

    print("script done")


def main():
    align(str(sys.argv[1]))


if __name__ == "__main__":
    main()
