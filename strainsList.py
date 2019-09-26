

strains_list = []
with open("assembly_summary.txt") as parent_dir:
    for line in parent_dir:
        if not line.startswith("#"):
            line_tab = line.split('\t')
            strain = line_tab[19].split(":")
            strain = strain[1][2:].split("/")[-1]
            strains_list.append(strain)

with open("strains_list.txt", 'w') as strains_file:
    strains_file.write("Number of strains: {}\n".format(len(strains_list)))
    strains_file.write('\n'.join(strains_list))
