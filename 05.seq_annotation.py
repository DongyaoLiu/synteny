import re
import sys
import os
import copy




#################################### input gtf_file into dict
def gtf2dict(gtf_location, species_name):
    exec(f"{species_name}_dict = {{}}" ,globals())
    gtf_file = open(f"{gtf_location}", "r")
    for line in gtf_file.readlines():
        line = line.strip()
        if re.match("#", line):
            pass
        else:
            line_elements = line.split("\t")
            scaff = line_elements[0]
            line_type = line_elements[2]
            if line_type == "gene" or "exon":
                start = line_elements[3]
                end = line_elements[4]
                strand = line_elements[6]
                gene_name = re.search(r"gene_id \"([-A-Za-z0-9._]+)\";", line_elements[8]).group(1)
                exec(f"{species_name}_dict[\"{scaff}\t{line_type}\t{start}\t{end}\t{strand}\"] = \"{gene_name}\"",globals())
    gtf_file.close()
    print(f"gtf2dict {species_name} done")



gtf_location = "/home/dongyao/synteny"
gtf_file_list = os.listdir(f"{gtf_location}")
for file in gtf_file_list:
    if re.search(r"caenorhabditis.+gtf", file):
        species = re.match(r"caenorhabditis_([a-z]+).",file).group(1)
        gtf2dict(f"{gtf_location}/{file}", species)

#################################### annotation start 
seq_file = open(f"/lustre1/g/sbs_cgz/maf_synteny/{sys.argv[1]}", "r")
annotation_synteny_file =open(f"/lustre1/u/dongyao/maf_exp/{sys.argv[2]}", "a")
line_number = 0
for line in seq_file.readlines():
    line_number = line_number + 1
    if line_number == 1:
        continue
    line = line.strip()
    line_elements = line.split("\t")
    seq_species = line_elements[0].split("_")[1]
    seq_scaff = line_elements[1]
    seq_start = int(line_elements[2])
    seq_end = int(line_elements[3])
    seq_strand = line_elements[4]
    exec(f"target_dict = {maf_species}_dict", globals())
    for key,val in target_dict.items():
        key_elements = key.split("\t")
        if key_elements[0] == seq_scaff:
            if seq_start <= key_elements[3] and seq_end >= key_elements[3]:
                print(key)
            elif key_elements[3] >= 250:
                print(key)
                break











