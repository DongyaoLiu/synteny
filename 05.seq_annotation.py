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
            if line_type == "gene" or line_type =="CDS":
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
ava_file = open(f"/lustre1/g/sbs_cgz/maf_synteny/{sys.argv[1]}", "r")
ana_file = open(f"/lustre1/g/sbs_cgz/maf_synteny/{sys.argv[2]}", "a")

seq_dict = {}
for line in ava_file.readlines():
    line = line.strip()
    line_elements = line.split("\t")
    ava_species = line_elements[0]
    ava_scaff = line_elements[1]
    ava_start = int(line_elements[2])
    ava_end = int(line_elements[3])
    ava_species2 = line_elements[5]
    ava_scaff2 = line_elements[6]
    ava_start2 = int(line_elements[7]) 
    ava_end2 = int(line_elements[8])
    if f"{ava_species}_{ava_scaff}" not in seq_dict.keys():
        seq_dict[f"{ava_species}_{ava_scaff}"] = [ava_start, ava_end]
    if f"{ava_species2}_{ava_scaff2}" not in seq_dict.keys():
        seq_dict[f"{ava_species2}_{ava_scaff2}"] = [ava_start, ava_end]
    if f"{ava_species}_{ava_scaff}" in seq_dict.keys():
        if seq_dict[f"{ava_species}_{ava_scaff}"][0] >= ava_start:
            seq_dict[f"{ava_species}_{ava_scaff}"][0] = ava_start
        if seq_dict[f"{ava_species}_{ava_scaff}"][1] <= ava_end:
            seq_dict[f"{ava_species}_{ava_scaff}"][1] = ava_end
    if f"{ava_species2}_{ava_scaff2}" in seq_dict.keys():
        if seq_dict[f"{ava_species2}_{ava_scaff2}"][0] >= ava_start2:
            seq_dict[f"{ava_species2}_{ava_scaff2}"][0] = ava_start2
        if seq_dict[f"{ava_species2}_{ava_scaff2}"][1] <= ava_end2:
            seq_dict[f"{ava_species2}_{ava_scaff2}"][1] = ava_end2

seq_file = open(f"{sys.argv[3]}", "a")
for ava_key, ava_val in seq_dict.items():
    ava_key_elements = ava_key.split("_", 1)
    seq_start = ava_val[0]
    seq_end = ava_val[1]
    seq_species = ava_key_elements[0]
    seq_scaff = ava_key_elements[1]
    seq_file.write(f"{seq_species}\t{seq_scaff}\t{seq_start}\t{seq_end}\n")
    exec(f"target_dict = {seq_species}_dict", globals())
    for key,val in target_dict.items():
        key_elements = key.split("\t")
        if key_elements[0] == seq_scaff:
            if "gene_coordinates" not in globals():
                if seq_start <= int(key_elements[3]) and seq_end >= int(key_elements[3]):
                    print(f"{seq_species}\t{key}\t{val}")
                    ana_file.write(f"{seq_species}\t{key}\t{val}\n")
                    if key_elements[1] == "CDS":
                        gene_coordinates = [int(key_elements[2]),int(key_elements[3])] 
                elif int(key_elements[3]) >= seq_end:
                    print(f"{seq_species}\t{key}\t{val}")
                    gene_name = val
                    ana_file.write(f"{seq_species}\t{key}\t{val}\n")
                    break
            else:
                if int(key_elements[2]) >= gene_coordinates[0] and int(key_elements[3]) <= gene_coordinates[1]:
                    print(f"{seq_species}\t{key}\t{val}")
                else:
                    del gene_coordinates

seq_file.close()





