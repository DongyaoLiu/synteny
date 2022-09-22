import re 
import os
import sys
import itertools


##########################################################
##########################################################
##  split function
##########################################################
##########################################################


def split_species_scaff(maf_line):
    species_scaff = re.search(r"caenorhabditis_([a-z]+)[.]{1}(\S*)", maf_line)
    species = species_scaff.group(1)
    scaff = species_scaff.group(2)
    species_name_list = [species, scaff, f"caenorhabditis_{species}.{scaff}"]
    return species_name_list

def split_start_end_strand(maf_line):  #####turn to gff 1 based coordinates
    maf_line_elements = maf_line.strip().split("\t")
    strand = maf_line_elements[4]
    if strand == "+":
        start = int(maf_line_elements[2]) + 1
        end = int(maf_line_elements[2]) + int(maf_line_elements[3])    ### conclude the end index ##full closed
        match_number =  int(maf_line_elements[3])
    else:
        full_len = int(maf_line_elements[5])
        last_nuc = int(maf_line_elements[5]) - 1
        neg_start = int(maf_line_elements[2])
        end = full_len - neg_start    ### full closed
        start = end - int(maf_line_elements[3]) + 1
    return [start, end, strand]

class split_speciesDICT_key_info:
    def __init__(self,dict_key):
        key_elements = dict_key.split("$")
        self.scaff = key_elements[0]
        self.start = int(key_elements[1])
        self.end = int(key_elements[2])
        self.strand = key_elements[3]
        self.block_number = key_elements[4]


def synteney_analysis(list_name_list, two_compare_object, ava_location, out_seq_dict):
    ava_file = open(f"{ava_location}", "a")
    Obj1 = two_compare_object[0]
    Obj2 = two_compare_object[1]
    out_ava_list = []
    for block_name in list_name_list:
        exec(f"Target_list = {block_name}",globals())
        for line in Target_list:
            start_end_strand_list = split_start_end_strand(line)
            species_scaff= split_species_scaff(line)
            if species_scaff[0] == Obj1 or species_scaff[0] == Obj2:
                out_ava_list = out_ava_list + [species_scaff[0], 
                        species_scaff[1], 
                        str(start_end_strand_list[0]), 
                        str(start_end_strand_list[1]), 
                        start_end_strand_list[2]]
                if species_scaff[2] not in out_seq_dict.keys():
                    out_seq_dict[f"{species_scaff[2]}"] = [start_end_strand_list[0]] ### record seq start
                elif len(out_seq_dict[f"{species_scaff[2]}"]) == 1:
                    if out_seq_dict[f"{species_scaff[2]}"][0] >= start_end_strand_list[0]: ##compare seq start
                        out_seq_dict[f"{species_scaff[2]}"][0] = start_end_strand_list[0]  ## change seq start
                    out_seq_dict[f"{species_scaff[2]}"].append(start_end_strand_list[1]) ### record seq end
                else:
                    if out_seq_dict[f"{species_scaff[2]}"][1] <= start_end_strand_list[1]:    #compare end 
                        out_seq_dict[f"{species_scaff[2]}"][1] = start_end_strand_list[1]### change seq end
                    if out_seq_dict[f"{species_scaff[2]}"][0] >= start_end_strand_list[0]: ##compare seq start
                        out_seq_dict[f"{species_scaff[2]}"][0] = start_end_strand_list[0]  ## change seq start
            if len(out_ava_list) >= 6:
                out_line = "\t".join(out_ava_list)
                ava_file.write(f"{out_line}\n")
                out_ava_list = out_ava_list[0:5]
        out_ava_list = []
    ava_file.close()
    return out_seq_dict
##########################################################
##########################################################
##                     Maf read
##########################################################
##########################################################
##seq_table rownames = track_names colnames = length of track 

maf = open(f"{sys.argv[1]}", "r")

block_number = 1
block_1 = []
block_name_list = ["block_1"]
species_name_list = []
s_number = 0
for line in maf.readlines():
    line = line.strip()
    if line == "":
        if s_number != 0:
            block_number = block_number + 1
            exec(f"block_{block_number} = []", globals())
            exec(f"block_name_list.append(block_{block_number})",globals())
        block_dict = {}
        s_number = 0
    elif re.match("#", line):
        print(line)
    elif re.match("s", line):
        s_number = s_number + 1 
        if re.search("caenorhabditis", line):
            species = re.search(r"caenorhabditis_([a-zA-Z]+)",line).group(1)
            if species not in species_name_list :
                species_name_list.append(species)
            exec(f"block_{block_number}.append(\"{line}\")", globals())


##########################################################
##########################################################
##########################################################
ava_location = sys.argv[2]
seq_location = sys.argv[3]

seq_start_end_dict = {}
for combination in itertools.combinations(species_name_list,2):
    print(combination, "start data extract")
    seq_start_end_dict = synteney_analysis(block_name_list, combination, ava_location, seq_start_end_dict)
    
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################

seq_file = open(seq_location, "a")
seq_file.write("species\tscaffold\tstart\tend\n")
for key,val in seq_start_end_dict.items():
    species = split_species_scaff(key)[0]
    scaff = split_species_scaff(key)[1]
    start = val[0]
    end = val[1]
    seq_file.write(f"{species}\t{scaff}\t{start}\t{end}\n")





