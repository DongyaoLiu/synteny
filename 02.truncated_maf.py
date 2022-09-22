import linecache
import sys
import re


def find_close_lines(index_location, maf_scaff, start, end):
    index_file = open(f"{index_location}", "r")
    line_out = []
    for line in index_file.readlines():
        line = line.strip()
        if re.match("#", line):
            pass
        else:
            line_elements = line.split("\t")
            scaff =line_elements[0]
            align_start = int(line_elements[1])
            align_end = int(line_elements[2])
            strand = line_elements[3]
            line_start = line_elements[4]
            line_end = line_elements[5]
            if maf_scaff == scaff:
                if len(line_out) == 0:
                    if align_start <= start and start <= align_end:
                        line_out.append(int(line_start))
                if end <= align_end:
                    line_out.append(int(line_end))
                    return line_out


##main loop
#################################################################################################################

maf_index_file = "/lustre1/u/dongyao/maf_exp/elegans_ref.maf.index"
maf_file = "/lustre1/u/dongyao/maf_exp/elegans_ref.maf"
maf_line_coordinates = find_close_lines(maf_index_file, str(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
print("###" + f"{maf_line_coordinates}" + " " + "line coordinates")
maf = open(maf_file, "r")
line_number = 0
for line in maf.readlines():
    line_number = line_number + 1
    if line_number >= maf_line_coordinates[0] and line_number <= maf_line_coordinates[1]:
        print(line.strip())
    elif line_number > maf_line_coordinates[1]:
        break
