import sys
import re


maf_file  = open(f"/lustre1/u/dongyao/maf_exp/{sys.argv[1]}", "r")

s_number = 0
line_num = 0
for line in maf_file.readlines():
    line = line.strip()
    line_num = line_num + 1
    if re.match("#",line) or line == "":
        if s_number != 0:
            ref_info.append(str(line_num))
            print("\t".join(ref_info))
        s_number = 0
        ref_info = []
        continue
    elif re.match("a", line):
        pass
    elif re.match("s", line):
        s_number = s_number + 1
        if s_number == 1:
           line_elements = line.split("\t")
           scaff = line_elements[1].split(".")[1]
           start = line_elements[2]
           end = int(line_elements[2]) + int(line_elements[3]) -1
           strand = line_elements[4]
           ref_info = [scaff, start, str(end), strand, str(line_num)]

