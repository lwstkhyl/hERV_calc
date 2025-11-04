#!/opt/apps/python3/bin/python3

import sys
import re


for ln in sys.stdin:
    items = ln.split("\t")
    read_name = items[0]
    chromosome = items[2]
    problematic_line = 0
    try:
        flag = "{0:{fill}12b}".format(int(items[1]), fill='0')
        strand_indic = [x for x in flag][6]
        if strand_indic == '0':
            strand = "+"
        else:
            strand = "-"
        start = int(items[3]) -1
        score = re.sub(r'.*:', '', [s for s in items if "NM:i" in s][0])
        flag = "{0:{fill}12b}".format(int(items[1]), fill='0')
        strand_indic = [x for x in flag][6]
        cigar = items[5]
        cigar=items[5]
        blocks = cigar.split("N")
        block_start = [0]
        block_size = []
        block_number = len(blocks)
        counter = 1
        total_size = 0
        for b in blocks:
            b_size = re.compile("[A-Z]").split(b)
            b_type = re.compile("\d").split(re.sub(r'\d',"", b))
            b_n = range(len(b_type))
            M_size = 0
            for i in b_n: 
                if b_type[i] == "M"or b_type[i] == "D":
                    M_size += int(b_size[i])
            block_size.append(M_size)
            if counter < block_number:
                total_size += M_size + int(b_size[-1])
                block_start.append(total_size)
            counter = counter +1
        end = start + block_start[-1] + block_size[-1]
        block_size_str = ",".join([str(x) for x in block_size])
        block_start_str = ",".join([str(x) for x in block_start])
        ln_items = [chromosome, str(start), str(end), read_name, score, strand, str(start), \
                    str(end), "0",  str(block_number), block_size_str, block_start_str]
        ln_output = "\t".join(ln_items)
        sys.stdout.write(ln_output + "\n")
    except:
        problematic_line += 1
        continue
