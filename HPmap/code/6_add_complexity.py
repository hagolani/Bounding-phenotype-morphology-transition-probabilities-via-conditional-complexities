import numpy as np
import sys
import os

def KC_LZ(string):
    n = len(string)
    s = '0' + string
    c = 1
    l = 1
    i = 0
    k = 1
    k_max = 1
    stop = 0

    while stop == 0:
        if s[i + k] != s[l + k]:
            if k > k_max:
                k_max = k

            i = i + 1

            if i == l:
                c = c + 1
                l = l + k_max

                if l + 1 > n:
                    stop = 1
                else:
                    i = 0
                    k = 1
                    k_max = 1
            else:
                k = 1
        else:
            k = k + 1

            if l + k > n:
                c = c + 1
                stop = 1

    return c

def calc_KC(s):
    L = len(s)
    if s == '0' * L or s == '1' * L:
        return np.log2(L)
    else:
        return np.log2(L) * (KC_LZ(s) + KC_LZ(s[::-1])) / 2.0

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        #header = infile.readline().strip()
        outfile.write("ParentPhenotyope MutantPhenotype count proportion KC_col1 KC_col4 Conditional_Complexity\n")

        for line in infile:
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # Skip lines that don't have enough columns

            col1 = parts[0]
            col4 = parts[1]

            KC_col1 = calc_KC(col1)
            KC_col4 = calc_KC(col4)

            k1 = round(KC_col1, 1)
            str_arrlist_con = [col1 + col4 for _ in range(len(col4))]  # As per the given example, using col4 multiple times to form concatenated strings
            Z = [round(calc_KC(s), 1) for s in str_arrlist_con]
            Z1 = np.array(Z) - k1

            conditional_complexity = np.mean(Z1)

            outfile.write(line.strip() + f" {KC_col1:.2f} {KC_col4:.2f} {conditional_complexity:.2f}\n")

# Replace 'input.txt' and 'output.txt' with your actual file names
argument = sys.argv[1]

# Split the input filename into name and extension
name, ext = os.path.splitext(argument)
#URDDDDDLDLULUURDRUULLUR  dataToAll
#UUURDDDDLLULLURULURURDDD HP_
#URDDDDDLLLULURRDRUULULD  real
#UUUUURDDDDDDLLULURUULURU 
# Create the output filename by appending '_KC' to the name
output_file = name + "_KC" + ext

#input_file = "HP_output_data/HP_25_possible_new_unique_str_count_prob_"+argument+".txt"
#output_file = "HP_output_data/dataNeigh"+argument+"_new_KC.txt"

process_file(argument, output_file)









