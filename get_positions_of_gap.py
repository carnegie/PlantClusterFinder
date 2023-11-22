import re
import os
from argparse import ArgumentParser
from operator import itemgetter
from itertools import groupby


def regex_search(string_pattern_list, string_input):
    string_pattern = re.compile("[" + ''.join(string_pattern_list) + "]")
    return [m.start() for m in re.finditer(string_pattern, string_input)]


# From BioPython
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)


def groupbyseq(input_list, output_path, amino_acid, chromosome):
    if os.path.exists(output_path):
        append_write = 'a'  # append if already exists
    else:
        append_write = 'w'  # make a new file if not

    with open(output_path, append_write) as op:
        for k, g in groupby(enumerate(input_list), lambda x: x[1] - x[0]):
            try:
                cur_group = list(map(itemgetter(1), g))
                start_index = str(cur_group[0] + 1)
                length = str(len(cur_group))
                op.write(' '.join([chromosome, amino_acid, start_index, length]) + "\n")
            except (TypeError or IndexError) as e:
                print("Error", e, list(map(itemgetter(1), g)))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', '--input_path', dest="input_path", required=True)
    parser.add_argument('-o', '--output_path', dest="output_path", required=True)

    args = parser.parse_args()

    try:
        os.remove(args.output_path)
    except OSError:
        pass
    with open(args.input_path, 'r') as fp, open(args.output_path, 'w') as op:
        for seq_name, seq_aa in read_fasta(fp):
            fasta_id = re.split(r'[\t\s|]', seq_name)[0].replace('>', '')
            indices = regex_search("N", seq_aa.upper())
            groupbyseq(indices, args.output_path, "N", fasta_id)
            indices = regex_search(["A", "T", "C", "G", "*"], seq_aa.upper())
            groupbyseq(indices, args.output_path, "A", fasta_id)
