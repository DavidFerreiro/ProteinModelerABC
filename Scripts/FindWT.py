#!/usr/bin/env python
import numpy as np
import sys
import argparse
from Bio import AlignIO
from collections import Counter
np.set_printoptions(threshold=sys.maxsize)
arr = []



def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("--input", default="/dev/stdin")

	return parser.parse_args()


def main():
    args = parse_args()
    aln = AlignIO.read(args.input, "fasta")
    
    sequence = []
    p = aln.get_alignment_length()
    for x in range(0, p):
        max = []
        for y in range(0, len(aln)):
            let = aln[y][x]
            max.append(let)

        counter = Counter(max)
        try:
            first, second, *_, last = counter.most_common()
            print(x+1,first, second, last)
            sequence.append(str(first)[2])
        except:
            first= counter.most_common()
            print(first)
            print(x+1,first)
            sequence.append(str(first)[3])
    print('\nWT sequence:')
    print(''.join(sequence))
    print('')


if __name__ == "__main__":
	main()