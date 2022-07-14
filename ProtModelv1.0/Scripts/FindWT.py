#!/usr/bin/env python
import numpy as np
import sys
import argparse
np.set_printoptions(threshold=sys.maxsize)
arr = []



def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("--input", default="/dev/stdin")

	return parser.parse_args()


def main():
    args = parse_args()
    with open(args.input, 'r') as sequencia:
	    seqs = sequencia.readlines()
    
    sequence = ''
    t = 0
    for x in seqs:
        if x[0] == '>':
            if sequence != '':
                arr.append(list(sequence))
                print(list(sequence))
                t = t + 1
                sequence = ''
        else:
            sequence = sequence + str(x[:-1])


    df = np.array( arr, dtype="object" )
    df = df.T

    i = 1
    final_seq = []
    for seq in df:
        seq_str = ''.join(str(x) for x in seq)
        print(i, max(set(seq_str), key=seq_str.count), {y: seq_str.count(y) for y in set(seq_str)})
        final_seq.append(max(set(seq_str), key=seq_str.count))
        i = i+1

    print('\nWT sequence:')
    print(''.join(final_seq))

if __name__ == "__main__":
	main()