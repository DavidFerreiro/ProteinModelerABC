from Bio import AlignIO
import sys
import argparse

#Example python Theta.py --input example.phy

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("--input", default="/dev/stdin")
	parser.add_argument("--output", default="/dev/stdout")

	return parser.parse_args()


def perc_identity(aln):
    # For every position, for every sequence, we take the aa and compare it with the aa of same position of the next sequence and if it is different sums 1 in a list
    total_dif = 0
    i = 0
    p = aln.get_alignment_length()
    for x in range (0, p):
        c = 1
        for a in range(0, (len(aln)-1)):
            s = aln[a][x]
            for b in range(c, len(aln)):
                t = aln[b][x]
                i = i + 1
                if s == t:
                    total_dif = total_dif + 1
            c = c + 1
    seqid = total_dif / i
    return '\n - Sequence identity: ' + str(seqid)


def main():
    args = parse_args()
    align = AlignIO.read(str(args.input),"phylip-sequential")
    print('Substitution rate prior should be chosen considering the alignment sequence identity to ensure that it encompass the real data. As a recomendation:')
    print('Highly conserved alignments (SeqID > 70%) usually works better with low limits in the substition rate prior distirbution (θ < 300)')
    print('Medium conserved alignments (SeqID < 70% >50%) usually works better with medium limits in the substition rate prior distirbution (θ > 300)')
    print('Poorly conserved alignments (SeqID < 50%) usually works better with high limits in the substition rate prior distirbution (θ > 500)')
    print(perc_identity(align) + '\n')

    try:
        file = open(args.input)
        file.close()
    except FileNotFoundError:
        print('Sorry the file we\'re looking for doesn\' exist')
        exit()

    try:
        Theta = input('Enter your desired Theta: ')
    except ValueError:
        print("Sorry, I didn't understand that.")

    try:
        if Theta.isnumeric() == True:
            print(Theta + ' seems OK')
        else:
            sys.exit()
    except SystemExit:
        print('Sorry we are looking for a number value as Theta')
        exit()

    try:
        Dip = input('Is your organism Diploid or Haplod: ')
    except ValueError:
        print('Sorry , I didn\'t understand that.')

    try:
        if Dip == 'Diploid' or Dip == 'Haploid' or Dip == 'diploid' or Dip == 'haploid':
            print(Dip + ' seems OK')
        else:
            sys.exit()
    except SystemExit:
        print('Sorry , I didn\'t understand that.')
        exit()

    if Dip == 'Diploid' or Dip == 'diploid':
        Dip = 2
    else:
        Dip = 1

    try:
        Ne = input('Enter your Population Size: ')
    except ValueError:
        print('Sorry, I didn\'t understand that.')

    try:
        if Ne.isnumeric() == True:
            print(Ne + ' seems OK')
        else:
            sys.exit()
    except SystemExit:
        print('Sorry we are looking for a number value as Population Size')
        exit()

    SR = float(Theta) / 2 / float(Dip) / len(align[0]) / float(Ne)

    print(' - Substitution Rate per site = ' + "{:e}".format(SR))

if __name__ == "__main__":
	main()