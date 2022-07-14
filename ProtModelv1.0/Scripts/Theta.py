from Bio import AlignIO
import sys
import argparse

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("--input", default="/dev/stdin")
	parser.add_argument("--output", default="/dev/stdout")

	return parser.parse_args()

def main():
    args = parse_args()


    try:
        file = open(args.input)
        print(args.input + ' exists')
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

    align = AlignIO.read(str(args.input),"phylip-sequential")

    i = 0
    p = 1
    q = 0
    SeqIDT = 0

    for x in range(0, (len(align) - 1)):
        for y in range(p, len(align)):
            q = q + 1
            for e in range(0, len(align[0])):
                if align[x][e] == "-":
                    i = i + 0
                elif align[y][e] == "-":
                    i = i + 0
                elif align[x][e] != align[y][e]:
                    i = i + 1
            SeqID = (len(align[0]) - i) / len(align[0]) * 100
            SeqIDT = SeqIDT + SeqID
            #print(SeqID)
            i = 0

        p = p + 1

    SeqID = SeqIDT/q



    SR = float(Theta) / 2 / float(Dip) / len(align[0]) / float(Ne)

    print('\nSeqID  = ' + str(SeqID))
    print('Substitution Rate per site = ' + "{:e}".format(SR))

if __name__ == "__main__":
	main()