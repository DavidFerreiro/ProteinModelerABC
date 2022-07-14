import argparse
from Bio import SeqIO
import os



def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1 # use start += 1 to find overlapping matches

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("--input", default="/dev/stdin")
	parser.add_argument("--temp", default="/dev/stdin")
	parser.add_argument("--chain", default="/dev/stdin")
	parser.add_argument("--output", default="/dev/stdout")

	return parser.parse_args()


def main():

	args = parse_args()

	with open(args.temp, 'r') as pdb_file:
		for record in SeqIO.parse(pdb_file, 'pdb-atom'):
			if record.annotations["chain"] == str(args.chain):
				PDB_Name = ('>' + record.id[:-2])
				PDB_Seq = (record.seq)
				PDB_Seq = str(PDB_Seq)
				PDB_Seq = PDB_Seq.replace('X', '')

	with open('T.fasta', 'w') as T:
		T.write(PDB_Name + '\n')
		T.write(PDB_Seq + '\n')
		with open(args.input) as a:
			for line in a:
				T.write(str(line))


	os.system('muscle -align T.fasta -output T-a.fasta')

	i=0
	t=0
	sequence=''

	with open('T-a.fasta', 'r') as T:
		for x in T:
			if x[0] == '>':
				i = i + 1
				if sequence != '':
					sequence = ''
				else:
					pass
			else:
				sequence = sequence + str(x[:-1])
		NumSeq = i

	with open('T-a.fasta', 'r') as T:
		with open('T-f.fasta', "w") as F:
			sequence = ''
			for line in T:
				if line[0] == '>':
					if sequence != '':
						#if sequence != PDB_Seq:
						F.write(sequence + '\n')
						t = t + 1
						sequence = ''
						#else:
							#sequence = ''
					F.write(line)
				else:
					sequence = sequence + str(line[:-1])

			if t == (int(NumSeq) - 1):
				F.write(sequence + '\n')

	i = 0
	t = 0
	NumAA = 0
	sequence = ''
	with open(args.input) as handle:
		for x in handle:
			if x[0] == '>':
				i = i + 1
				if sequence != '':
					sequence = ''
				else:
					pass
			else:
				sequence = sequence + str(x[:-1])
		NumSeq = i


		with open(args.output, "w") as output_handle:
			output_handle.write(str(NumSeq) + ' ' + str(len(PDB_Seq)) + '\n')


	q = 0
	ip = 0
	ind=[]
	with open('T-f.fasta') as Q:
		for x in Q:
			q = q + 1
			if x[:-1] == PDB_Name:
				ip = q +1
			if ip == q:
				for i in range(len(x)):
					if x.startswith('-', i):
						ind.append(i)
	#print(ind)
	q = 0
	ip = 0
	m = 0
	with open('T-f.fasta', 'r') as Q:
		with open('T-l.fasta', 'w') as L:
			for x in Q:
				q = q + 1
				if x[0] == '>':
					if x[:-1] == PDB_Name:
						ip = q + 1
					else:
						L.write(str(x))
				else:
					if ip == q:
						pass
					else:
						if str(len(ind)) != '0':
							for y in range(0, len(ind)):
								#print(len(ind))
								if len(ind) == 1:
									L.write(str(x[0:int(ind[y])]))
									#print('[0:' + str(ind[y]) + ']')
									m = int(ind[y]) + 1
									L.write(str(x[m:]))
									#print('[' + str(m) + ':]')
								else:
									if y == 0:
										if ind[y] == 0:
											#print('1:' + str(ind[1]))
											L.write(str(x[1:int(ind[1])]))
										else:
											#print('0:' + str(ind[y]))
											L.write(str(x[0:int(ind[y])]))
											m = int(ind[y]) +1
											L.write(str(x[m:int(ind[y+1])]))
									elif y == (len(ind)-1):
										m = int(ind[y]) + 1
										#print( str(m)+ ':]')
										L.write(str(x[m:]))
									else:
										m = int(ind[y]) +1
										#print(m)
										#print(str(m) + ':' + str(ind[y+1]))
										L.write(str(x[m:int(ind[y+1])]))
						else:
							L.write(str(x))


	with open('T-l.fasta') as handle:
		with open(args.output, "a") as output_handle:
			sequence = ''
			for line in handle:
				if line[0] == '>':
					if sequence != '':
						output_handle.write(sequence + '\n')
						t = t + 1
						sequence = ''
					output_handle.write(line[1:9] + '  ')
				else:
					sequence = sequence + str(line[:-1])

			if t == (int(NumSeq) - 1):
				output_handle.write(sequence + '\n')


	os.system('rm T-f.fasta')
	os.system('rm T-a.fasta')
	os.system('rm T-l.fasta')
	os.system('rm T.fasta')




if __name__ == "__main__":
	main()
