import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="The 1_Summary file of an IMGT zip file")
parser.add_argument("--fasta", help="The output fasta file")

args = parser.parse_args()

infile = args.input
fasta = args.fasta

with open(infile, 'r') as i, open(fasta, 'w') as o:
	first = True
	id_col = 0
	seq_col = 0
	no_results = 0
	no_seqs = 0
	passed = 0
	for line in i:
		splt = line.split("\t")
		if first:
			id_col = splt.index("Sequence ID")
			seq_col = splt.index("Sequence")
			first = False
			continue
		if len(splt) < 5:
			no_results += 1
			continue
		
		ID = splt[id_col]
		seq = splt[seq_col]
		
		if not len(seq) > 0:
			no_seqs += 1
			continue
		
		o.write(">" + ID + "\n" + seq + "\n")
		passed += 1
			
	print("No results:", no_results)
	print("No sequences:", no_seqs)
	print("Written to fasta file:", passed)
