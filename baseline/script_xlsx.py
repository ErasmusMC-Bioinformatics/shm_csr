import xlrd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Excel input file containing one or more sheets where column G has the gene annotation, H has the sequence id and J has the sequence")
parser.add_argument("--ref", help="Reference file")
parser.add_argument("--output", help="Output file")

args = parser.parse_args()

gene_column = 6
id_column = 7
seq_column = 8
LETTERS = [x for x in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]


refdic = dict()
with open(args.ref, 'r') as ref:
	currentSeq = ""
	currentId = ""
	for line in ref.readlines():
		if line[0] is ">":
			if currentSeq is not "" and currentId is not "":
				refdic[currentId[1:]] = currentSeq
			currentId = line.rstrip()
			currentSeq = ""
		else:
			currentSeq += line.rstrip()
	refdic[currentId[1:]] = currentSeq
	
currentSeq = ""
currentId = ""
with xlrd.open_workbook(args.input, 'r') as wb:
	with open(args.output, 'a') as o:
		for sheet in wb.sheets():
			if sheet.cell(1,gene_column).value.find("IGHV") < 0:
				print "Genes not in column " + LETTERS[gene_column] + ", skipping sheet " + sheet.name
				continue
			o.write(">>>" + sheet.name + "\n")
			outputdic = dict()
			for rowindex in range(1, sheet.nrows):
				ref = sheet.cell(rowindex, gene_column).value.replace(">", "")
				if ref in outputdic:
					outputdic[ref] += [(sheet.cell(rowindex, id_column).value.replace(">", ""), sheet.cell(rowindex, seq_column).value)]
				else:
					outputdic[ref] = [(sheet.cell(rowindex, id_column).value.replace(">", ""), sheet.cell(rowindex, seq_column).value)]
			#print outputdic
			
			for k in outputdic.keys():
				if k in refdic:
					o.write(">>" + k + "\n")
					o.write(refdic[k] + "\n")
					for seq in outputdic[k]:
						#print seq
						o.write(">" + seq[0] + "\n")
						o.write(seq[1] + "\n")
				else:
					print k + " not in reference, skipping " + k
