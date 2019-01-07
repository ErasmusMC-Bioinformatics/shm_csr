#import xlrd #avoid dep
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Excel input file containing one or more sheets where column G has the gene annotation, H has the sequence id and J has the sequence")
parser.add_argument("--ref", help="Reference file")
parser.add_argument("--output", help="Output file")
parser.add_argument("--id", help="ID to be used at the '>>>' line in the output")

args = parser.parse_args()

print "script_imgt.py"
print "input:", args.input
print "ref:", args.ref
print "output:", args.output
print "id:", args.id

refdic = dict()
with open(args.ref, 'rU') as ref:
	currentSeq = ""
	currentId = ""
	for line in ref:
		if line.startswith(">"):
			if currentSeq is not "" and currentId is not "":
				refdic[currentId[1:]] = currentSeq
			currentId = line.rstrip()
			currentSeq = ""
		else:
			currentSeq += line.rstrip()
	refdic[currentId[1:]] = currentSeq

print "Have", str(len(refdic)), "reference sequences"

vPattern = [r"(IGHV[0-9]-[0-9ab]+-?[0-9]?D?\*\d{1,2})"]#,
#						r"(TRBV[0-9]{1,2}-?[0-9]?-?[123]?)",
#						r"(IGKV[0-3]D?-[0-9]{1,2})",
#						r"(IGLV[0-9]-[0-9]{1,2})",
#						r"(TRAV[0-9]{1,2}(-[1-46])?(/DV[45678])?)",
#						r"(TRGV[234589])",
#						r"(TRDV[1-3])"]

#vPattern = re.compile(r"|".join(vPattern))
vPattern = re.compile("|".join(vPattern))

def filterGene(s, pattern):
    if type(s) is not str:
        return None
    res = pattern.search(s)
    if res:
        return res.group(0)
    return None



currentSeq = ""
currentId = ""
first=True
with open(args.input, 'r') as i:
	with open(args.output, 'a') as o:
		o.write(">>>" + args.id + "\n")
		outputdic = dict()
		for line in i:
			if first:
				first = False
				continue
			linesplt = line.split("\t")
			ref = filterGene(linesplt[1], vPattern)
			if not ref or not linesplt[2].rstrip():
				continue
			if ref in outputdic:
				outputdic[ref] += [(linesplt[0].replace(">", ""), linesplt[2].replace(">", "").rstrip())]
			else:
				outputdic[ref] = [(linesplt[0].replace(">", ""), linesplt[2].replace(">", "").rstrip())]
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
