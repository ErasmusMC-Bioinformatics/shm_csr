import re
import argparse
import time
starttime= int(time.time() * 1000)

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="The 1_Summary file from an IMGT zip file")
parser.add_argument("--output", help="The annotated output file to be merged back with the summary file")

args = parser.parse_args()

infile = args.input
#infile = "test_VH-Ca_Cg_25nt/1_Summary_test_VH-Ca_Cg_25nt_241013.txt"
output = args.output
#outfile = "identified.txt"

dic = dict()
total = 0


first = True
IDIndex = 0
seqIndex = 0

with open(infile, 'r') as f: #read all sequences into a dictionary as key = ID, value = sequence
	for line in f:
		total += 1
		linesplt = line.split("\t")
		if first:
			print "linesplt", linesplt
			IDIndex = linesplt.index("Sequence ID")
			seqIndex = linesplt.index("Sequence")
			first = False
			continue
		
		ID = linesplt[IDIndex]
		if len(linesplt) < 28: #weird rows without a sequence
			dic[ID] = ""
		else:
			dic[ID] = linesplt[seqIndex]
			
print "Number of input sequences:", len(dic)

#old cm sequence: gggagtgcatccgccccaacccttttccccctcgtctcctgtgagaattccc
#old cg sequence: ctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctgggggcacagcggccctgggctgcctggtcaaggactacttccccgaaccggtgacggtgtcgtggaactcaggcgccctgaccag

#lambda/kappa reference sequence
searchstrings = {"ca": "catccccgaccagccccaaggtcttcccgctgagcctctgcagcacccagccagatgggaacgtggtcatcgcctgcctgg",
                 "cg": "ctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctgggggcacagcggcc",
                 "ce": "gcctccacacagagcccatccgtcttccccttgacccgctgctgcaaaaacattccctcc",
                 "cm": "gggagtgcatccgccccaacc"} #new (shorter) cm sequence

compiledregex = {"ca": [],
                 "cg": [],
                 "ce": [],
                 "cm": []}

#lambda/kappa reference sequence variable nucleotides
ca1 = {38: 't', 39: 'g', 48: 'a', 49: 'g', 51: 'c', 68: 'a', 73: 'c'}
ca2 = {38: 'g', 39: 'a', 48: 'c', 49: 'c', 51: 'a', 68: 'g', 73: 'a'}
cg1 = {0: 'c', 33: 'a', 38: 'c', 44: 'a', 54: 't', 56: 'g', 58: 'g', 66: 'g', 132: 'c'}
cg2 = {0: 'c', 33: 'g', 38: 'g', 44: 'g', 54: 'c', 56: 'a', 58: 'a', 66: 'g', 132: 't'}
cg3 = {0: 't', 33: 'g', 38: 'g', 44: 'g', 54: 't', 56: 'g', 58: 'g', 66: 'g', 132: 'c'}
cg4 = {0: 't', 33: 'g', 38: 'g', 44: 'g', 54: 'c', 56: 'a', 58: 'a', 66: 'c', 132: 'c'}

#remove last snp for shorter cg sequence --- note, also change varsInCG
del cg1[132]
del cg2[132]
del cg3[132]
del cg4[132]

#reference sequences are cut into smaller parts of 'chunklength' length, and with 'chunklength' / 2 overlap
chunklength = 8

#create the chunks of the reference sequence with regular expressions for the variable nucleotides
for i in range(0, len(searchstrings["ca"]) - chunklength, chunklength / 2):
  pos = i
  chunk = searchstrings["ca"][i:i+chunklength]
  result = ""
  varsInResult = 0
  for c in chunk:
    if pos in ca1.keys():
      varsInResult += 1
      result += "[" + ca1[pos] + ca2[pos] + "]"
    else:
      result += c
    pos += 1
  compiledregex["ca"].append((re.compile(result), varsInResult))

for i in range(0, len(searchstrings["cg"]) - chunklength, chunklength / 2):
  pos = i
  chunk = searchstrings["cg"][i:i+chunklength]
  result = ""
  varsInResult = 0
  for c in chunk:
    if pos in cg1.keys():
      varsInResult += 1
      result += "[" + "".join(set([cg1[pos], cg2[pos], cg3[pos], cg4[pos]])) + "]"
    else:
      result += c
    pos += 1
  compiledregex["cg"].append((re.compile(result), varsInResult))

for i in range(0, len(searchstrings["cm"]) - chunklength, chunklength / 2):
  compiledregex["cm"].append((re.compile(searchstrings["cm"][i:i+chunklength]), False))

for i in range(0, len(searchstrings["ce"]) - chunklength + 1, chunklength / 2):
  compiledregex["ce"].append((re.compile(searchstrings["ce"][i:i+chunklength]), False))

def removeAndReturnMaxIndex(x): #simplifies a list comprehension
  m = max(x)
  index = x.index(m)
  x[index] = 0
  return index
  

start_location = dict()
hits = dict()
alltotal = 0
for key in compiledregex.keys(): #for ca/cg/cm/ce
	regularexpressions = compiledregex[key] #get the compiled regular expressions
	for ID in dic.keys()[0:]: #for every ID
		if ID not in hits.keys(): #ensure that the dictionairy that keeps track of the hits for every gene exists
			hits[ID] = {"ca_hits": 0, "cg_hits": 0, "cm_hits": 0, "ce_hits": 0, "ca1": 0, "ca2": 0, "cg1": 0, "cg2": 0, "cg3": 0, "cg4": 0}
		currentIDHits = hits[ID]
		seq = dic[ID]
		lastindex = 0
		start_zero = len(searchstrings[key]) #allows the reference sequence to start before search sequence (start_locations of < 0)
		start = [0] * (len(seq) + start_zero)
		for i, regexp in enumerate(regularexpressions): #for every regular expression
			relativeStartLocation = lastindex - (chunklength / 2) * i
			if relativeStartLocation >= len(seq):
				break
			regex, hasVar = regexp
			matches = regex.finditer(seq[lastindex:])
			for match in matches: #for every match with the current regex, only uses the first hit because of the break at the end of this loop
				lastindex += match.start()
				start[relativeStartLocation + start_zero] += 1
				if hasVar: #if the regex has a variable nt in it
					chunkstart = chunklength / 2 * i #where in the reference does this chunk start
					chunkend = chunklength / 2 * i + chunklength #where in the reference does this chunk end
					if key == "ca": #just calculate the variable nt score for 'ca', cheaper
						currentIDHits["ca1"] += len([1 for x in ca1 if chunkstart <= x < chunkend and ca1[x] == seq[lastindex + x - chunkstart]])
						currentIDHits["ca2"] += len([1 for x in ca2 if chunkstart <= x < chunkend and ca2[x] == seq[lastindex + x - chunkstart]])
					elif key == "cg": #just calculate the variable nt score for 'cg', cheaper
						currentIDHits["cg1"] += len([1 for x in cg1 if chunkstart <= x < chunkend and cg1[x] == seq[lastindex + x - chunkstart]])
						currentIDHits["cg2"] += len([1 for x in cg2 if chunkstart <= x < chunkend and cg2[x] == seq[lastindex + x - chunkstart]])
						currentIDHits["cg3"] += len([1 for x in cg3 if chunkstart <= x < chunkend and cg3[x] == seq[lastindex + x - chunkstart]])
						currentIDHits["cg4"] += len([1 for x in cg4 if chunkstart <= x < chunkend and cg4[x] == seq[lastindex + x - chunkstart]])
					else: #key == "cm" #no variable regions in 'cm' or 'ce'
						pass
				break #this only breaks when there was a match with the regex, breaking means the 'else:' clause is skipped
			else: #only runs if there were no hits
				continue
			#print "found ", regex.pattern , "at", lastindex, "adding one to", (lastindex - chunklength / 2 * i), "to the start array of", ID, "gene", key, "it's now:", start[lastindex - chunklength / 2 * i]
			currentIDHits[key + "_hits"] += 1
		start_location[ID + "_" + key] = str([(removeAndReturnMaxIndex(start) + 1 - start_zero) for x in range(5) if len(start) > 0 and max(start) > 1])
		#start_location[ID + "_" + key] = str(start.index(max(start)))


varsInCA = float(len(ca1.keys()) * 2)
varsInCG = float(len(cg1.keys()) * 2) - 2 # -2 because the sliding window doesn't hit the first and last nt twice
varsInCM = 0
varsInCE = 0

def round_int(val):
	return int(round(val))

first = True
seq_write_count=0
with open(infile, 'r') as f: #read all sequences into a dictionary as key = ID, value = sequence
	with open(output, 'w') as o:
		for line in f:
			total += 1
			if first:
				o.write("Sequence ID\tbest_match\tnt_hit_percentage\tchunk_hit_percentage\tstart_locations\n")
				first = False
				continue
			linesplt = line.split("\t")
			if linesplt[2] == "No results":
				pass
			ID = linesplt[1]
			currentIDHits = hits[ID]
			possibleca = float(len(compiledregex["ca"]))
			possiblecg = float(len(compiledregex["cg"]))
			possiblecm = float(len(compiledregex["cm"]))
			possiblece = float(len(compiledregex["ce"]))
			cahits = currentIDHits["ca_hits"]
			cghits = currentIDHits["cg_hits"]
			cmhits = currentIDHits["cm_hits"]
			cehits = currentIDHits["ce_hits"]
			if cahits >= cghits and cahits >= cmhits and cahits >= cehits: #its a ca gene
				ca1hits = currentIDHits["ca1"]
				ca2hits = currentIDHits["ca2"]
				if ca1hits >= ca2hits:
					o.write(ID + "\tIGA1\t" + str(round_int(ca1hits / varsInCA * 100)) + "\t" + str(round_int(cahits / possibleca * 100)) + "\t" + start_location[ID + "_ca"] + "\n")
				else:
					o.write(ID + "\tIGA2\t" + str(round_int(ca2hits / varsInCA * 100)) + "\t" + str(round_int(cahits / possibleca * 100)) + "\t" + start_location[ID + "_ca"] + "\n")
			elif cghits >= cahits and cghits >= cmhits and cghits >= cehits: #its a cg gene
				cg1hits = currentIDHits["cg1"]
				cg2hits = currentIDHits["cg2"]
				cg3hits = currentIDHits["cg3"]
				cg4hits = currentIDHits["cg4"]
				if cg1hits >= cg2hits and cg1hits >= cg3hits and cg1hits >= cg4hits: #cg1 gene
					o.write(ID + "\tIGG1\t" + str(round_int(cg1hits / varsInCG * 100)) + "\t" + str(round_int(cghits / possiblecg * 100)) + "\t" + start_location[ID + "_cg"] + "\n")
				elif cg2hits >= cg1hits and cg2hits >= cg3hits and cg2hits >= cg4hits: #cg2 gene
					o.write(ID + "\tIGG2\t" + str(round_int(cg2hits / varsInCG * 100)) + "\t" + str(round_int(cghits / possiblecg * 100)) + "\t" + start_location[ID + "_cg"] + "\n")
				elif cg3hits >= cg1hits and cg3hits >= cg2hits and cg3hits >= cg4hits: #cg3 gene
					o.write(ID + "\tIGG3\t" + str(round_int(cg3hits / varsInCG * 100)) + "\t" + str(round_int(cghits / possiblecg * 100)) + "\t" + start_location[ID + "_cg"] + "\n")
				else: #cg4 gene
					o.write(ID + "\tIGG4\t" + str(round_int(cg4hits / varsInCG * 100)) + "\t" + str(round_int(cghits / possiblecg * 100)) + "\t" + start_location[ID + "_cg"] + "\n")
			else: #its a cm or ce gene
				if cmhits >= cehits:
					o.write(ID + "\tIGM\t100\t" + str(round_int(cmhits / possiblecm * 100)) + "\t" + start_location[ID + "_cm"] + "\n")
				else:
					o.write(ID + "\tIGE\t100\t" + str(round_int(cehits / possiblece * 100)) + "\t" + start_location[ID + "_ce"] + "\n")
			seq_write_count += 1

print "Time: %i" % (int(time.time() * 1000) - starttime)

print "Number of sequences written to file:", seq_write_count





