import argparse
import logging
import sys
import os
import typing
from typing import Optional

from collections import defaultdict


class Mutation(typing.NamedTuple):
	"""Represent a mutation type as a tuple"""
	frm: str  # 'from' is a reserved python keyword.
	where: int
	to: str
	frmAA: Optional[str] = None
	whereAA: Optional[int] = None
	toAA: Optional[str] = None
	thing: Optional[str] = None  # '(---)' or '(+-+)' etc. No idea

	@classmethod
	def from_string(cls, string: str):
		# Complete mutation example: a88>g,I30>V(+ - +)
		# Only nucleotide example: g303>t
		if ',' in string:
			nucleotide_change, aa_change = string.split(',', maxsplit=1)  # type: str, Optional[str]
		else:
			nucleotide_change = string
			aa_change = None
		frm_part, to = nucleotide_change.split('>', maxsplit=1)
		frm = frm_part[0]
		where = int(frm_part[1:])

		if aa_change is None:
			return cls(frm, where, to)

		frmAA_part, toAA_part = aa_change.split('>', maxsplit=1)  # type: str, str
		frmAA = frmAA_part[0]
		whereAA = int(frmAA_part[1:])
		brace_start = toAA_part.index('(')
		toAA = toAA_part[:brace_start]
		thing = toAA_part[brace_start:]
		return cls(frm, where, to, frmAA, whereAA, toAA, thing)


class Hotspot(typing.NamedTuple):
	start: int
	end: int
	region: str

	@classmethod
	def from_string(cls, string):
		# Example: aa,40-41(FR1)
		sequence, rest = string.split(',')  # type: str, str
		brace_pos = rest.index('(')
		numbers = rest[:brace_pos]
		start, end = numbers.split('-')
		region = rest[brace_pos + 1:-1]  # Remove the braces
		return cls(int(start), int(end), region)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", help="The '7_V-REGION-mutation-and-AA-change-table' and '10_V-REGION-mutation-hotspots' merged together, with an added 'best_match' annotation")
	parser.add_argument("--genes", help="The genes available in the 'best_match' column")
	parser.add_argument("--empty_region_filter", help="Where does the sequence start?", choices=['leader', 'FR1', 'CDR1', 'FR2'])
	parser.add_argument("--output", help="Output file")

	args = parser.parse_args()

	infile = args.input
	genes = str(args.genes).split(",")
	empty_region_filter = args.empty_region_filter
	outfile = args.output

	genedic = dict()

	mutationdic = dict()
	NAMatchResult = (None, None, None, None, None, None, '')
	linecount = 0

	IDIndex = 0
	best_matchIndex = 0
	fr1Index = 0
	cdr1Index = 0
	fr2Index = 0
	cdr2Index = 0
	fr3Index = 0
	first = True
	IDlist = []
	mutationList = []
	mutationListByID = {}
	cdr1LengthDic = {}
	cdr2LengthDic = {}

	fr1LengthDict = {}
	fr2LengthDict = {}
	fr3LengthDict = {}

	cdr1LengthIndex = 0
	cdr2LengthIndex = 0

	fr1SeqIndex = 0
	fr2SeqIndex = 0
	fr3SeqIndex = 0

	tandem_sum_by_class = defaultdict(int)
	expected_tandem_sum_by_class = defaultdict(float)

	with open(infile, 'r') as i:
		for line in i:
			if first:
				linesplt = line.split("\t")
				IDIndex = linesplt.index("Sequence.ID")
				best_matchIndex = linesplt.index("best_match")
				fr1Index = linesplt.index("FR1.IMGT")
				cdr1Index = linesplt.index("CDR1.IMGT")
				fr2Index = linesplt.index("FR2.IMGT")
				cdr2Index = linesplt.index("CDR2.IMGT")
				fr3Index = linesplt.index("FR3.IMGT")
				cdr1LengthIndex = linesplt.index("CDR1.IMGT.length")
				cdr2LengthIndex = linesplt.index("CDR2.IMGT.length")
				fr1SeqIndex = linesplt.index("FR1.IMGT.seq")
				fr2SeqIndex = linesplt.index("FR2.IMGT.seq")
				fr3SeqIndex = linesplt.index("FR3.IMGT.seq")
				first = False
				continue
			linecount += 1
			linesplt = line.split("\t")
			ID = linesplt[IDIndex]
			genedic[ID] = linesplt[best_matchIndex]
			
			mutationdic[ID + "_FR1"] = []
			if len(linesplt[fr1Index]) > 5 and empty_region_filter == "leader":
				mutationdic[ID + "_FR1"] = [Mutation.from_string(x) for x in linesplt[fr1Index].split("|") if x]

			mutationdic[ID + "_CDR1"] = []
			if len(linesplt[cdr1Index]) > 5 and empty_region_filter in ["leader", "FR1"]:
				mutationdic[ID + "_CDR1"] = [Mutation.from_string(x) for x in linesplt[cdr1Index].split("|") if x]

			mutationdic[ID + "_FR2"] = []
			if len(linesplt[fr2Index]) > 5 and empty_region_filter in ["leader", "FR1", "CDR1"]:
				mutationdic[ID + "_FR2"] = [Mutation.from_string(x) for x in linesplt[fr2Index].split("|") if x]

			mutationdic[ID + "_CDR2"] = []
			if len(linesplt[cdr2Index]) > 5:
				mutationdic[ID + "_CDR2"] = [Mutation.from_string(x) for x in linesplt[cdr2Index].split("|") if x]
			
			mutationdic[ID + "_FR2-CDR2"] = mutationdic[ID + "_FR2"] + mutationdic[ID + "_CDR2"]

			mutationdic[ID + "_FR3"] = []
			if len(linesplt[fr3Index]) > 5:
				mutationdic[ID + "_FR3"] = [Mutation.from_string(x) for x in linesplt[fr3Index].split("|") if x]
				
			mutationList += mutationdic[ID + "_FR1"] + mutationdic[ID + "_CDR1"] + mutationdic[ID + "_FR2"] + mutationdic[ID + "_CDR2"] + mutationdic[ID + "_FR3"]
			mutationListByID[ID] = mutationdic[ID + "_FR1"] + mutationdic[ID + "_CDR1"] + mutationdic[ID + "_FR2"] + mutationdic[ID + "_CDR2"] + mutationdic[ID + "_FR3"]

			try:
				cdr1Length = int(linesplt[cdr1LengthIndex])
			except:
				cdr1Length = 0
			
			try:
				cdr2Length = int(linesplt[cdr2LengthIndex])
			except:
				cdr2Length = 0

			#print linesplt[fr2SeqIndex]
			fr1Length = len(linesplt[fr1SeqIndex]) if empty_region_filter == "leader" else 0
			fr2Length = len(linesplt[fr2SeqIndex]) if empty_region_filter in ["leader", "FR1", "CDR1"] else 0
			fr3Length = len(linesplt[fr3SeqIndex])

			cdr1LengthDic[ID] = cdr1Length
			cdr2LengthDic[ID] = cdr2Length

			fr1LengthDict[ID] = fr1Length
			fr2LengthDict[ID] = fr2Length
			fr3LengthDict[ID] = fr3Length

			IDlist += [ID]
	print("len(mutationdic) =", len(mutationdic))

	with open(os.path.join(os.path.dirname(os.path.abspath(infile)), "mutationdict.txt"), 'w') as out_handle:
		for ID, lst in mutationdic.items():
			for mut in lst:
				out_handle.write("{0}\t{1}\n".format(ID, "\t".join([str(x) for x in mut])))

	#tandem mutation stuff
	tandem_frequency = defaultdict(int)
	mutation_frequency = defaultdict(int)
	
	mutations_by_id_dic = {}
	first = True
	mutation_by_id_file = os.path.join(os.path.dirname(outfile), "mutation_by_id.txt")
	with open(mutation_by_id_file, 'r') as mutation_by_id:
		for l in mutation_by_id:
			if first:
				first = False
				continue
			splt = l.split("\t")
			mutations_by_id_dic[splt[0]] = int(splt[1])
    
	tandem_file = os.path.join(os.path.dirname(outfile), "tandems_by_id.txt")
	with open(tandem_file, 'w') as o:
		highest_tandem_length = 0

		o.write("Sequence.ID\tnumber_of_mutations\tnumber_of_tandems\tregion_length\texpected_tandems\tlongest_tandem\ttandems\n")
		for ID in IDlist:
			mutations = mutationListByID[ID]
			if len(mutations) == 0:
				continue
			last_mut = max(mutations, key=lambda x: int(x[1]))

			last_mut_pos = int(last_mut[1])

			mut_positions = [False] * (last_mut_pos + 1)

			for mutation in mutations:
				frm, where, to, frmAA, whereAA, toAA, thing = mutation
				where = int(where)
				mut_positions[where] = True

			tandem_muts = []
			tandem_start = -1
			tandem_length = 0
			for i in range(len(mut_positions)):
				if mut_positions[i]:
					if tandem_start == -1:
						tandem_start = i
					tandem_length += 1
					#print "".join(["1" if x else "0" for x in mut_positions[:i+1]])
				else:
					if tandem_length > 1:
						tandem_muts.append((tandem_start, tandem_length))
						#print "{0}{1} {2}:{3}".format(" " * (i - tandem_length), "^" * tandem_length, tandem_start, tandem_length)
					tandem_start = -1
					tandem_length = 0
			if tandem_length > 1:  # if the sequence ends with a tandem mutation
				tandem_muts.append((tandem_start, tandem_length))

			if len(tandem_muts) > 0:
				if highest_tandem_length < len(tandem_muts):
					highest_tandem_length = len(tandem_muts)

			region_length = fr1LengthDict[ID] + cdr1LengthDic[ID] + fr2LengthDict[ID] + cdr2LengthDic[ID] + fr3LengthDict[ID]
			longest_tandem = max(tandem_muts, key=lambda x: x[1]) if len(tandem_muts) else (0, 0)
			num_mutations = mutations_by_id_dic[ID] # len(mutations)
			f_num_mutations = float(num_mutations)
			num_tandem_muts = len(tandem_muts)
			expected_tandem_muts = f_num_mutations * (f_num_mutations - 1.0) / float(region_length)
			o.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(ID,
																str(num_mutations),
																str(num_tandem_muts),
																str(region_length),
																str(round(expected_tandem_muts, 2)),
																str(longest_tandem[1]),
																str(tandem_muts)))
			gene = genedic[ID]
			if gene.find("unmatched") == -1:
				tandem_sum_by_class[gene] += num_tandem_muts
				expected_tandem_sum_by_class[gene] += expected_tandem_muts

				tandem_sum_by_class["all"] += num_tandem_muts
				expected_tandem_sum_by_class["all"] += expected_tandem_muts

				gene = gene[:3]
				if gene in ["IGA", "IGG"]:
					tandem_sum_by_class[gene] += num_tandem_muts
					expected_tandem_sum_by_class[gene] += expected_tandem_muts
			else:
				tandem_sum_by_class["unmatched"] += num_tandem_muts
				expected_tandem_sum_by_class["unmatched"] += expected_tandem_muts


			for tandem_mut in tandem_muts:
				tandem_frequency[str(tandem_mut[1])] += 1
			#print "\t".join([ID, str(len(tandem_muts)), str(longest_tandem[1]) , str(tandem_muts)])

	tandem_freq_file = os.path.join(os.path.dirname(outfile), "tandem_frequency.txt")
	with open(tandem_freq_file, 'w') as o:
		for frq in sorted([int(x) for x in list(tandem_frequency.keys())]):
			o.write("{0}\t{1}\n".format(frq, tandem_frequency[str(frq)]))

	tandem_row = []
	genes_extra = list(genes)
	genes_extra.append("all")
	for x, y, in zip([tandem_sum_by_class[x] for x in genes_extra], [expected_tandem_sum_by_class[x] for x in genes_extra]):
		if y != 0:
			tandem_row += [x, round(y, 2), round(x / y, 2)]
		else:
			tandem_row += [x, round(y, 2), 0]

	tandem_freq_file = os.path.join(os.path.dirname(outfile), "shm_overview_tandem_row.txt")
	with open(tandem_freq_file, 'w') as o:
		o.write("Tandems/Expected (ratio),{0}\n".format(",".join([str(x) for x in tandem_row])))

	#print mutationList, linecount

	AALength = (int(max(mutationList, key=lambda i: int(i[4]) if i[4] and i[5] != ";" else 0)[4]) + 1)  # [4] is the position of the AA mutation, None if silent
	if AALength < 60:
		AALength = 64

	AA_mutation = [0] * AALength
	AA_mutation_dic = {"IGA": AA_mutation[:], "IGG": AA_mutation[:], "IGM": AA_mutation[:], "IGE": AA_mutation[:], "unm": AA_mutation[:], "all": AA_mutation[:]}
	AA_mutation_empty = AA_mutation[:]

	print("AALength:", AALength)
	aa_mutations_by_id_file = outfile[:outfile.rindex("/")] + "/aa_id_mutations.txt"
	with open(aa_mutations_by_id_file, 'w') as o:
		o.write("ID\tbest_match\t" + "\t".join([str(x) for x in range(1,AALength)]) + "\n")
		for ID in list(mutationListByID.keys()):
			AA_mutation_for_ID = AA_mutation_empty[:]
			for mutation in mutationListByID[ID]:
				if mutation[4] and mutation[5] != ";":
					AA_mutation_position = int(mutation[4])
					try:
						AA_mutation[AA_mutation_position] += 1
						AA_mutation_for_ID[AA_mutation_position] += 1
					except Exception as e:
						print(e)
						print(mutation)
						sys.exit()
					clss = genedic[ID][:3]
					AA_mutation_dic[clss][AA_mutation_position] += 1
			o.write(ID + "\t" + genedic[ID] + "\t" + "\t".join([str(x) for x in AA_mutation_for_ID[1:]]) + "\n")



	#absent AA stuff
	absentAACDR1Dic = defaultdict(list)
	absentAACDR1Dic[5] = list(range(29,36))
	absentAACDR1Dic[6] = list(range(29,35))
	absentAACDR1Dic[7] = list(range(30,35))
	absentAACDR1Dic[8] = list(range(30,34))
	absentAACDR1Dic[9] = list(range(31,34))
	absentAACDR1Dic[10] = list(range(31,33))
	absentAACDR1Dic[11] = [32]

	absentAACDR2Dic = defaultdict(list)
	absentAACDR2Dic[0] = list(range(55,65))
	absentAACDR2Dic[1] = list(range(56,65))
	absentAACDR2Dic[2] = list(range(56,64))
	absentAACDR2Dic[3] = list(range(57,64))
	absentAACDR2Dic[4] = list(range(57,63))
	absentAACDR2Dic[5] = list(range(58,63))
	absentAACDR2Dic[6] = list(range(58,62))
	absentAACDR2Dic[7] = list(range(59,62))
	absentAACDR2Dic[8] = list(range(59,61))
	absentAACDR2Dic[9] = [60]

	absentAA = [len(IDlist)] * (AALength-1)
	for k, cdr1Length in cdr1LengthDic.items():
		for c in absentAACDR1Dic[cdr1Length]:
			absentAA[c] -= 1

	for k, cdr2Length in cdr2LengthDic.items():
		for c in absentAACDR2Dic[cdr2Length]:
			absentAA[c] -= 1


	aa_mutations_by_id_file = outfile[:outfile.rindex("/")] + "/absent_aa_id.txt"
	with open(aa_mutations_by_id_file, 'w') as o:
		o.write("ID\tcdr1length\tcdr2length\tbest_match\t" + "\t".join([str(x) for x in range(1,AALength)]) + "\n")
		for ID in IDlist:
			absentAAbyID = [1] * (AALength-1)
			cdr1Length = cdr1LengthDic[ID]
			for c in absentAACDR1Dic[cdr1Length]:
				absentAAbyID[c] -= 1

			cdr2Length = cdr2LengthDic[ID]
			for c in absentAACDR2Dic[cdr2Length]:
				absentAAbyID[c] -= 1
			o.write(ID + "\t" + str(cdr1Length) + "\t" + str(cdr2Length) + "\t" + genedic[ID] + "\t" + "\t".join([str(x) for x in absentAAbyID]) + "\n")

	if linecount == 0:
		print("No data, exiting")
		with open(outfile, 'w') as o:
			o.write("RGYW (%)," + ("0,0,0\n" * len(genes)))
			o.write("WRCY (%)," + ("0,0,0\n" * len(genes)))
			o.write("WA (%)," + ("0,0,0\n" * len(genes)))
			o.write("TW (%)," + ("0,0,0\n" * len(genes)))
		sys.exit()

	RGYWCount = {}
	WRCYCount = {}
	WACount = {}
	TWCount = {}

	#IDIndex = 0
	ataIndex = 0
	tatIndex = 0
	aggctatIndex = 0
	atagcctIndex = 0
	first = True
	with open(infile, 'r') as i:
		for line in i:
			if first:
				linesplt = line.split("\t")
				ataIndex = linesplt.index("X.a.t.a")
				tatIndex = linesplt.index("t.a.t.")
				aggctatIndex = linesplt.index("X.a.g.g.c.t..a.t.")
				atagcctIndex = linesplt.index("X.a.t..a.g.c.c.t.")
				first = False
				continue
			linesplt = line.split("\t")
			gene = linesplt[best_matchIndex]
			ID = linesplt[IDIndex]
			RGYW = [Hotspot.from_string(x) for x in linesplt[aggctatIndex].split("|") if x]
			WRCY = [Hotspot.from_string(x) for x in linesplt[atagcctIndex].split("|") if x]
			WA = [Hotspot.from_string(x) for x in linesplt[ataIndex].split("|") if x]
			TW = [Hotspot.from_string(x) for x in linesplt[tatIndex].split("|") if x]
			RGYWCount[ID], WRCYCount[ID], WACount[ID], TWCount[ID] = 0, 0, 0, 0

			with open(os.path.join(os.path.dirname(os.path.abspath(infile)), "RGYW.txt"), 'a') as out_handle:
				for hotspot in RGYW:
					out_handle.write("{0}\t{1}\n".format(ID, "\t".join([str(x) for x in hotspot])))

			mutationList = mutationdic[ID + "_FR1"] + mutationdic[ID + "_CDR1"] + mutationdic[ID + "_FR2"] + mutationdic[ID + "_CDR2"] + mutationdic[ID + "_FR3"]
			for mutation in mutationList:
				frm, where, to, AAfrm, AAwhere, AAto, junk = mutation
				mutation_in_RGYW = any(((start <= int(where) <= end) for (start, end, region) in RGYW))
				mutation_in_WRCY = any(((start <= int(where) <= end) for (start, end, region) in WRCY))
				mutation_in_WA = any(((start <= int(where) <= end) for (start, end, region) in WA))
				mutation_in_TW = any(((start <= int(where) <= end) for (start, end, region) in TW))

				in_how_many_motifs = sum([mutation_in_RGYW, mutation_in_WRCY, mutation_in_WA, mutation_in_TW])

				if in_how_many_motifs > 0:
					RGYWCount[ID] += (1.0 * int(mutation_in_RGYW)) / in_how_many_motifs
					WRCYCount[ID] += (1.0 * int(mutation_in_WRCY)) / in_how_many_motifs
					WACount[ID] += (1.0 * int(mutation_in_WA)) / in_how_many_motifs
					TWCount[ID] += (1.0 * int(mutation_in_TW)) / in_how_many_motifs
			
			mutations_in_motifs_file = os.path.join(os.path.dirname(os.path.abspath(infile)), "mutation_in_motifs.txt")
			if not os.path.exists(mutation_by_id_file):
				with open(mutations_in_motifs_file, 'w') as out_handle:
					out_handle.write("{0}\n".format("\t".join([
						"Sequence.ID",
						"mutation_position",
						"region",
						"from_nt",
						"to_nt",
						"mutation_position_AA",
						"from_AA",
						"to_AA",
						"motif",
						"motif_start_nt",
						"motif_end_nt",
						"rest"
					])))

			with open(mutations_in_motifs_file, 'a') as out_handle:
				motif_dic = {"RGYW": RGYW, "WRCY": WRCY, "WA": WA, "TW": TW}
				for mutation in mutationList:
					frm, where, to, AAfrm, AAwhere, AAto, junk = mutation
					for motif in list(motif_dic.keys()):
							
						for start, end, region in motif_dic[motif]:
							if start <= int(where) <= end:
								out_handle.write("{0}\n".format(
									"\t".join([
										ID,
										str(where),
										region,
										frm,
										to,
										str(AAwhere),
										str(AAfrm),
										str(AAto),
										motif,
										str(start),
										str(end),
										str(junk)
									])
								))



	def mean(lst):
		return (float(sum(lst)) / len(lst)) if len(lst) > 0 else 0.0


	def median(lst):
		lst = sorted(lst)
		l = len(lst)
		if l == 0:
			return 0
		if l == 1:
			return lst[0]
			
		l = int(l / 2)
		
		if len(lst) % 2 == 0:
			return float(lst[l] + lst[(l - 1)]) / 2.0
		else:
			return lst[l]

	funcs = {"mean": mean, "median": median, "sum": sum}

	directory = outfile[:outfile.rfind("/") + 1]
	value = 0
	valuedic = dict()

	for fname in list(funcs.keys()):
		for gene in genes:
			with open(directory + gene + "_" + fname + "_value.txt", 'r') as v:
				valuedic[gene + "_" + fname] = float(v.readlines()[0].rstrip())
		with open(directory + "all_" + fname + "_value.txt", 'r') as v:
			valuedic["total_" + fname] = float(v.readlines()[0].rstrip())
		

	def get_xyz(lst, gene, f, fname):
		x = round(round(f(lst), 1))
		y = valuedic[gene + "_" + fname]
		z = str(round(x / float(y) * 100, 1)) if y != 0 else "0"
		return (str(x), str(y), z)

	dic = {"RGYW": RGYWCount, "WRCY": WRCYCount, "WA": WACount, "TW": TWCount}
	arr = ["RGYW", "WRCY", "WA", "TW"]

	for fname in list(funcs.keys()):
		func = funcs[fname]
		foutfile = outfile[:outfile.rindex("/")] + "/hotspot_analysis_" + fname + ".txt"
		with open(foutfile, 'w') as o:
			for typ in arr:
				o.write(typ + " (%)")
				curr = dic[typ]
				for gene in genes:
					if valuedic[gene + "_" + fname] is 0:
						o.write(",0,0,0")
					else:
						x, y, z = get_xyz([curr[x] for x in [y for y, z in genedic.items() if z.startswith(gene)]], gene, func, fname)
						o.write("," + x + "," + y + "," + z)
				x, y, z = get_xyz([y for x, y in curr.items() if not genedic[x].startswith("unmatched")], "total", func, fname)
				#x, y, z = get_xyz([y for x, y in curr.iteritems()], "total", func, fname)
				o.write("," + x + "," + y + "," + z + "\n")


	# for testing
	seq_motif_file = outfile[:outfile.rindex("/")] + "/motif_per_seq.txt"
	with open(seq_motif_file, 'w') as o:
		o.write("ID\tRGYW\tWRCY\tWA\tTW\n")
		for ID in IDlist:
			#o.write(ID + "\t" + str(round(RGYWCount[ID], 2)) + "\t" + str(round(WRCYCount[ID], 2)) + "\t" + str(round(WACount[ID], 2)) + "\t" + str(round(TWCount[ID], 2)) + "\n")
			o.write(ID + "\t" + str(RGYWCount[ID]) + "\t" + str(WRCYCount[ID]) + "\t" + str(WACount[ID]) + "\t" + str(TWCount[ID]) + "\n")

if __name__ == "__main__":
	main()
