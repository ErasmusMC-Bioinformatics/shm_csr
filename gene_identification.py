#!/usr/bin/env python3

import argparse
import re
from typing import Dict, Iterator, List, Tuple


def generate_sequence_and_id_from_summary(summary_file: str
                                          ) -> Iterator[Tuple[str, str]]:
    with open(summary_file, "rt") as summary:
        header = next(summary)
        column_names = header.strip("\n").split("\t")
        id_column = column_names.index("Sequence ID")
        sequence_column = column_names.index("Sequence")
        for line in summary:
            values = line.strip("\n").split("\t")
            id = values[id_column]
            try:
                sequence = values[sequence_column]
            except IndexError:  # weird rows without a sequence
                sequence = ""
            yield id, sequence


# old cm sequence: gggagtgcatccgccccaacccttttccccctcgtctcctgtgagaattccc
# old cg sequence: ctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctg
# ggggcacagcggccctgggctgcctggtcaaggactacttccccgaaccggtgacggtgtcgtggaactcagg
# cgccctgaccag
SEARCHSTRINGS = {"ca": "catccccgaccagccccaaggtcttcccgctgagcctctgcagcacccagccag"
                       "atgggaacgtggtcatcgcctgcctgg",
                 "cg": "ctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctc"
                       "tgggggcacagcggcc",
                 "ce": "gcctccacacagagcccatccgtcttccccttgacccgctgctgcaaaaacatt"
                       "ccctcc",
                 "cm": "gggagtgcatccgccccaacc"} #new (shorter) cm sequence

#lambda/kappa referesearchstringsnce sequence variable nucleotides
CA1_MUTATIONS = {38: 't', 39: 'g', 48: 'a', 49: 'g', 51: 'c', 68: 'a', 73: 'c'}
CA2_MUTATIONS = {38: 'g', 39: 'a', 48: 'c', 49: 'c', 51: 'a', 68: 'g', 73: 'a'}
CG1_MUTATIONS = {0: 'c', 33: 'a', 38: 'c', 44: 'a', 54: 't', 56: 'g', 58: 'g', 66: 'g', 132: 'c'}
CG2_MUTATIONS = {0: 'c', 33: 'g', 38: 'g', 44: 'g', 54: 'c', 56: 'a', 58: 'a', 66: 'g', 132: 't'}
CG3_MUTATIONS = {0: 't', 33: 'g', 38: 'g', 44: 'g', 54: 't', 56: 'g', 58: 'g', 66: 'g', 132: 'c'}
CG4_MUTATIONS = {0: 't', 33: 'g', 38: 'g', 44: 'g', 54: 'c', 56: 'a', 58: 'a', 66: 'c', 132: 'c'}

#remove last snp for shorter cg sequence --- note, also change varsInCG
del CG1_MUTATIONS[132]
del CG2_MUTATIONS[132]
del CG3_MUTATIONS[132]
del CG4_MUTATIONS[132]

# reference sequences are cut into smaller parts of 'chunklength' length,
# and with 'chunklength' / 2 overlap
CHUNK_LENGTH = 8


def create_compiled_regexes() -> Dict[str, List[Tuple[re.Pattern, int]]]:

    compiledregex: Dict[str, List[Tuple[re.Pattern, int]]] = {
        "ca": [],
        "cg": [],
        "ce": [],
        "cm": []
    }

    for i in range(0, len(SEARCHSTRINGS["ca"]) - CHUNK_LENGTH, CHUNK_LENGTH // 2):
      pos = i
      chunk = SEARCHSTRINGS["ca"][i:i + CHUNK_LENGTH]
      result = ""
      varsInResult = 0
      for c in chunk:
        if pos in list(CA1_MUTATIONS.keys()):
          varsInResult += 1
          result += "[" + CA1_MUTATIONS[pos] + CA2_MUTATIONS[pos] + "]"
        else:
          result += c
        pos += 1
      compiledregex["ca"].append((re.compile(result), varsInResult))

    for i in range(0, len(SEARCHSTRINGS["cg"]) - CHUNK_LENGTH, CHUNK_LENGTH // 2):
      pos = i
      chunk = SEARCHSTRINGS["cg"][i:i + CHUNK_LENGTH]
      result = ""
      varsInResult = 0
      for c in chunk:
        if pos in list(CG1_MUTATIONS.keys()):
          varsInResult += 1
          result += "[" + "".join(set([CG1_MUTATIONS[pos], CG2_MUTATIONS[pos], CG3_MUTATIONS[pos], CG4_MUTATIONS[pos]])) + "]"
        else:
          result += c
        pos += 1
      compiledregex["cg"].append((re.compile(result), varsInResult))

    for i in range(0, len(SEARCHSTRINGS["cm"]) - CHUNK_LENGTH, CHUNK_LENGTH // 2):
      compiledregex["cm"].append((re.compile(SEARCHSTRINGS["cm"][i:i + CHUNK_LENGTH]), 0))

    for i in range(0, len(SEARCHSTRINGS["ce"]) - CHUNK_LENGTH + 1, CHUNK_LENGTH // 2):
      compiledregex["ce"].append((re.compile(SEARCHSTRINGS["ce"][i:i + CHUNK_LENGTH]), 0))

    return compiledregex


def removeAndReturnMaxIndex(x): #simplifies a list comprehension
  m = max(x)
  index = x.index(m)
  x[index] = 0
  return index


def match_sequence(seq, compiledregex):
    currentIDHits = {"ca_hits": 0, "cg_hits": 0, "cm_hits": 0, "ce_hits": 0,
            "ca1": 0, "ca2": 0, "cg1": 0, "cg2": 0, "cg3": 0, "cg4": 0}
    alltotal = 0
    start_location = dict()
    for key in compiledregex:  # for ca/cg/cm/ce
        regularexpressions = compiledregex[key]
        lastindex = 0
        start_zero = len(SEARCHSTRINGS[key]) #allows the reference sequence to start before search sequence (start_locations of < 0)
        start = [0] * (len(seq) + start_zero)
        for i, regexp in enumerate(regularexpressions): #for every regular expression
            relativeStartLocation = lastindex - (CHUNK_LENGTH // 2) * i
            if relativeStartLocation >= len(seq):
                break
            regex, hasVar = regexp
            matches = regex.finditer(seq[lastindex:])
            for match in matches: #for every match with the current regex, only uses the first hit because of the break at the end of this loop
                lastindex += match.start()
                start[relativeStartLocation + start_zero] += 1
                if hasVar: #if the regex has a variable nt in it
                    chunkstart = CHUNK_LENGTH // 2 * i #where in the reference does this chunk start
                    chunkend = CHUNK_LENGTH // 2 * i + CHUNK_LENGTH #where in the reference does this chunk end
                    if key == "ca": #just calculate the variable nt score for 'ca', cheaper
                        currentIDHits["ca1"] += len([1 for x in CA1_MUTATIONS if chunkstart <= x < chunkend and CA1_MUTATIONS[x] == seq[lastindex + x - chunkstart]])
                        currentIDHits["ca2"] += len([1 for x in CA2_MUTATIONS if chunkstart <= x < chunkend and CA2_MUTATIONS[x] == seq[lastindex + x - chunkstart]])
                    elif key == "cg": #just calculate the variable nt score for 'cg', cheaper
                        currentIDHits["cg1"] += len([1 for x in CG1_MUTATIONS if chunkstart <= x < chunkend and CG1_MUTATIONS[x] == seq[lastindex + x - chunkstart]])
                        currentIDHits["cg2"] += len([1 for x in CG2_MUTATIONS if chunkstart <= x < chunkend and CG2_MUTATIONS[x] == seq[lastindex + x - chunkstart]])
                        currentIDHits["cg3"] += len([1 for x in CG3_MUTATIONS if chunkstart <= x < chunkend and CG3_MUTATIONS[x] == seq[lastindex + x - chunkstart]])
                        currentIDHits["cg4"] += len([1 for x in CG4_MUTATIONS if chunkstart <= x < chunkend and CG4_MUTATIONS[x] == seq[lastindex + x - chunkstart]])
                    else: #key == "cm" #no variable regions in 'cm' or 'ce'
                        pass
                break #this only breaks when there was a match with the regex, breaking means the 'else:' clause is skipped
            else: #only runs if there were no hits
                continue
            #print "found ", regex.pattern , "at", lastindex, "adding one to", (lastindex - chunklength / 2 * i), "to the start array of", ID, "gene", key, "it's now:", start[lastindex - chunklength / 2 * i]
            currentIDHits[key + "_hits"] += 1
        start_location[key] = str([(removeAndReturnMaxIndex(start) + 1 - start_zero) for x in range(5) if len(start) > 0 and max(start) > 1])

    cahits = currentIDHits["ca_hits"]
    cghits = currentIDHits["cg_hits"]
    cmhits = currentIDHits["cm_hits"]
    cehits = currentIDHits["ce_hits"]
    if cahits >= cghits and cahits >= cmhits and cahits >= cehits:  # its a ca gene
        ca1hits = currentIDHits["ca1"]
        ca2hits = currentIDHits["ca2"]
        if ca1hits >= ca2hits:
            return "IGA1", ca1hits, cahits, start_location["ca"]
        else:
            return "IGA2", ca2hits, cahits, start_location["ca"]
    elif cghits >= cahits and cghits >= cmhits and cghits >= cehits:  # its a cg gene
        cg1hits = currentIDHits["cg1"]
        cg2hits = currentIDHits["cg2"]
        cg3hits = currentIDHits["cg3"]
        cg4hits = currentIDHits["cg4"]
        if cg1hits >= cg2hits and cg1hits >= cg3hits and cg1hits >= cg4hits:  # cg1 gene
            return "IGG1", cg1hits, cghits, start_location["cg"]
        elif cg2hits >= cg1hits and cg2hits >= cg3hits and cg2hits >= cg4hits:  # cg2 gene
            return "IGG2", cg2hits, cghits, start_location["cg"]
        elif cg3hits >= cg1hits and cg3hits >= cg2hits and cg3hits >= cg4hits:  # cg3 gene
            return "IGG3", cg3hits, cghits, start_location["cg"]
        else:  # cg4 gene
            return "IGG4", cg4hits, cghits, start_location["cg"]
    else:  # its a cm or ce gene
        if cmhits >= cehits:
            return "IGM", 0, cmhits, start_location["cm"]
        else:
            return "IGE", 0, cehits, start_location["ce"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        help="The 1_Summary file from an IMGT zip file")
    parser.add_argument("--output",
                        help="The annotated output file to be merged back "
                             "with the summary file")
    args = parser.parse_args()
    varsInCA = float(len(list(CA1_MUTATIONS.keys())) * 2)
    varsInCG = float(len(list(
        CG1_MUTATIONS.keys())) * 2) - 2  # -2 because the sliding window doesn't hit the first and last nt twice
    subclass_vars = {
        "IGA1": varsInCA, "IGA2": varsInCA,
        "IGG1": varsInCG, "IGG2": varsInCG, "IGG3": varsInCG, "IGG4": varsInCG,
        "IGE": 0,
        "IGM": 0,
    }
    compiledregex = create_compiled_regexes()
    possibleca = float(len(compiledregex["ca"]))
    possiblecg = float(len(compiledregex["cg"]))
    possiblecm = float(len(compiledregex["cm"]))
    possiblece = float(len(compiledregex["ce"]))
    class_chunks = {
        "IGA1": possibleca, "IGA2": possibleca,
        "IGE": possiblece,
        "IGG1": possiblecg, "IGG2": possiblecg, "IGG3": possiblecg,
        "IGG4": possiblecg,
        "IGM": possiblecm
    }
    with open(args.output, "wt") as output:
        output.write("Sequence ID\tbest_match\tnt_hit_percentage\t"
                     "chunk_hit_percentage\tstart_locations\n")
        for id, sequence in generate_sequence_and_id_from_summary(args.input):
            best_match, subclass_hits, class_hits, start_locations = \
                match_sequence(sequence, compiledregex)
            variable_nucs = subclass_vars[best_match]
            if variable_nucs:
                subclass_percentage = round(subclass_hits * 100 /
                                            variable_nucs)
            else:
                subclass_percentage = 100
            class_percentage = round(class_hits * 100 / class_chunks[best_match])
            output.write(f"{id}\t{best_match}\t{subclass_percentage}\t"
                         f"{class_percentage}\t{start_locations}\n")


if __name__ == "__main__":
    main()
