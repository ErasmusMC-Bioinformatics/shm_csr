#!/usr/bin/env python3

"""
Find naive mutations (< 2% mutated) for IGM genes
"""

import argparse


def find_naive_mutations(mutation_file, naive_file, naive_memory_file,
                         percentage_cutoff=0.02):
    with (open(mutation_file, "rt") as mutations,
          open(naive_file, "wt") as naive,
          open(naive_memory_file) as naive_memory,):
        for line in mutations:
            sequence_id, best_match, mutation_no, region_length, _ = \
                line.strip('\n').split('\t')
            if best_match != "IGM":
                continue
            mutation_no = int(mutation_no)
            region_length = int(region_length)
            if (mutation_no / region_length) < percentage_cutoff:
                naive.write(line)
            else:
                naive_memory.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mutation_file", help="scatter.txt")
    parser.add_argument("naive_file")
    parser.add_argument("naive_memory_file")


if __name__ == "__main__":
    main()