#!/usr/bin/env python3

"""
Find naive mutations (< 2% mutated) for IGM genes
"""

import argparse
import contextlib


def find_naive_mutations(mutation_file, naive_file, naive_memory_file,
                         percentage_cutoff=0.02):
    # A compound with statement throws a syntax error with the included python
    # 3.7.1 in the container, so use an exit stack instead.
    with contextlib.ExitStack() as stack:
        mutations = stack.enter_context(open(mutation_file, "rt"))
        naive = stack.enter_context(open(naive_file, "wt"))
        naive_memory = stack.enter_context(open(naive_memory_file, "wt"))
        header = next(mutations)
        naive.write(header)
        naive_memory.write(header)
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
    args = parser.parse_args()
    find_naive_mutations(args.mutation_file, args.naive_file,
                         args.naive_memory_file)


if __name__ == "__main__":
    main()