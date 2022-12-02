#!/usr/bin/env python3

"""Script to split IMGT file into several archives for each of the genes"""

import argparse
import io
import os
import tarfile
import tempfile
from typing import Iterator, Tuple


def merged_txt_to_match_dict(merged: str):
    with open(merged, "rt") as f:
        header = next(f).strip("\n")
        column_names = header.split("\t")
        best_match_index = column_names.index("best_match")
        sequence_id_index = column_names.index("Sequence.ID")
        match_dict = {}
        for line in f:
            values = line.strip().split("\t")
            sequence_id = values[sequence_id_index]
            best_match = values[best_match_index]
            if "unmatched" in best_match:
                # For some reason the table has values such as: unmatched, IGA2
                continue
            match_dict[sequence_id] = best_match
    return match_dict


def imgt_to_tables(imgt_file: str) -> Iterator[Tuple[str, io.TextIOWrapper]]:
    with tarfile.open(imgt_file, "r") as archive:
        while True:
            member = archive.next()
            if member is None:
                return
            if member.name in {"README.txt"}:
                continue
            if member.name.startswith("11_"):
                continue
            f = archive.extractfile(member)
            f_text = io.TextIOWrapper(f)
            yield member.name, f_text
            f_text.close()


def split_imgt(imgt_file, merged_file, outdir):
    match_dict = merged_txt_to_match_dict(merged_file)
    genes = ["", "IGA", "IGA1", "IGA2", "IGG", "IGG1", "IGG2", "IGG3", "IGG4",
             "IGM", "IGE"]
    gene_tarfiles = []
    os.makedirs(outdir, exist_ok=True)
    for gene in genes:
        new_filename = f"new_IMGT_{gene}.txz" if gene else "new_IMGT.txz"
        gene_tarfiles.append(
            tarfile.open(os.path.join(outdir, new_filename), mode="w:xz")
        )
    for name, table in imgt_to_tables(imgt_file):
        # [0] select the file descriptor for the open function
        gene_files = []
        for gene in genes:
            fp, fname = tempfile.mkstemp()
            f = open(fp, mode="wt")
            gene_files.append((gene, f, fname))
        header = next(table)
        column_names = header.strip("\n").split("\t")
        fr1_columns = [index for index, column in enumerate(column_names)
                       if column.startswith("FR1")]
        sequence_id_index = column_names.index("Sequence ID")
        for _, gene_file, _ in gene_files:
            gene_file.write(header)
        for line in table:
            values = line.strip("\n").split("\t")
            sequence_id = values[sequence_id_index]
            match = match_dict.get(sequence_id)
            if match is None:
                continue
            if name.startswith("8_"):
                # change the FR1 columns to 0 in the "8_..." file
                for index in fr1_columns:
                    values[index] = "0"
                line = "\t".join(values) + "\n"
            for gene, gene_file, _ in gene_files:
                if gene in match:
                    gene_file.write(line)
        for gene_tarfile, (_, gene_file, fname) in zip(gene_tarfiles, gene_files):
            gene_file.flush()
            gene_tarfile.add(fname, name)
            gene_file.close()
    for gene_tarfile in gene_tarfiles:
        gene_tarfile.close()


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("imgt_file", help="The original IMGT FILE")
    parser.add_argument("merged", help="merged.txt file")
    parser.add_argument("--outdir", help="output directory")
    return parser


def main():
    args = argument_parser().parse_args()
    split_imgt(args.imgt_file, args.merged, args.outdir)


if __name__ == "__main__":
    main()
