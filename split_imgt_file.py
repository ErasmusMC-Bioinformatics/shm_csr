#!/usr/bin/env python3

"""
Script to split IMGT file into several archives for each of the genes

Rather than creating each new archive individually this script will read
the input files only once and as such enormously shorten processing time.
"""

import argparse
import io
import os
import tarfile
from typing import Iterator, List, Tuple


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


def split_imgt(imgt_file: str, merged_file: str, outdir: str, genes: List[str],
               prefix: str):
    """
    This function creates a separate tar file for each of the gene matches
    based on the merged file. Unmatched genes are left out.
    :param imgt_file: The original IMGT file
    :param merged_file: The merged data file generated by SHM&CSR pipeline
    :param outdir: The output directory.
    :param genes: The genes to split out. Use '-' for all identified genes.
    :return:
    """
    match_dict = merged_txt_to_match_dict(merged_file)
    gene_outdirs = []
    gene_tarfiles = []
    os.makedirs(outdir, exist_ok=True)
    for gene in genes:
        gene_dir = f"{prefix}_{gene}" if gene else f"{prefix}"
        gene_dirpath = os.path.join(outdir, gene_dir)
        os.mkdir(gene_dirpath)
        gene_outdirs.append(gene_dirpath)
        new_filename = f"new_IMGT_{gene}.txz" if gene else "new_IMGT.txz"
        gene_tarfiles.append(
            tarfile.open(os.path.join(outdir, new_filename), mode="w:xz")
        )
    for name, table in imgt_to_tables(imgt_file):
        # Read each table one by one and per line select in which output
        # files it should go.
        gene_files = []
        for gene, gene_dir in zip(genes, gene_outdirs):
            fname = os.path.join(gene_dir, name)
            f = open(fname, "wt")
            gene_files.append((gene, f, fname))
        header = next(table)
        header_number_of_tabs = header.count('\t')
        column_names = header.strip("\n").split("\t")
        fr1_columns = [index for index, column in enumerate(column_names)
                       if column.startswith("FR1")]
        sequence_id_index = column_names.index("Sequence ID")
        for _, gene_file, _ in gene_files:
            gene_file.write(header)
        for line in table:
            # IMGT sometimes delivers half-empty rows.
            row_number_of_tabs = line.count("\t")
            missing_tabs = header_number_of_tabs - row_number_of_tabs
            if missing_tabs:
                line = line.strip("\n") + missing_tabs * "\t" + "\n"
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
    parser.add_argument(
        "genes",
        nargs="+",
        help="The genes to split out. Use '-' for all identified genes.")
    parser.add_argument("--prefix", help="Prefix for the archives and "
                                         "directories")
    return parser


def main():
    args = argument_parser().parse_args()
    genes = ["" if gene == "-" else gene for gene in args.genes]
    split_imgt(args.imgt_file, args.merged, args.outdir, genes,
               args.prefix)


if __name__ == "__main__":
    main()
