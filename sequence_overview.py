#!/usr/bin/env/python3

"""Create a HTML sequence overview"""

import argparse
import os
import typing
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List


class SequenceTableRow(typing.NamedTuple):
    sequence_id: str
    sequence: str
    best_match: str
    functionality: str


class SequenceStats:
    __slots__ = ("counts", "table_rows")

    def __init__(self):
        self.counts: Dict[str, int] = {
            "IGA1": 0,
            "IGA2": 0,
            "IGE": 0,
            "IGG1": 0,
            "IGG2": 0,
            "IGG3": 0,
            "IGG4": 0,
            "IGM": 0,
            "unmatched": 0}
        self.table_rows: List[SequenceTableRow] = []


def get_sequence_stats(before_unique: str,
                       sequence_columns: List[str]):
    sequence_statistics = defaultdict(SequenceStats)
    with open(before_unique, "rt") as table:
        header = next(table)
        header_columns = header.strip("\n").split("\t")
        for line in table:
            values = line.strip("\n").split("\t")
            row_dict = dict(zip(header_columns, values))
            sequence = "".join(row_dict[column] for column in sequence_columns)
            best_match = row_dict["best_match"]
            sequence_statistics[sequence].counts[best_match] += 1
            functionality = row_dict["Functionality"]
            sequence_statistics[sequence].table_rows.append(
                SequenceTableRow(row_dict["Sequence ID"], sequence,
                                 best_match, functionality))
    return sequence_statistics


def get_background_color(value: str):
    if value in ("TRUE", "T"):
        return "#eafaf1"
    elif value in ("FALSE", "F"):
        return "#f9ebea"
    try:
        flt = float(value)
    except ValueError:
        return "white"
    if flt > 0:
        return "#eaecee"
    return "white"


def td(val):
    return f"<td bgcolor='{get_background_color(val)}'>{val}</td>"


def tr(val: Iterable[str]):
    return f"<tr>{''.join(td(v) for v in val)}</tr>"


def make_link(link, val):
    return f"<a href='{link}'>{val}</a>"


def tbl(df: Iterable[Iterable[str]]):
    return f"<table border='1'>{''.join(tr(v) for v in df)}</table>"


def to_bool_str(cond):
    return "TRUE" if cond else "FALSE"


def sequence_overview(before_unique: str,
                      merged: str,
                      outdir: str,
                      gene_classes: List[str],
                      hotspot_analysis_sum: str,
                      empty_region_filter: str):
    os.makedirs(outdir, exist_ok=True)
    sequence_columns = [
        "FR1.IMGT.seq", "CDR1.IMGT.seq", "FR2.IMGT.seq", "CDR2.IMGT.seq",
        "FR3.IMGT.seq", "CDR3.IMGT.seq"]
    if empty_region_filter == "leader":
        sequence_columns = sequence_columns
    elif empty_region_filter == "FR1":
        sequence_columns = sequence_columns[1:]
    elif empty_region_filter == "CDR1":
        sequence_columns = sequence_columns[2:]
    elif empty_region_filter == "FR2":
        sequence_columns = sequence_columns[3:]
    else:
        raise ValueError(f"Unknown region filter: {empty_region_filter}")
    main_html_file = os.path.join(outdir, "index.html")
    by_id_file = os.path.join(outdir, "by_id.html")
    with open(main_html_file, "wt") as main_html, open(by_id_file, "wt") as by_id:
        main_html.write("<center><img src='data:image/png;base64,"
                        "iVBORw0KGgoAAAANSUhEUgAAAA8AAAAPCAYAAAA71pVKAAAAzElEQ"
                        "VQoka2TwQ2CQBBFpwTshw4ImW8ogJMlUIMmhNCDxgasAi50oSXA8X"
                        "lAjCG7aqKTzGX/vsnM31mzR0gk7tTudO5MEizpzvQ4ryUSe408J3X"
                        "n+grE0p1rnpOamVmWsZG4rS+dzzAMsN8Hi9yyjI1JNGtxu4VxBJgL"
                        "RLpoTKIPiW0LlwtUVRTubW2OBGUJu92cZRmdfbKQMAw8o+vi5v0fL"
                        "orZ7Y9waGYJjsf38DJz0O1PsEQffOcv4Sa6YYfDDJ5Obzbsp93+5Vf"
                        "dATueO1fdLdI0AAAAAElFTkSuQmCC'"
                        "> Please note that this tab is based on all "
                        "sequences before filter unique sequences and the "
                        "remove duplicates based on filters are applied. In "
                        "this table only sequences occuring more than once "
                        "are included. </center>")
        main_html.write("<table border='1' class='pure-table pure-table-striped'>")
        main_html.write(f"<caption>{'+'.join(sequence_columns)} sequences "
                        f"that show up more than once</caption>")
        main_html.write("<tr>")
        main_html.write("<th>Sequence</th><th>Functionality</th><th>IGA1</th>"
                        "<th>IGA2</th><th>IGG1</th><th>IGG2</th><th>IGG3</th>"
                        "<th>IGG4</th><th>IGM</th><th>IGE</th><th>UN</th>")
        main_html.write("<th>total IGA</th><th>total IGG</th><th>total IGM</th>"
                        "<th>total IGE</th><th>number of subclasses</th>"
                        "<th>present in both IGA and IGG</th>"
                        "<th>present in IGA, IGG and IGM</th>"
                        "<th>present in IGA, IGG and IGE</th>"
                        "<th>present in IGA, IGG, IGM and IGE</th>"
                        "<th>IGA1+IGA2</th>")
        main_html.write("</tr>")
        sequence_stats = get_sequence_stats(before_unique, sequence_columns)
        sorted_sequences = sorted(sequence_stats.keys())

        single_sequences = 0  # sequence only found once, skipped
        in_multiple = 0  # same sequence across multiple subclasses
        multiple_in_one = 0  # same sequence multiple times in one subclass
        unmatched = 0  # all the sequences are unmatched
        some_unmatched = 0  # one or more sequences in a clone are unmatched
        matched = 0  # should be the same als matched sequences

        for i, sequence in enumerate(sorted_sequences):
            sequence_stat: SequenceStats = sequence_stats[sequence]
            count_dict = sequence_stat.counts
            class_sum = sum(count_dict.values())
            if class_sum == 1:
                single_sequences += 1
                continue
            if count_dict["unmatched"] == class_sum:
                unmatched += 1
                continue
            in_classes = len([key for key in count_dict.keys() if key != "unmatched"])
            matched += in_classes
            if any(value == class_sum for value in count_dict.values()):
                multiple_in_one += 1
            elif count_dict["unmatched"] > 0:
                some_unmatched += 1
            else:
                in_multiple += 1
            functionality = ",".join(row.functionality
                                     for row in sequence_stat.table_rows)
            links: Dict[str, str] = {}
            for key, value in count_dict.items():
                name_key = "un" if key == "unmatched" else key
                html_file = f"{name_key}_{i}.html"
                links[key] = html_file
                if value > 0:
                    rows = [row for row in sequence_stat.table_rows
                            if row.best_match == key]
                    Path(outdir, html_file).write_text(tbl(rows))
                    for row in rows:
                        by_id.write(make_link(html_file, row.sequence_id) + "<br />")
            iga_count = count_dict["IGA1"] + count_dict["IGA2"]
            igg_count =  count_dict["IGG1"] + count_dict["IGG2"] + \
                count_dict["IGG3"] + count_dict["IGG4"],

            contained_classes = set(key for key, value in count_dict if value > 0)
            if iga_count:
                contained_classes.add("IGA")
            if igg_count:
                contained_classes.add("IGG")
            main_row = [
                sequence, functionality,
                make_link(links["IGA1"], count_dict["IGA1"]),
                make_link(links["IGA2"], count_dict["IGA2"]),
                make_link(links["IGG1"], count_dict["IGG1"]),
                make_link(links["IGG2"], count_dict["IGG2"]),
                make_link(links["IGG3"], count_dict["IGG3"]),
                make_link(links["IGG4"], count_dict["IGG4"]),
                make_link(links["IGM"], count_dict["IGM"]),
                make_link(links["IGE"], count_dict["IGE"]),
                make_link(links["unmatched"], count_dict["unmatched"]),
                make_link(links["IGA1"], count_dict["IGA1"]),
                iga_count,
                igg_count,
                count_dict["IGM"],
                count_dict["IGE"],
                in_classes,
                to_bool_str({"IGA", "IGG"}.issubset(contained_classes)),
                to_bool_str({"IGA", "IGG", "IGM"}.issubset(contained_classes)),
                to_bool_str({"IGA", "IGG", "IGE"}.issubset(contained_classes)),
                to_bool_str({"IGA", "IGG", "IGM", "IGE"}.issubset(contained_classes)),
                to_bool_str({"IGA1", "IGA2"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG2"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG3"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG4"}.issubset(contained_classes)),
                to_bool_str({"IGG2", "IGG3"}.issubset(contained_classes)),
                to_bool_str({"IGG2", "IGG4"}.issubset(contained_classes)),
                to_bool_str({"IGG3", "IGG4"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG2", "IGG3"}.issubset(contained_classes)),
                to_bool_str({"IGG2", "IGG3", "IGG4"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG2", "IGG4"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG3", "IGG4"}.issubset(contained_classes)),
                to_bool_str({"IGG1", "IGG2", "IGG3", "IGG4"}.issubset(contained_classes)),
            ]
            main_html.write(tr(main_row))
        main_html.write("</table>")


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--before-unique", help="File with the overview before unique filters")
    parser.add_argument("--merged", help="Merged.txt")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--gene-classes", help="Comma-separated list of gene classes")
    parser.add_argument("--hotspot-analysis-sum", help="hotspot-analysis-sum.txt")
    parser.add_argument("--nt-overview")
    parser.add_argument("--nt-sum")
    parser.add_argument("--empty-region-filter")
    return parser


def main():
    args = argument_parser().parse_args()
    sequence_overview(args.before_unique,
                      args.merged,
                      args.outdir,
                      args.gene_classes.split(","),
                      args.hotspot_analysis_sum,
                      args.empty_region_filter,
                      )
