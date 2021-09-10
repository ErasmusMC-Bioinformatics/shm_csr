#!/usr/bin/env python3

"""Small script to profile bash scripts that have been run with the following
code inside:

    exec 5> debug_output.txt
    BASH_XTRACEFD="5"
    PS4='$(date +%s.%N) $LINENO: '
    set -x


"""
import calendar
import time
import sys

import re

SECONDS_FINDER = re.compile(r"^(\d+.\d+).*")


def file_to_timestamped_lines(input_file):
    with open(input_file, "rt") as file_h:
        for line in file_h:
            time_since_epoch = float(SECONDS_FINDER.search(line).group(1))
            yield time_since_epoch, line


def time_delta_lines(input_file):
    timestamped_lines = file_to_timestamped_lines(input_file)
    current_time, current_line = next(timestamped_lines)
    for next_time, next_line in timestamped_lines:
        time_since = next_time - current_time
        yield time_since, current_line
        current_time = next_time
        current_line = next_line


if __name__ == "__main__":
    input_file = sys.argv[1]
    # Sort by time ascending order.
    sorted_time = sorted(time_delta_lines(input_file), key=lambda tup: tup[0])
    for time_since, line in sorted_time:
        if time_since > 60*60*24*365:
            # big times are probably nonsensical parsing errors.
            continue
        print(time_since, line.strip())
