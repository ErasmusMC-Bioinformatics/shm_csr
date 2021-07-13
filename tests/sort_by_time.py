#!/usr/bin/env python3

"""Small script to profile bash scripts that have been run with the following
code inside:

    exec 5> debug_output.txt
    BASH_XTRACEFD="5"
    PS4='\t $LINENO: '
    set -x


"""
import calendar
import time
import sys


def file_to_timestamped_lines(input_file):
    with open(input_file, "rt") as file_h:
        for line in file_h:
            try:
                time_object = time.strptime(line[0:8], "%H:%M:%S")
            except ValueError:
                # Account for misformed lines
                continue
            time_since_epoch = calendar.timegm(time_object)
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
        if time_since == 0:
            continue
        print(time_since, line.strip())
