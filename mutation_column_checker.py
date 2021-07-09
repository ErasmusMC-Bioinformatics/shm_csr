import re

mutationMatcher = re.compile("^([nactg])(\d+).([nactg]),?[ ]?([A-Z])?(\d+)?[>]?([A-Z;])?(.*)?")

with open("7_V-REGION-mutation-and-AA-change-table.txt", 'r') as file_handle:
    first = True
    fr3_index = -1
    for i, line in enumerate(file_handle):
        line_split = line.split("\t")
        if first:
            fr3_index = line_split.index("FR3-IMGT")
            first = False
            continue

        if len(line_split) < fr3_index:
            continue
        
        fr3_data = line_split[fr3_index]
        if len(fr3_data) > 5:
            try:
                test = [mutationMatcher.match(x).groups() for x in fr3_data.split("|") if x]
            except:
                print((line_split[1]))
                print(("Something went wrong at line {line} with:".format(line=line_split[0])))
                #print([x for x in fr3_data.split("|") if not mutationMatcher.match(x)])
        if i % 100000 == 0:
            print(i)
