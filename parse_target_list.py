
import sys
import os
import fnmatch
import os.path
import fastaIO


in_handle = open(sys.argv[1], "r")
info = in_handle.readlines()
match = []
query_len = ''
query_start = ''
query_end = ''
c = 1
for line in info:
    line = line.rstrip()
    line_split = line.split("  ")
    #print "line split:"
    #print line_split
    query_len = line_split[0].split(" ")[2]
    query_start = line_split[1].split(" ")[0]
    query_end = line_split[-1].split(" ")[1]
    #print "query_len:", query_len, "query_start:", query_start, "query_end:", query_end
    if query_start == "1" and query_end == query_len:
        line = str(c) + "_" + line.split("_", 1)[1]
        match.append(line + "  ")
        c += 1
in_handle.close()

filter_out = open(sys.argv[1], "w")
for item in match:
    print>>filter_out, item
filter_out.close()
