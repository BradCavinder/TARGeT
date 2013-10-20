
import sys
import os
import fnmatch
import os.path
import fastaIO


in_handle = open(sys.argv[1], "r")
info = in_handle.readlines()
window = int(sys.argv[2])

match = []
query_len = ''
query_start = ''
query_end = ''
title_dict = {}
wanted = []


for line in info:
    line = line.rstrip()
    copy = line.split(" ")[0]
    line_split = line.split("  ")
    #print "line split:"
    #print line_split
    query_len = line_split[0].split(" ")[2]
    query_start = int(line_split[1].split(" ")[0])
    query_end = line_split[-1].split(" ")[1]
    #print "query_len:", query_len, "query_start:", query_start, "query_end:", query_end
    if query_start <= (1 + window) and query_end >= (query_len - window):
        match.append(line + "  ")
        title_dict[copy] = 1

in_handle.close()

filter_out = open(sys.argv[1], "w")
for item in match:
    print>>filter_out, item
filter_out.close()

flank_path = os.path.splitext(sys.argv[1])[0] + ".flank"
in_handle = open(flank_path, "r")
for title, seq in fastaIO.FastaGeneralIterator(in_handle):
    title_check = title.split(" ")[0]
    if title_check in title_dict:
        wanted.append([title, seq])
in_handle.close()

out_handle = open(flank_path, "w")
for item in wanted:
    print>>out_handle, ">" + item[0] + "\n" + item[1]
out_handle.close()
