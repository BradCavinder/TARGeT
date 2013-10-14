
import sys
import os
import fnmatch
import os.path
import fastaIO


in_handle = open(sys.argv[1], "r")
info = in_handle.readlines()
filter_out = open(sys.argv[1] + "_filtered", "w")

header = 1
match = []
query_len = ''
query_start = ''
query_end = ''

for line in info:
    if header == 1:
        if line[0] != ">":
            #print "First match."
            print>>filter_out, line
            if " letters)" in line:
                line = line.rstrip()
                query_len = line.split(" letters)")[0].split("(")[1]
                #print "query length:", query_len
        else:
            header = 0
            match.append(line)
    elif query_start == '':
        #print "Need Query start."
        line = line.rstrip()
        match.append(line)
        if "Query: " in line:
            #print "Found Query start."
            query_start = line.split("Query: ")[1].split(" ")[0]
            #print "query start:", query_start
            line_split = line.split(" ")
            query_end = line_split[-1]
            #print "query end:", query_end
    else:
        #print "Still parsing hit, looking for new query end."
        if "Query: " in line:
            #print "Found next Query line, parsing for new end."
            line = line.rstrip()
            match.append(line)
            line_split = line.split(" ")
            query_end = line_split[-1]
            #print "query end:", query_end
        elif "Score =" in line or line[0] == ">":
            #print "Found next hit. Print previous hit if good and continue."
            line = line.rstrip()
            if query_start == "1" or query_end == query_len:
                for line2 in match:
                    print>>filter_out, line2
            match = []
            query_len = ''
            query_start = ''
            query_end = ''
            match.append(line)
        else:
            line = line.rstrip()
            match.append(line)
            

if query_start == "1" or query_end == query_len:
    for line2 in match:
        print>>filter_out, line2

in_handle.close()
filter_out.close()
            
            
    
            
            
