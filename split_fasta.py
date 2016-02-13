import sys
import os
import re
import fnmatch
import os.path
import subprocess as subp
import fastaIO

in_name = sys.argv[1]
in_file = os.path.split(in_name)[1]
in_trim = os.path.splitext(in_file)[0]
out_dir = os.path.split(in_name)[0]
in_handle = open(in_name, "r")
c = 1
for title, seq in fastaIO.FastaTitleStandardization(in_handle):
    out_handle = open(os.path.join(out_dir, in_trim + "_split" + str(c) + ".fa"), 'w')
    print>>out_handle, ">" + title, "\n", seq
    c += 1
    out_handle.close()
in_handle.close()
