#!/usr/bin/env python

"""
This is the Python wrapper script for running the command line version of TARGeT. Please see the README for dependencies, installation instructions, and change log. 

From the directory containing TARGeT, run 
python target.py -h
for help.

Brad Cavinder
"""

import os
import os.path
import sys
import datetime
import glob
import subprocess as subp
import argparse
import fnmatch
import re
import fastaIO

path = sys.path[0]
path = str(path) + "/"

#-----------Define functions-----------------------------------------------------

def frange(x, y, jump):
    frange_list = []
    while x < y:
        frange_list.append(x)
        x += jump
        continue
    return frange_list
def frange_convert(frange_list):
    convert_list = []
    for i in frange_list:
        i = str(i)
        i = i[:4]
        convert_list.append(float(i))
    return convert_list

#setup available floating point options
ident_list1 = frange(0.01, 1, 0.01)
ident_list2 = frange_convert(ident_list1)
ident_list3 = frange(0, 10, 0.01)
ident_list4 = frange_convert(ident_list3)
    

def BLASTN(query, blast_file_out):
    subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.P)])

def TBLASTN(query, blast_file_out):
    subp.call(["blastall", "-p", "tblastn", "-F", "F", "-m", "0", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.P)])

def Blast_draw(blast_file_out):
    subp.call(["perl", path + "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])

def img_convert(in_file, out_file):
    subp.call(["convert", in_file, out_file])

def PHI(blast_in, PHI_out):
    subp.call(["perl", path + "PHI_2.4.pl", "-i", blast_in, "-q", query, "-D", args.genome, "-o", PHI_out, "-e", str(args.p_e), "-M", str(args.p_M), "-P", Type, "-d", str(args.p_d), "-g", str(args.p_g), "-n", str(args.p_n), "-c", str(args.p_c), "-G", str(args.p_G), "-t", str(args.p_t), "-f", str(args.p_f), "-p", args.p_p, "-R", realign])

def PHI_draw(PHI_out, Type):
    subp.call(["perl", path + "PHI_drawer2.pl", "-i", str(PHI_out) + ".list", "-o", str(PHI_out) + ".tcf_drawer", "-m", args.p_t, "-P", Type, "-n", str(args.p_n)])

def MAFFT_NT(in_path, out_path):
    subp.call(["mafft", "--ep", "0.75", "--op", "3.0", "--thread",  str(args.P), "--localpair", "--maxiterate",  "4", "--out", out_path, in_path])

def MAFFT_P(in_path, out_path):
    subp.call(["mafft", "--thread",  str(args.P), "--maxiterate",  "100", "--out", out_path, in_path])

def runTarget(query, blast_out, blast_file_out):
    #make output directory
    subp.call(["mkdir", blast_out])

    #use blastn if DNA
    if args.Type == 'nucl':
        print "Using BLASTN\n"
        BLASTN(query, blast_file_out)
        
    #use tblastn if protein
    elif args.Type == 'prot':
        print "Using TBLASTN\n"
        TBLASTN(query, blast_file_out)

    #make svg drawing(s)
    print "Making svg image of blast results\n"
    Blast_draw(blast_file_out)

    #convert svg image to jpg
    print "Converting svg to jpg\n"
    for svg_file in glob.glob(str(blast_file_out) + "*.svg"): 
        jpg_file = os.path.splitext(svg_file)
        jpg_file = jpg_file[0] + ".jpg"
        img_convert(svg_file, jpg_file)
        
    if args.S == 'Blast':
        return
        
    blast_in = str(blast_file_out) + ".blast"
    PHI_out = str(blast_file_out) + ".tcf"
    print "Running PHI"
    PHI(blast_in, PHI_out)
    print "PHI finished!\n"

    #make svg image of PHI homologs
    print "Making svg image of homologs\n"
    PHI_draw(PHI_out, Type)

    #convert svg to pdf
    print "Coverting svg image to pdf\n"
    img_convert(str(PHI_out) + ".tcf_drawer.svg", str(PHI_out) + ".tcf_drawer.pdf")
    
    if args.S == 'PHI':
        return

    #get query length
    query_in = open(query, "r")
    query_len = 0
    for title, seq in fastaIO.FastaGeneralIterator(query_in):
        query_len = len(seq)
    query_in.close()
    
    #check that two or more copies were found
    copies = 0
    dna_copies_in = open(PHI_out + ".dna", "r")
    for title, seq in fastaIO.FastaGeneralIterator(dna_copies_in):
        copies+=1
    dna_copies_in.close()
    
    if copies >= 2: 
        filter_list = []
        in_list = []
        if args.a == 'hits' or args.a == 'both':
            print "hits will be aligned"
            if args.f > 0:
                print "Filtering flagged for hits\n"
                if args.Type == 'nucl':
                    filter_list.append([PHI_out + ".dna", PHI_out + ".dna_filter-" + str(args.f)])
                else:
                    in_list.append(str(PHI_out) + ".aa")
            else:
                if args.Type == 'nucl':
                    in_list.append(str(PHI_out) + ".dna")
                else:
                    in_list.append(str(PHI_out) + ".aa")

        if args.a == 'flanks' or args.a == 'both':
            print "Flanks will be aligned"
            if args.f > 0:
                print "Filtering flagged for flanks\n"
                if args.Type == 'nucl':
                    filter_list.append([PHI_out + ".flank", PHI_out + ".flank_filter-" + str(args.f)])
                else:
                    in_list.append(str(PHI_out) + ".flank")
            else:
                in_list.append(str(PHI_out) + ".flank")
        print "Entries in filter list:  ", len(filter_list), "\n"
        if len(filter_list) != 0:
            #print "in_list =\n", in_list
            for in_path, out_path in filter_list:
                in_file = open(in_path, "r")
                out_file = open(out_path, "w")
                for title, seq in fastaIO.FastaGeneralIterator(in_file):
                    copy_len = len(seq) - (int(args.p_f) * 2)
                    if copy_len <= (query_len * args.f):
                        print>>out_file, ">" + title + "\n" + seq
                in_list.append(out_path)
                in_file.close()
                out_file.close()

        #Run Mafft
        copies = 0
        in_count = len(in_list)
        processed = 0
        for in_path in in_list:
            in_file = open(in_path, "r")
            for title, seq in fastaIO.FastaGeneralIterator(in_file):
                copies += 1
            print str(copies) + " copies in " + in_path, "\n"
            msa_out = in_path + ".msa"
            print "Running Mafft"
            if args.Type == 'nucl':
                MAFFT_NT(in_path, msa_out)
            else:
                MAFFT_P(in_path, msa_out)
            processed += 1
            
            #print "in_count - processed = ", in_count - processed, "\n"
            if args.S == 'MSA':
                if (in_count - processed) == 0:
                    return
                else:
                    continue

            #Run FastTreeMP
            print "Running FastTreeMP\n"
            tree_out = msa_out + ".nw"
            print "Output tree path: ", tree_out, "\n\n"
            
            #Can only limit FastTree processor use through OMP_NUM_THREADS. Otherwise, it will use all processors available when compiled.
            current_env = os.environ.copy()
            #print "current_env before change:  ", current_env
            current_env['OMP_NUM_THREADS'] = str(args.P)
            #print "OMP_NUM_THREADS after change:  ", current_env['OMP_NUM_THREADS'], "\n\n"

            
            if args.Type == 'nucl':
                proc = subp.Popen(["FastTreeMP", "-nt", "-gamma", "-out", tree_out, msa_out], env = current_env)
                proc.wait()
                print "\nFastTreeMP finished.\n"
            else:
                proc = subp.Popen(["FastTreeMP", "-gamma", "-out", tree_out, msa_out], env = current_env)
                proc.wait()
                print "\nFastTreeMP finished.\n"

            print "Converting output tree file to eps image\n"
            out = open(tree_out + ".eps", "w") #open output file for redirected stdout
            #print "Eps image out path: ", out
            
            if copies > 45:
                height = copies * 12
                width = round(height/3)
                print "Image height: ", height, "\twidth: ", width, "\n"
                subp.call(["treebest",  "export", "-y", str(height), "-x", str(width), "-b", "4.5", "-f", "13", tree_out], stdout=out)
            else:
                subp.call(["treebest",  "export", tree_out], stdout=out)
            out.close() #close output file

            print "Coverting eps image to pdf\n"
            subp.call(["convert", tree_out + ".eps", tree_out + ".pdf"])
    else:
        print "Less than two copies found. Multiple alignment and tree building will not be performed.\n"
            
    
#-----------Make command line argument parser------------------------------------

parser = argparse.ArgumentParser(prog='target.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="This is the command line implementation of TARGeT: Tree Analysis of Related Genes and Transposons (http://target.iplantcollaborative.org/). Please read the README file for program dependencies.")

input_group = parser.add_mutually_exclusive_group(required=True)

input_group.add_argument("-q", metavar="Input file", help="Path to single input query sequence file (may contain multiple sequences of the same type [DNA or protein]).")

input_group.add_argument("-d", metavar="Input directory", help="Path to input directory with multiple query sequence files (each may contain multiple sequences). All sequences must be of the same type (DNA or protein). The directory can not contain non-query files.")

parser.add_argument("-t", dest="Type", metavar="Search type", choices=("nucl", "prot"), default="prot", help="Type of input query sequence(s): DNA = 'nucl'  protein = 'prot'")

parser.add_argument("genome", metavar="Genome_sequence_file", help="Path of the genome sequence file to be searched (or other sequence to be searched).")

parser.add_argument("Run_Name", help="Name for overall run")

parser.add_argument("-o", metavar="Output path", default="Current working directory", help="Path to an existing directory for output. A subdirectory for the run will be created with further subdirectories containing the output files for each search sequence.")

parser.add_argument("-i", metavar="Input query type", choices=("s", "mi", "g"), default="s", help="Number of input queries per file(s): s = single query sequence, mi = multiple individual query sequences, g = FASTA multiple sequence alignment or list of sequences; All input files/sequences must be the same type of sequence as set by -t.")

parser.add_argument("-a", metavar="Alignments to perform", choices=("hits", "flanks", "both"), default="hits", help="Input sequences for multiple sequence alignent: hits = each copy only contains sequence that matches the query (recomended for most protein queries), flanks = each copy contains the sequences that matches the query plus flanking sequence on each side of the match for which the length is set by -p_f ")

parser.add_argument("-f", metavar="Filter length (query length * X)", type=float, default=0, choices=ident_list4, help="Multiple of query length used as maximum length of copies to be included in the multiple sequence alignment(s). Valid values are: 0-10 by 0.01 increments. Default of 0 means no filtering. Values between 0 and 1 return copies smaller than the input query sequence.")

parser.add_argument("-P", metavar="Processors", type=int, default=1, help="The number of processors to use for Blast, Mafft, and FastTree steps. All other steps use 1 processor. The programs are not multi-node ready, so the number of processors is limited to that available to one computer/node.")

parser.add_argument("-S", metavar="Stopping point of program", type=str, choices=("Blast", "PHI", "MSA", "Tree"), default="Tree", help="The stage, after completion, to stop the program. By default, all stages (Blast, PHI, MSA, Tree) are run. For example if you want to stop the program after Blast and PHI, exiting before the MSA stage, enter PHI.")

parser.add_argument("-v", action='version', version='TARGeT-2.00', help="Version information")

#BLAST arguments
parser_blast = parser.add_argument_group("BLAST")

parser_blast.add_argument("-b_e", metavar="BLAST E-value", default="0.001", type=float, help="Maximum Expect value to include a match in the BLAST results")

parser_blast.add_argument("-b_a", metavar="# alignments", default="500", type=int, help="The maximum number of hit alignments to include per query sequence.")

parser_blast.add_argument("-b_d", metavar="# descriptions", default="100", help="The maximum number of descriptions of hits included in BLAST output per query sequence.")


#PHI arguments
parser_phi = parser.add_argument_group("PHI")

parser_phi.add_argument("-p_e", metavar="PHI E-value", default="0.01", type=float, help="Maximum Expect value for PHI.")

parser_phi.add_argument("-p_M", metavar="Minimum query match", type=float, default=0.7, choices=ident_list2, help="Minimum decimal fraction of match length.")

parser_phi.add_argument("-p_d", metavar="Max intron length", default="8000", type=int, help="Maximum intron length.")

parser_phi.add_argument("-p_g", metavar="Min intron length", type=int,default=10, help="Minimum intron length.")

parser_phi.add_argument("-p_n", metavar="Max copies", type=int, default=100, help="Maximum number of homolgs in output.")

parser_phi.add_argument("-p_R", action='store_true', help="Perform realignments to find small exons.")

parser_phi.add_argument("-p_c", metavar="Composition stats", default="0", choices=('0', '1', '2', '3'), help="Use composition-based statistics during realignment. 0 = off, 1 = use 2001-based stats, 2 = use 2005-based stats cibditioned on sequence properties, 3 = use 2005-based stats unconditionally. See BLAST+ help, man page, or web resources for further explanation.")

parser_phi.add_argument("-p_G", metavar="Add GenBank accession #", default="0", choices=('0', '1'), help="Add GenBank accession numbers. Use only if genome being searched is in GenBank format. 0 = no, 1 = yes.")

parser_phi.add_argument("-p_t", metavar="Substitution matrix", default="62", choices=('50', '62', '80'), help="Substitution matrix used to score alignments. 50 = BLOSUM50, 62 = BLOSUM62, 80 = BLOSUM80. DNA queries will automatically use an identity matrix.")

parser_phi.add_argument("-p_f", metavar="Flanking sequence", type=int, default=100, help="Length of flanking sequences to be included.")

parser_phi.add_argument("-p_p", action='store_true', help="Filter out psuedogenes (matches with internal stop codons).")


#-----------Parse command line input---------------------------------------------

args = parser.parse_args()

#change arguments to usable forms or make new assignment
realign = '0'
#print "realign before = ", realign
if args.Type == 'prot':
    Type = '1'
elif args.Type == 'nucl':
    Type = '0'
args.p_t = path + 'BLOSUM' + str(args.p_t) + '.txt'
if args.p_R == 1:
    realign = '1'
if args.p_p == 'True':
    args.p_p = '1'
else:
    args.p_p = '0'



#-----------Create main output directory-----------------------------------------

if args.o != "Current working directory":
    if args.o[-1] != "/":
        args.o = args.o + "/"
else:
    args.o = ""
#make the main output directory using current date and time
now = datetime.datetime.now()


if args.q and args.i == 's' or args.i == 'g':
    query_no_ext = os.path.splitext(args.q)[0]
    basedir, query_name = os.path.split(query_no_ext)
    out_dir = args.o + args.Run_Name + "_" + query_name + "_" + now.strftime("%Y_%m_%d_%H%M%S")
else:
    out_dir = args.o + args.Run_Name + "_" + now.strftime("%Y_%m_%d_%H%M%S")
    try:
        subp.call(["mkdir", out_dir])
    except:
        pass


#-----------Make BLAST database -------------------------------------------------
#print "Genome:  ", args.genome, "\n"
name = args.genome.split('/')
name = name[-1]
print "\nRunning formatdb\n"
subp.call(["formatdb", "-i", args.genome,"-p", "F", "-o", "T", "-n", args.genome])


#-----------Run custom indexer---------------------------------------------------
print "Running custom indexer\n"
subp.call(["perl", path + "reads_indexer.pl", "-i", args.genome])



#-----------Single input file----------------------------------------------------

#for single indiviual query or grouped query
if args.q and args.i == 's' or args.i == 'g':
    print "Single input file, single or group input\n"
    #set output filename
    blast_file_out = os.path.join(out_dir, query_name)
    runTarget(args.q, out_dir, blast_file_out)
    print "TARGeT has finished!"


#for multiple individual queries-------------------------------------------------
elif args.q and args.i == 'mi':
    print "Single input file, multiple individual inputs.\n"
    subp.call(["python", "split_fasta.py", args.q])
    
    #count the number of new input files
    query_no_ext = os.path.splitext(args.q)[0]
    basedir, query_name = os.path.split(query_no_ext)
    new_files = glob.glob(query_no_ext + "_*.fa")
    print "new files:\n", new_files
    print str(query_no_ext).strip() + "_*.fa\n"
    count = 0
    for f in new_files:
        count += 1        
    print count, " files to be processed"

    #counter to keep track of the number of files that have been processed
    p = 0

    #Run pipeline on each file with it's own output directory in the main output directory
    for fasta2 in new_files:
        filename = os.path.splitext(fasta2)[0]
        file_name = os.path.split(filename)[1]
        #set output directory
        blast_out = os.path.normpath(os.path.join(out_dir, file_name))
        #set output filename
        blast_file_out = os.path.join(blast_out, file_name + ".blast")
        runTarget(fasta2, blast_out, blast_file_out)
        p += 1
        print "TARGeT has processed ", p, " of ", count, " subfiles"

    print "TARGeT has finished!"

#-----------Directory input------------------------------------------------------

#for single individual queries or grouped queries
elif args.d and args.i == 's' or args.i == 'g':
    print "Directory input, each file has a single or group input.\n"
    
    files = os.listdir(args.d) #get all files in the directory

    #count the number of input files
    count = 0
    for f in files:
        count += 1        
    print count, " files to be processed"

    #counter to keep track of the number of files that have been processed
    p = 0
 
    #Run pipeline on each file with it's own output directory in the main output directory
    for f in files:
        query = os.path.normpath(os.path.join(args.d, f))
        query_name = os.path.split(os.path.splitext(query)[0])[1]
        blast_out = os.path.normpath(os.path.join(out_dir, query_name)) #output directory
        blast_file_out = os.path.join(blast_out, query_name) #output files basename
        runTarget(query, blast_out, blast_file_out)
        p +=1
        print "TARGeT has processed ", p, " of ", count, " files"
    print "TARGeT has finished!"

#for multiple individual queries-------------------------------------------------
elif args.d and args.i == 'mi':
    print "Directory input, each file has multiple individual queries.\n"
    
    files = os.listdir(args.d) #get all files in the directory

    #count the number of input files
    count = 0
    for f in files:
        count += 1
        print f        
    print count, " files to be processed"

    #counter to keep track of the number of files that have been processed
    p = 0

    for f in files:
        #setup name for output subdirectory for pre-split file
        base_dir = os.path.normpath(os.path.join(out_dir, os.path.splitext(f)[0]))
        
        #make the subdirectory
        subp.call(["mkdir", base_dir])

        #setup full path to fasta file for splitting
        in_full = os.path.normpath(os.path.join(args.d, f))
        in_no_ext = os.path.splitext(in_full)[0]
        
        #split the fasta file
        subp.call(["python", "split_fasta.py", in_full])

        #count the number of new input files
        new_files = glob.glob(in_no_ext + "_*.fa")
        count2 = len(new_files)
        print count2, " subfiles to be processed"

        #setup conter for subfiles processed
        p2 = 0

        #Run pipeline on each split file, output directory in subdirectory for pre-split file that's in main output directory
        for fasta2 in new_files:
            #set query to split file path
            blast_out = os.path.normpath(os.path.join(base_dir, os.path.split(os.path.splitext(fasta2)[0])[1]))
            blast_file_out = os.path.normpath(os.path.join(blast_out, os.path.split(fasta2)[1]))
            runTarget(fasta2, blast_out, blast_file_out)
            p2 += 1
            print "TARGeT has processed ", p2, " of ", count2, " subfiles from file ", p + 1, " of ", count, " files"
            p += 1
        print "TARGeT has fully processed ", p, " of ", count, " files"
    print "TARGeT has finished!"
