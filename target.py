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
    subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.b_p)])

def TBLASTN(query, blast_file_out):
    subp.call(["blastall", "-p", "tblastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.b_p)])

def Blast_draw(blast_file_out):
    subp.call(["perl", path + "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])

def img_convert(in_file, out_file):
    subp.call(["convert", in_file, out_file])

def PHI(blast_in, PHI_out):
    subp.call(["perl", path + "PHI_2.4.pl", "-i", blast_in, "-q", query, "-D", args.genome, "-o", PHI_out, "-e", str(args.p_e), "-M", str(args.p_M), "-P", Type, "-d", str(args.p_d), "-g", str(args.p_g), "-n", str(args.p_n), "-c", str(args.p_c), "-G", str(args.p_G), "-t", str(args.p_t), "-f", str(args.p_f), "-p", args.p_p, "-R", realign])

def PHI_draw(PHI_out, Type):
    subp.call(["perl", path + "PHI_drawer2.pl", "-i", str(PHI_out) + ".list", "-o", str(PHI_out) + ".tcf_drawer", "-m", args.p_t, "-P", Type, "-n", str(num_hom)])

def MUSCLE(in_path, out_path):
    subp.call(["muscle", "-in", in_path, "-out", out_path, "-maxiters", str(args.m_m), str(m_args)])

def runTarget(query, blast_out, blast_file_out):
    #make output directory
    subp.call(["mkdir", blast_out])

    #use blastn if DNA
    if args.Type == 'nucl':
        print "Using BLASTN"
        BLASTN(query, blast_file_out)
        
    #use tblastn if protein
    elif args.Type == 'prot':
        print "Using TBLASTN"
        TBLASTN(query, blast_file_out)

    #make svg drawing(s)
    print "Making svg image of blast results"
    Blast_draw(blast_file_out)

    #convert svg image to jpg
    print "Converting svg to jpg"
    for svg_file in glob.glob(str(blast_file_out) + "*.svg"): 
        jpg_file = os.path.splitext(svg_file)
        jpg_file = jpg_file[0] + ".jpg"
        img_convert(svg_file, jpg_file)
        
    blast_in = str(blast_file_out) + ".blast"
    PHI_out = str(blast_file_out) + ".tcf"
 
    print "Running PHI"
    PHI(blast_in, PHI_out)
    print "PHI finished!"

    #make svg image of PHI homologs
    print "Making svg image of homologs"
    num_hom = args.p_n
    if num_hom > 500:
        num_hom = 500
    PHI_draw(PHI_out, Type)

    #convert svg to jpg
    print "Coverting svg image to jpg"
    img_convert(str(PHI_out) + ".tcf_drawer.svg", str(PHI_out) + ".tcf_drawer.jpg")

    #Count the number of homologs to determine whether to run MUSCLE and TreeBest and possibly adjust tree image height
    query_in = open(query, "r")
    query_len = 0
    filter_list = []
    in_list = []
    for title, seq in fastaIO.FastaGeneralIterator(query_in):
        query_len = len(seq)
    query_in.close()
    
    seq1 = 0
    seq2 = 0
    seq3 = 0
    seq4 = 0
    if str(args.a) == 'hits' or str(args.a) == 'both':
        if args.f != 0:
            filter_list.append([str(PHI_out) + ".dna", str(PHI_out) + ".dna_filter-" + str(args.f)])
        else:
            in_list.append([str(PHI_out) + ".dna"])

    if str(args.a) == 'flanks' or str(args.a) == 'both':
        if args.f != 0:
            filter_list.append([str(PHI_out) + ".flank", str(PHI_out) + ".flank_filter-" + str(args.f)])
        else:
            in_list.append([str(PHI_out) + ".flank"])
    
    if int(args.f) != 0:        
        for in_path, out_path in filter_list:
            in_file = open(in_path, "r")
            out_file = open(out_path, "w")
            for title, seq in fastaIO.FastaGeneralIterator(in_file):
                copy_len = len(seq) - (int(args.p_f) * 2)
                if copy_len <= query_len * args.f:
                    print>>out_file, ">" + title + "\n" + seq
                    if '.dna' in in_path:
                        seq1 += 1
                    elif '.flank' in in_path:
                        seq2 += 1
            in_list.append(out_path)
            in_file.close()
            out_file.close()
    else:
        for in_path in in_list:
            if '.dna' in in_path:
                seq3 += 1
            elif '.flank' in in_path:
                seq4 += 1

    if seq1 > 1 or seq2 > 1 or seq3 > 1 or seq4 > 1:

        #setup Muscle arguments
        m_args = ''
        if args.m_d == 'diags':
            m_args = "-diags"
        if args.m_g == 'input':
            m_args = m_args + "-input"
        if args.m_c == 'True':
            m_args = m_args + "-clwstrict"
    

        #Run Muscle
        for in_path in in_list:
            print in_path + "\n"
            in_file2 = open(in_path, "r")
            copy_count = 0
            for title, seq in fastaIO.FastaGeneralIterator(in_file2):
                copy_count += 1
            print str(copy_count) + " copies in " + in_path
            muscle_out = in_path + ".msa"
            print "Running Muscle on filtered hits"
            MUSCLE(in_path, muscle_out)

            #Run TreeBest
            tree_out = muscle_out + ".nw"
            out = open(tree_out, "w") #open output file for redirected stdout
            subp.call(["treebest", "nj", muscle_out], stdout=out)
            out.close() #close output file
            print "TreeBest finished, converting output file to eps image"

            out = open(tree_out + ".eps", "w") #open output file for redirected stdout
            if seq > 45:
                height = seq * 10
                subp.call(["treebest",  "export", "-y", str(height), tree_out], stdout=out)
            else:
                subp.call(["treebest",  "export", tree_out], stdout=out)
            out.close() #close output file

            print "Coverting eps image to jpg"
            subp.call(["convert", tree_out + ".eps", tree_out + ".jpg"])
    else:
        print "One or less copies found"
    
#-----------Make command line argument parser------------------------------------

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="This is the command line implementation of TARGeT: Tree Analysis of Related Genes and Transposons (http://target.iplantcollaborative.org/). Please read the README file for program dependencies.")

input_group = parser.add_mutually_exclusive_group(required=True)

input_group.add_argument("-q", metavar="Input file", help="Path to single input query sequence file (may contain multiple sequences of the same type [DNA or protein]).")

input_group.add_argument("-d", metavar="Input directory", help="Path to input directory with multiple query sequence files (each may contain multiple sequences). All sequences must be of the same type (DNA or protein). The directory can not contain non-query files.")

parser.add_argument("-t", dest="Type", metavar="Search type", choices=("nucl", "prot"), default="prot", help="Type of input query sequence(s): DNA = 'nucl'  protein = 'prot'")

parser.add_argument("genome", metavar="Genome sequence file", help="Path of the genome sequence file to be searched (or other sequence to be searched).")

parser.add_argument("Run_Name", help="Name for overall run")

parser.add_argument("-o", metavar="Output path", default="Current working directory", help="Path to an existing directory for output. A subdirectory for the run will be created with further subdirectories containing the output files for each search sequence.")

parser.add_argument("-i", metavar="Input query type", choices=("s", "mi", "g"), default="s", help="Number of input queries per file(s): s = single query sequence, mi = multiple individual query sequences, g = FASTA multiple sequence alignment or list of sequences; All input files/sequences must be the same type of sequence as set by -t.")

parser.add_argument("-a", metavar="Alignments to perform", choices=("hits", "flanks", "both"), default="hits", help="Input sequences for multiple sequence alignent: hits = each copy only contains sequence that matches the query (recomended for most protein queries), flanks = each copy contains the sequences that matches the query plus flanking sequence on each side of the match for which the length is set by -p_f ")

parser.add_argument("-f", metavar="Filter length (query length * X)", type=float, default=0, choices=ident_list4, help="Multiple of query length used as maximum length of copies to be included in the multiple sequence alignment(s). Valid values are: 0-10 by 0.01 increments. Default of 0 means no filtering. Values between 0 and 1 return copies smaller than the input query sequence.")



#BLAST arguments
parser_blast = parser.add_argument_group("BLAST")

parser_blast.add_argument("-b_e", metavar="BLAST E-value", default="0.001", type=float, help="Maximum Expect value to include a match in the BLAST results")

parser_blast.add_argument("-b_a", metavar="# alignments", default="500", type=int, help="The maximum number of hit alignments to include per query sequence.")

parser_blast.add_argument("-b_d", metavar="# descriptions", default="100", help="The maximum number of descriptions of hits included in BLAST output per query sequence.")

parser_blast.add_argument("-b_o", metavar="Other BLAST options", help="Other options: you may use any other valid Blast options for the search type being used except those that change the output format from default (including the use of GenBank accession numbers). Place them all here as you would for running standalone BLAST")

parser_blast.add_argument("-b_p", metavar="# processors", default="1", type=int, help="The number of processors to use for BLAST stage. All other stages currently use only 1 processor.")



#PHI arguments
parser_phi = parser.add_argument_group("PHI")

parser_phi.add_argument("-p_e", metavar="PHI E-value", default="0.01", type=float, help="Maximum Expect value for PHI.")

parser_phi.add_argument("-p_M", metavar="Minimum query match", type=float, default=0.7, choices=ident_list2, help="Minimum decimal fraction of match length.") #choices=frange(0,1,0.01)

parser_phi.add_argument("-p_d", metavar="Max intron length", default="8000", type=int, help="Maximum intron length.")

parser_phi.add_argument("-p_g", metavar="Min intron length", type=int,default=10, help="Minimum intron length.")

parser_phi.add_argument("-p_n", metavar="Max copies", type=int, default=100, help="Maximum number of homolgs in output.")

parser_phi.add_argument("-p_R", action='store_true', default="False", help="Perform realignments to find small exons.")

parser_phi.add_argument("-p_c", metavar="Composition stats", default="0", choices=('0', '1', '2', '3'), help="Use composition-based statistics during realignment. 0 = off, 1 = use 2001-based stats, 2 = use 2005-based stats cibditioned on sequence properties, 3 = use 2005-based stats unconditionally. See BLAST+ help, man page, or web resources for further explanation.")

parser_phi.add_argument("-p_G", metavar="Add GenBank accession #", default="0", choices=('0', '1'), help="Add GenBank accession numbers. Use only if genome being searched is in GenBank format. 0 = no, 1 = yes.")

parser_phi.add_argument("-p_t", metavar="Substitution matrix", default="62", choices=('50', '62', '80'), help="Substitution matrix used to score alignments. 50 = BLOSUM50, 62 = BLOSUM62, 80 = BLOSUM80. DNA queries will automatically use an identity matrix.")

parser_phi.add_argument("-p_f", metavar="Flanking sequence", type=int, default=100, help="Length of flanking sequences to be included.")

parser_phi.add_argument("-p_p", action='store_true', default="False", help="Filter out psuedogenes (matches with internal stop codons).")



#Muscle arguments
parser_muscle = parser.add_argument_group("Muscle")

parser_muscle.add_argument("-m_m", metavar="Max alignment iterations", type=int, default=8, help="Maximum number of alignment iterations.")

parser_muscle.add_argument("-m_d", metavar="Find diagonals", choices=('diags', ''), default='diags', help="Find diagonals (faster for similar sequences)")

parser_muscle.add_argument("-m_g", metavar="Grouping", choices=('group', 'input'), default="group", help="Output sequences in input or group (by similarity) order.")

parser_muscle.add_argument("-m_c", action='store_true', default="False", help="Write output in Clustal W format (default format is FASTA) with 'CLUSTAL W (1.81)' header.")


#-----------Parse command line input---------------------------------------------

args = parser.parse_args()

#change arguments to usable forms or make new assignment
if args.Type == 'prot':
    Type = '1'
elif args.Type == 'nucl':
    Type = '0'
args.p_t = path + 'BLOSUM' + str(args.p_t) + '.txt'
if args.p_R == 'False':
    realign = '0'
elif args.p_R == 'True':
    realign = '1'
if args.p_p == 'False':
    args.p_p = '0'
else:
    args.p_p = '1'
m_args = ''
if args.m_d == 'diags':
    m_args = "-diags"
if args.m_g == 'input':
    m_args = m_args + "-input"
if args.m_c == 'True':
    m_args = m_args + "-clwstrict"

#limit the number of homologs to be drawn by PHI drawer but allow more for other steps
num_hom = args.p_n    
if num_hom > 500:
    num_hom = 500

#-----------Create main output directory-----------------------------------------

if args.o != "Current working directory":
    if args.o[-1] != "/":
        args.o = args.o + "/"
else:
    args.o = ""
#make the main output directory using current date and time
now = datetime.datetime.now()
out_dir = args.o + args.Run_Name + "_" + now.strftime("%Y_%m_%d_%H%M%S")
subp.call(["mkdir", out_dir])


#-----------Make BLAST database -------------------------------------------------

name = args.genome.split('/')
name = name[-1]

#don't remake the database if it's already present
if not (os.path.exists(str(args.genome) + ".nin") or os.path.exists(str(args.genome) + ".nsd") or os.path.exists(str(args.genome) + ".nhr") or os.path.exists(str(args.genome) + ".nsi") or os.path.exists(str(args.genome) + ".nsq")):
    subp.call(["formatdb", "-i", args.genome,"-p", "F", "-o", "T", "-n", args.genome])


#-----------Run custom indexer---------------------------------------------------

#don't rerun the custom indexer if it's already been done
if not os.path.exists(str(args.genome) + ".index"):
    subp.call(["perl", path + "reads_indexer.pl", "-i", args.genome])

if args.b_o == 'None':
    args.b_o = ''


#-----------Single input file----------------------------------------------------

#for single indiviual query or grouped query
if args.q and args.i == 's' or args.i == 'g':
    print "Single input file, single or group input."
    query = args.q
    query_name = os.path.split(args.q)[1]
    #set output directory
    blast_out = os.path.normpath(os.path.join(out_dir, query_name))
    #set output filename
    blast_file_out = os.path.join(blast_out, query_name)
    runTarget(query, blast_out, blast_file_out)
    print "TARGeT has finished!"


#for multiple individual queries-------------------------------------------------
elif args.q and args.i == 'mi':
    print "Single input file, multiple individual inputs."
    basequery = os.path.split(args.q)
    basedir = basequery[0]
    #print str(args.q) + '\n', basedir + '\n'
    subp.call(["python", "split_fasta.py", args.q])
    
    #count the number of new input files
    new_files = glob.glob(str(basedir) + "/*.fix.fa")
    count = 0
    for f in new_files:
        count += 1        
    print count, " files to be processed"

    #counter to keep track of the number of files that have been processed
    p = 0

    #Run pipeline on each file with it's own output directory in the main output directory
    for fasta2 in new_files:
        file_name = os.path.split(fasta2)[1]
        query = fasta2
        #set output directory
        blast_out = os.path.normpath(os.path.join(out_dir, file_name))
        #set output filename
        blast_file_out = os.path.join(blast_out, file_name +".blast")
        runTarget(query, blast_out, blast_file_out)
        p += 1
        print "TARGeT has processed ", p, " of ", count, " subfiles"

    print "TARGeT has finished!"

#-----------Directory input------------------------------------------------------

#for single individual queries or grouped queries
elif args.d and args.i == 's' or args.i == 'g':
    print "Directory input, each file has a single or group input."
    
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
        query_name = f
        blast_out = os.path.normpath(os.path.join(out_dir, query_name)) #output directory
        blast_file_out = os.path.join(blast_out, query_name) #output files basename
        runTarget(query, blast_out, blast_file_out)
        p +=1
        print "TARGeT has processed ", p, " of ", count, " files"
    print "TARGeT has finished!"

#for multiple individual queries-------------------------------------------------
elif args.d and args.i == 'mi':
    print "Directory input, each file has multiple individual queries."
    
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
        basename = os.path.normpath(os.path.join(out_dir, os.path.splitext(f)[0]))
        #make the subdirectory
        subp.call(["mkdir", basename])

        #setup full path to fasta file for splitting
        full = os.path.normpath(os.path.join(args.d, f))
        #split the fasta file
        subp.call(["python", "split_fasta.py", full])

    #count the number of new input files
    files2 = os.listdir(args.d) #get all files in the directory
    count2 = 0
    for fasta2 in new_files:
        if fasta2 in files:
            continue
        else:
            count2 += 1
        #print fasta2        
    print count2, " subfiles to be processed"

    #setup conter for subfiles processed
    p2 = 0

    #Run pipeline on each split file, output directory in subdirectory for pre-split file that's in main output directory
    for fasta2 in new_files:
        #set query to split file path
        query = fasta2
        blast_out = os.path.normpath(os.path.join(basename, os.path.split(os.path.splitext(fasta2)[0])[1]))
        blast_file_out = os.path.normpath(os.path.join(blast_out, os.path.split(fasta2)[1]))
        runTarget(query, blast_out, blast_file_out)
        p2 += 1
        print "TARGeT has processed ", p2, " of ", count2, " subfiles from file ", p + 1, " of ", count, " files"
        p += 1
        print "TARGeT has fully processed ", p, " of ", count, " files"
    print "TARGeT has finished!"
