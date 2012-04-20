#!/opt/stajichlab/stajichlab-python/2.7.2/bin/python

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
import decimal

#-----------Define frange function for floating point ranges------------------------

def frange(x, y, jump):
    frange_list = []
    while x < y:
        frange_list.append(x)
        x += jump
        continue
    return frange_list
#setup available floating point options for PHI minimum decimal fraction of match length argument
ident_list1 = frange(0.01, 1, 0.01)
#convert list of numbers to comma separated string surrounded by double quatation marks
ident_list2 = []
for i in ident_list1:
    i = str(i)
    i = i[:4]
    ident_list2.append(float(i))
    
path = sys.path[0]
path = str(path) + "/"
#print path

""" The following functions will eventually be finished and implemented to reduce the repetition of code that the script currently uses. The PHI function is already implemented but may be modified further. Also, the creation and implementation of argument groups (one long string for all arguments stored in an object) for each function will allow the proper use of optional arguments that are either present or absent, which are not properly implemented in the current code, as they can be added to the string if present and but won't change it if absent. 

#-----------Define BLASTN function-----------------------------------------------

def BLASTN():
    subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.b_p)])


#-----------Deinfe TBLASTN function----------------------------------------------
def TBLASTN():
    subp.call(["blastall", "-p", "tblastn", "-d", args.genome, "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), + "-b", str(args.b_a), "-v", str(args.b_d) + str(args.b_o), "-a", str(args.b_p)])


#-----------Define BLAST Drawer function-----------------------------------------
def Blast_draw():
    subp.call(["perl", "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])


#-----------Define image converting function-------------------------------------
def img_convert():
    subp.call([" """


#-----------Define PHI function--------------------------------------------------

def PHI(blast_in, PHI_out):
    subp.call(["perl", path + "PHI_2.4.pl", "-i", blast_in, "-q", query, "-D", args.genome, "-o", PHI_out, "-e", str(args.p_e), "-M", str(args.p_M), "-P", Type, "-d", str(args.p_d), "-g", str(args.p_g), "-n", str(args.p_n), "-c", str(args.p_c), "-G", str(args.p_G), "-t", args.p_t, "-f", str(args.p_f), "-p", args.p_p, "-R", realign])

""" . . . more unfinished functions
#-----------Define PHI Drawer function-------------------------------------------

def PHI_draw():
    subp.call(["


#-----------Define MUSCLE function-----------------------------------------------

def MUSCLE():
    subp.call(["muscle", "-in", PHI_out + ".dna", "-out", PHI_out + ".dna.msa", "-maxiters", str(args.m_m), str(m_args)])
"""

#-----------Make command line argument parser------------------------------------

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="This is the command line implementation of TARGeT: Tree Analysis of Related Genes and Transposons (http://target.iplantcollaborative.org/). Please read the README file for program dependencies.")

input_group = parser.add_mutually_exclusive_group(required=True)

input_group.add_argument("-q", metavar="Input file", help="Path to single input query sequence file (may contain multiple sequences of the same type [DNA or protein]).")

input_group.add_argument("-d", metavar="Input directory", help="Path to input directory with multiple query sequence files (each may contain multiple sequences). All sequences must be of the same type (DNA or protein). The directory can not contain non-query files. Include a '/' at the end of the directory name.")

parser.add_argument("-t", dest="Type", metavar="Search type", choices=("nucl", "prot"), default="prot", help="Type of input query sequence(s): DNA = 'nucl'  protein = 'prot'")

parser.add_argument("genome", metavar="Genome sequence file", help="Path of the genome sequence file to be searched (or other sequence to be searched).")

parser.add_argument("Run_Name", help="Name for overall run")

parser.add_argument("-o", metavar="Output path", default="Current working directory", help="Path to an existing directory for output. A subdirectory for the run will be created with further subdirectories containing the output fies for each search sequence.")

parser.add_argument("-i", metavar="Input query type", choices=("s", "mi", "g"), default="s", help="Number of input queries per file(s): s = single query sequence; mi = multiple individual query sequences, g = FASTA multiple sequence alignment or list of sequences; All input files must be the same.")



#BLAST arguments
parser_blast = parser.add_argument_group("BLAST+")

parser_blast.add_argument("-b_e", metavar="BLAST E-value", default="0.001", type=float, help="Maximum Expect value to include a match in the BLAST results")

parser_blast.add_argument("-b_a", metavar="# alignments", default="500", type=int, help="The maximum number of hit alignments to include per query sequence.")

parser_blast.add_argument("-b_d", metavar="# descriptions", default="100", help="The maximum number of descriptions of hits included in BLAST output per query sequence.")

parser_blast.add_argument("-b_o", metavar="Other BLAST options", help="Other options: you may use any other valid Blast+ options for the search type being used except those that change the output format from default (including the use of GenBank accession numbers). Place them all here as you would for running standalone BLASTas one string in single quotes")

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

parser_phi.add_argument("-p_f", metavar="Flanking sequence", type=int, default=0, help="Length of flanking sequences to be included.")

parser_phi.add_argument("-p_p", action='store_true', default="False", help="Filter out psuedogenes (matches with internal stop codons).")



#Muscle arguments
parser_muscle = parser.add_argument_group("Muscle")

parser_muscle.add_argument("-m_m", metavar="Max alignment iterations", type=int, default=1, help="Maximum number of alignment iterations.")

parser_muscle.add_argument("-m_d", metavar="Find diagonals", choices=('diags', ''), default="", help="Find diagonals (faster for similar sequences)")

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
if not os.path.exists(str(args.genome) + ".nin"):
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
    query_name = args.q.split('/')
    #set output directory
    blast_out = os.path.normpath(os.path.join(out_dir, query_name[-1])) 
    #set output filename
    blast_file_out = os.path.join(blast_out, query_name[-1])
    #make output directory
    subp.call(["mkdir", blast_out])

    #use blastn if DNA
    if args.Type == 'nucl':
        print "Using BLASTN"
        subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.b_p)]) #removed args.b_o for now

    #use tblastn if protein
    elif args.Type == 'prot':
        print "Using TBLASTN" 
        subp.call(["blastall", "-p", "tblastn", "-d", args.genome, "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), + "-b", str(args.b_a), "-v", str(args.b_d) + str(args.b_o), "-a", str(args.b_p)])

    #make svg drawing(s)
    print "Making svg image of blast results"
    subp.call(["perl", path + "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])

    #convert svg image to jpg
    print "Converting svg to jpg"
    for svg_file in glob.glob(str(blast_file_out) + "*.svg"): 
        jpg_file = os.path.splitext(svg_file)
        jpg_file = jpg_file[0] + ".jpg"
        subp.call(["convert", svg_file, jpg_file])

    blast_in = str(blast_file_out) + ".blast"
    PHI_out = str(blast_file_out) + ".tcf"
 
    print "Running PHI"
    PHI(blast_in, PHI_out) #call PHI function
    print "PHI finished!"

    #make svg image of PHI homologs
    print "Making svg image of homologs"
    num_hom = args.p_n
    if num_hom > 500:
        num_hom = 500
    subp.call(["perl", path + "PHI_drawer2.pl", "-i", str(PHI_out) + ".list", "-o", str(PHI_out) + ".tcf_drawer", "-m", args.p_t, "-P", Type, "-n", str(num_hom)])

    #convert svg to jpg
    print "Coverting svg image to jpg"
    subp.call(["convert", str(PHI_out) + ".tcf_drawer.svg", str(PHI_out) + ".tcf_drawer.jpg"])

    #Count the number of homologs to determine whether to run MUSCLE and TreeBest and possibly adjust tree image height
    in_file = open(str(PHI_out) + ".dna")
    info = in_file.readlines()
    seq = 0
    for i in info:
        if i[0] == '>':
            seq += 1
    print seq, "homologs"

    if seq > 1:

        #setup Muscle arguments
        m_args = ''
        if args.m_d == 'diags':
            m_args = "-diags"
        if args.m_g == 'input':
            m_args = m_args + "-input"
        if args.m_c == 'True':
            m_args = m_args + "-clwstrict"
    

        #Run Muscle
        print "Running Muscle"
        subp.call(["muscle", "-in", PHI_out + ".dna", "-out", PHI_out + ".dna.msa", "-maxiters", str(args.m_m), str(m_args)])
        print "Muscle Finished, running Treebest"


        #Run TreeBest
        out = open(PHI_out + ".dna.msa.nw", "w") #open output file for redirected stdout
        subp.call(["treebest", "nj", PHI_out + ".dna.msa"], stdout=out)
        out.close() #close output file
        print "TreeBest finished, converting output file to eps image"

        out = open(PHI_out + ".dna.msa.eps", "w") #open output file for redirected stdout
        if seq > 45:
            height = seq * 12
            subp.call(["treebest",  "export", "-y", str(height), PHI_out + ".dna.msa.nw"], stdout=out)
        else:
            subp.call(["treebest",  "export", PHI_out + ".dna.msa.nw"], stdout=out)
        out.close() #close output file

        print "Coverting eps image to jpg"
        subp.call(["convert", PHI_out + ".dna.msa.eps", PHI_out + ".dna.msa.jpg"])
    
    print "TARGeT has finished!"


#for multiple individual queries-------------------------------------------------
elif args.q and args.i == 'mi':
    print "Single input file, multiple individual inputs."
    basequery = os.path.split(args.q)
    basedir = basequery[0]
    print str(args.q) + '\n', basedir + '\n'
    subp.call(["pyfasta", "split", "--header", "%(fasta)s-%(seqid)s.fasta2", args.q])
    
    #count the number of new input files
    new_files = glob.glob(str(basedir) + "*.fasta2")
    count = 0
    for f in new_files:
        count += 1        
    print count, " files to be processed"

    #counter to keep track of the number of files that have been processed
    p = 0

    #Run pipeline on each file with it's own output directory in the main output directory
    for fasta2 in new_files:
        #set output directory
        blast_out = os.path.normpath(os.path.join(out_dir, fasta2)) 
        #set output filename
        blast_file_out = os.path.join(blast_out, fasta2)
        #make output directory
        subp.call(["mkdir", blast_out])
    
        #use blastn if DNA
        if args.Type == 'nucl':
            print "Using BLASTN"
            subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.b_p)]) #removed args.b_o for now

        #use tblastn if protein
        elif args.Type == 'prot':
            print "Using TBLASTN" 
            subp.call(["blastall", "-p", "tblastn", "-d", args.genome, "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), + "-b", str(args.b_a), "-v", str(args.b_d) + str(args.b_o), "-a", str(args.b_p)])

        #make svg drawing(s)
        print "Making svg image of blast results"
        subp.call(["perl", path + "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])

        #convert svg image to jpg
        print "Converting svg to jpg"
        for svg_file in glob.glob(str(blast_file_out) + "*.svg"): 
            jpg_file = os.path.splitext(svg_file)
            jpg_file = jpg_file[0] + ".jpg"
            subp.call(["convert", svg_file, jpg_file])

        blast_in = str(blast_file_out) + ".blast"
        PHI_out = str(blast_file_out) + ".tcf"
 
        print "Running PHI"
        PHI(blast_in, PHI_out) #call PHI function
        print "PHI finished!"

        #make svg image of PHI homologs
        print "Making svg image of homologs"
        num_hom = args.p_n
        if num_hom > 500:
            num_hom = 500
        subp.call(["perl", path + "PHI_drawer2.pl", "-i", str(PHI_out) + ".list", "-o", str(PHI_out) + ".tcf_drawer", "-m", args.p_t, "-P", Type, "-n", str(num_hom)])

        #convert svg to jpg
        print "Coverting svg image to jpg"
        subp.call(["convert", str(PHI_out) + ".tcf_drawer.svg", str(PHI_out) + ".tcf_drawer.jpg"])

        #Count the number of homologs to determine whether to run MUSCLE and TreeBest and possibly adjust tree image height
        in_file = open(str(PHI_out) + ".dna")
        info = in_file.readlines()
        seq = 0
        for i in info:
            if i[0] == '>':
                seq += 1
        print seq, "homologs"

        if seq > 1:
            #setup Muscle arguments
            m_args = ''
            if args.m_d == 'diags':
                m_args = "-diags"
            if args.m_g == 'input':
                m_args = m_args + "-input"
            if args.m_c == 'True':
                m_args = m_args + "-clwstrict"
    
            #Run Muscle
            print "Running Muscle"
            subp.call(["muscle", "-in", PHI_out + ".dna", "-out", PHI_out + ".dna.msa", "-maxiters", str(args.m_m), str(m_args)])
            print "Muscle Finished, running Treebest"


            #Run TreeBest
            out = open(PHI_out + ".dna.msa.nw", "w") #open output file for redirected stdout
            subp.call(["treebest", "nj", PHI_out + ".dna.msa"], stdout=out)
            out.close() #close output file
            print "TreeBest finished, converting output file to eps image"

            out = open(PHI_out + ".dna.msa.eps", "w") #open output file for redirected stdout
            if seq > 45:
                height = seq * 12
                subp.call(["treebest",  "export", "-y", str(height), PHI_out + ".dna.msa.nw"], stdout=out)
            else:
                subp.call(["treebest",  "export", PHI_out + ".dna.msa.nw"], stdout=out)
            out.close() #close output file

            print "Coverting eps image to jpg"
            subp.call(["convert", PHI_out + ".dna.msa.eps", PHI_out + ".dna.msa.jpg"])
    
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

        subp.call(["mkdir", blast_out]) #make output directory

        #use blastn if DNA
        if args.Type == 'nucl':
            print "Using BLASTN"
            subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-b", str(args.b_d), "-a", str(args.b_p)]) #removed args.b_o for now

        #use tblastn if protein
        elif args.Type == 'prot':
            print "Using TBLASTN" 
            subp.call(["tblastn", "-db", args.genome, "-query", query, "-out", blast_out + ".blast", "-evalue", str(args.b_e), + "-num_alignments", str(args.b_a), "-num_descriptions", str(args.b_d), "-a", str(args.b_p)]) #removed args.b_o for now

        #make svg drawing(s)
        print "Making svg image of blast results"
        subp.call(["perl", path + "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])

        #convert all svg images in output folder to jpg
        print "Converting svg to jpg"
        for svg_file in glob.glob(str(blast_file_out) + "*.svg"): 
            jpg_file = os.path.splitext(svg_file)
            jpg_file = jpg_file[0] + ".jpg"
            subp.call(["convert", svg_file, jpg_file])

        #setup PHI input and output files
        blast_in = str(blast_file_out) + ".blast"
        PHI_out = str(blast_file_out) + ".tcf" 

        #call PHI function
        print "Running PHI"
        PHI(blast_in, PHI_out) 
        print "PHI finished!"

        #make svg image of PHI homologs
        print "Making svg image of homologs"
        num_hom = args.p_n
        if num_hom > 500:
            num_hom = 500
        subp.call(["perl", path + "PHI_drawer2.pl", "-i", str(PHI_out) + ".list", "-o", str(PHI_out) + ".tcf_drawer", "-m", args.p_t, "-P", Type, "-n", str(num_hom)])

        #convert svg to jpg
        print "Converting svg image to jpg"
        subp.call(["convert", str(PHI_out) + ".tcf_drawer.svg", str(PHI_out) + ".tcf_drawer.jpg"])

        #Count the number of homologs to determine whether to run MUSCLE and TreeBest and possibly adjust tree image height
        in_file = open(str(PHI_out) + ".dna")
        info = in_file.readlines()
        seq = 0
        for i in info:
            if i[0] == '>':
                seq += 1
        print seq, "homologs"

        if seq > 1:
            #setup Muscle arguments
            m_args = ''
            if args.m_d == 'diags':
                m_args = "-diags"
            if args.m_g == 'input':
                m_args = m_args + "-input"
            if args.m_c == 'True':
                m_args = m_args + "-clwstrict"
        
            #call Muscle
            print "Running Muscle"
            subp.call(["muscle", "-in", PHI_out + ".dna", "-out", PHI_out + ".dna.msa", "-maxiters", str(args.m_m), str(m_args)])
            print "Muscle Finished, running Treebest"

            #call TreeBest
            out = open(PHI_out + ".dna.msa.nw", "w") #open output file for redirected stdout
            subp.call(["treebest", "nj", PHI_out + ".dna.msa"], stdout=out)
            out.close() #close output file
            print "TreeBest finished, converting output file to eps image"

            out = open(PHI_out + ".dna.msa.eps", "w") #open output file for redirected stdout
            if seq > 45: #adjust tree image height if over 45 homologs
                height = seq * 12
                subp.call(["treebest",  "export", "-y", str(height), PHI_out + ".dna.msa.nw"], stdout=out)
            else: #use default tree image height if 45 or fewer homologs
                subp.call(["treebest",  "export", PHI_out + ".dna.msa.nw"], stdout=out)
            out.close() #close the output file

            print "Coverting eps image to jpg"
            subp .call(["convert", PHI_out + ".dna.msa.eps", PHI_out + ".dna.msa.jpg"])
        
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
        subp.call(["pyfasta", "split", "--header", "%(fasta)s-%(seqid)s.fasta2", full])

        #count the number of new input files
        new_files = glob.glob(full + "*.fasta2")
        count2 = 0
        for fasta2 in new_files:
            count2 += 1
            print fasta2        
        print count2, " subfiles to be processed"

        #setup conter for subfiles processed
        p2 = 0

        #Run pipeline on each split file, output directory in subdirectory for pre-split file that's in main output directory
        for fasta2 in new_files:
            #set query to split file path
            query = fasta2

            #set split file directory
            blast_out = os.path.normpath(os.path.join(basename, os.path.split(os.path.splitext(fasta2)[0])[1])) 
            #make output directory
            subp.call(["mkdir", blast_out])
 
            #set output filename
            blast_file_out = os.path.normpath(os.path.join(blast_out, os.path.split(fasta2)[1]))
            
            
            #use blastn if DNA
            if args.Type == 'nucl':
                print "Using BLASTN"
                subp.call(["blastall", "-p", "blastn", "-d", str(args.genome), "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), "-b", str(args.b_a), "-v", str(args.b_d), "-a", str(args.b_p)]) #removed args.b_o for now

            #use tblastn if protein
            elif args.Type == 'prot':
                print "Using TBLASTN" 
                subp.call(["blastall", "-p", "tblastn", "-d", args.genome, "-i", query, "-o", str(blast_file_out) + ".blast", "-e", str(args.b_e), + "-b", str(args.b_a), "-v", str(args.b_d) + str(args.b_o), "-a", str(args.b_p)])

            #make svg drawing(s)
            print "Making svg image of blast results"
            subp.call(["perl", path + "v3_blast_drawer.pl", "-i", str(blast_file_out) + ".blast", "-o", str(blast_file_out)])

            #convert svg image to jpg
            print "Converting svg to jpg"
            for svg_file in glob.glob(str(blast_file_out) + "*.svg"): 
                jpg_file = os.path.splitext(svg_file)
                jpg_file = jpg_file[0] + ".jpg"
                subp.call(["convert", svg_file, jpg_file])

            blast_in = str(blast_file_out) + ".blast"
            PHI_out = str(blast_file_out) + ".tcf"
 
            print "Running PHI"
            PHI(blast_in, PHI_out) #call PHI function
            print "PHI finished!"

            #make svg image of PHI homologs
            print "Making svg image of homologs"
            num_hom = args.p_n
            if num_hom > 500:
                num_hom = 500
            subp.call(["perl", path + "PHI_drawer2.pl", "-i", str(PHI_out) + ".list", "-o", str(PHI_out) + ".tcf_drawer", "-m", args.p_t, "-P", Type, "-n", str(num_hom)])


            #convert svg to jpg
            print "Coverting svg image to jpg"
            subp.call(["convert", str(PHI_out) + ".tcf_drawer.svg", str(PHI_out) + ".tcf_drawer.jpg"])

            #Count the number of homologs to determine whether to run MUSCLE and TreeBest and possibly adjust tree image height
            in_file = open(str(PHI_out) + ".dna")
            info = in_file.readlines()
            seq = 0
            for i in info:
                if i[0] == '>':
                    seq += 1
            print seq, "homologs"

            if seq > 1:
                #setup Muscle arguments
                m_args = ''
                if args.m_d == 'diags':
                    m_args = "-diags"
                if args.m_g == 'input':
                    m_args = m_args + "-input"
                if args.m_c == 'True':
                    m_args = m_args + "-clwstrict"
    

                #Run Muscle
                print "Running Muscle"
                subp.call(["muscle", "-in", PHI_out + ".dna", "-out", PHI_out + ".dna.msa", "-maxiters", str(args.m_m), str(m_args)])
                print "Muscle Finished, running Treebest"

                #Run TreeBest
                out = open(PHI_out + ".dna.msa.nw", "w") #open output file for redirected stdout
                subp.call(["treebest", "nj", PHI_out + ".dna.msa"], stdout=out)
                out.close() #close output file
                print "TreeBest finished, converting output file to eps image"

                out = open(PHI_out + ".dna.msa.eps", "w") #open output file for redirected stdout
                if seq > 45:
                    height = seq * 12
                    subp.call(["treebest",  "export", "-y", str(height), PHI_out + ".dna.msa.nw"], stdout=out)
                else:
                    subp.call(["treebest",  "export", PHI_out + ".dna.msa.nw"], stdout=out)
                out.close() #close output file

                print "Coverting eps image to jpg"
                subp.call(["convert", PHI_out + ".dna.msa.eps", PHI_out + ".dna.msa.jpg"])
            
            p2 += 1
            print "TARGeT has processed ", p2, " of ", count2, " subfiles from file ", p + 1, " of ", count, " files"
        
        p += 1

        print "TARGeT has fully processed ", p, " of ", count, " files"

    print "TARGeT has finished!"
