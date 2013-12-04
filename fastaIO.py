"""This script contains several function for working with fasta formated sequence files. The main component is a fasta file iterator, FastaGeneralIterator, that is a modified version of the FastqGeneralIterator in the QualityIO.py script by Peter Cock that is is part of the Biopython distribution and is governed by its license. Please see the Biopython_License file that is included as part of this package for detail.

There is also code, FastaTitleStandardization, to standardize fasta record title lines output from third party software.

FastqGeneralIterator was modified to handle fasta files and FastaTitleStandardization was written by Brad Cavinder""" 


def FastaGeneralIterator(handle):
    """Iterate over Fasta records as string tuples (not as Biopython SeqRecord objects).
    """

    #Skip any text before the first record (e.g. blank lines, comments?)
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0]!=">":
            raise ValueError("Records in Fasta files should start with '>' character")
        title_line = line[1:].rstrip()
        #Will now be at least one line of sequence data - string concatenation is used (if needed)
        #rather using than the "".join(...) trick just in case it is multiline:
        seq_string = handle.readline().rstrip()
        #There may now be more sequence lines, or the ">" for the next record:
        while True:
            line = handle.readline()
            if not line:
                break
            if len(line) == 0:
                continue
            if line[0] == ">":               
                break
            seq_string += line.rstrip() #removes trailing newlines
        #Assuming whitespace isn't allowed, any should be caught here:
        if " " in seq_string or "\t" in seq_string:
            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)

        #Return the record and then continue...
        yield (title_line, seq_string)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"
    
def FastaTitleStandardization(handle):
    """Use FastaGeneralIterator to iterate through records in a fasta file while standardizing the title of each record (third party programs are inconsitant to each other in title line naming), removing whitespace and special characters other than '-'. Also, ':' and '=' are changed to '-'. Title lines longer than 60 characters wil be shortened to 60.
    """
    import re
    pattern = re.compile('[^\w-]')
    for title, seq in FastaGeneralIterator(handle):
        #Remove whitespace and special characters except '-' from title and shorten it to no more than 100 characters
        title = title.replace(":", "-")
        title = title.replace("(", "-")
        title = title.replace(")", "-")
        title = title.replace("=", "-")
        title = pattern.sub("_", title)
        title = title.replace('___', '_')
        title = title.replace('__', '_')
        title = title.replace('_-_', '-')
        if len(title) >= 99:
            title = title[:99]
        title = title.replace(" ", "")
        yield (title, seq)

def median(in_list):
    #if list has even number of elements, take the avg of middle two
    #otherwise return middle elemenent
    in_list.sort()
    mid = len(in_list)/2
    if len(alist) % 2 == 0:  
        return (srtd[mid-1] + srtd[mid]) / 2.0
    else:
        return in_list[mid]            

def reverse_complement(seq):
    import string
    trans_table = string.maketrans("ATGCatgcNn", "TACGtacgNn")
    rev_seq = seq[::-1]
    rev_comp = rev_seq.translate(trans_table)
    return(rev_comp)
    
def individual_seq_len (in_path):
    in_handle = open(in_path, 'r')
    out_path = in_path + ".lengths"
    out_handle = open(out_path, 'w')

    for title, seq in FastaGeneralIterator(in_handle):
        seq_len = len(seq)
        print>>out_handle, title, "\t", seq_len

    in_handle.close()
    out_handle.close()
    
def total_seq_len (in_path, write_mode, out_path, species_name):
    
    in_handle = open(in_path, "r")
    seq_len = 0
    for title, seq in FastaGeneralIterator(in_handle):
        seq_len +=  len(seq)
            
    in_handle.close()
    if write_mode == "append" or write_mode == "a" or write_mode == "A":
        out_handle = open(out_path, "a")
        print>>out_handle, "\t".join([species_name , str(seq_len)])
    else:
        out_handle = open(out_path, "w")
        print>>out_handle, "species\tgenome_length\n" + species_name + "\t" + str(seq_len)
        
    out_handle.close()

def sequence_retriever(contig, start, end, flank, genome_dict3):
    needed_left = 0
    needed_right = 0
    wanted_seq = ''
    add_left = ''
    add_right = ''
    left_coord = ''
    right_coord = ''
    if contig in genome_dict3:
        seq = genome_dict3[contig]
        contig_seq_len = len(seq)
        if int(flank) < int(start):
            left_coord = (int(start)-int(flank))
        else:
            needed_left = (int(flank) - int(start))
            left_coord = 0
        if (contig_seq_len - int(flank)) >= int(end):
            right_coord = (int(end+1) + int(flank))
        else:
            needed_right = int(end) - ((contig_seq_len - int(flank)))
            right_coord = contig_seq_len
        if needed_left > 0:
            add_left = "N" * needed_left
        if needed_right > 0:
            add_right = "N" * needed_right
        wanted_seq = add_left + seq[left_coord:right_coord] + add_right
    else:
        print "Didn't find " + contig + " in sequence dictionary."
    return(wanted_seq)

def SplitLongString(s, w):
    temp = ''.join(s[x:x+w] + '\n' for x in xrange(0, len(s), w))
    temp = temp.rstrip("\n")
    return temp

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value
