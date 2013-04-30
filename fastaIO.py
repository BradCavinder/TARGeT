"""This script contains a fasta file iterator, FastaGeneralIterator, that is a modified version of the FastqGeneralIterator in the QualityIO.py script by Peter Cock that is is part of the Biopython distribution and is governed by its license. Please see the Biopython_License file that is included as part of this package.

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
        title = title.replace("=", "-")
        title = pattern.sub("_", title)
        title = title.replace('___', '_')
        title = title.replace('__', '_')
        title = title.replace('_-_', '-')
        if len(title) >= 99:
            title = title[:99]
        title = title.replace(" ", "")
        yield (title, seq)
            
