#S: serine, Y:tyrosine, T:threoine
#Nalani

def findPhosSites(sequence):
    return

def isHeader(line):
    if (line[0] == '>'):
        return True
    else:
        return False

#######################################################
# convertFileToSequences - takes a FASTA file with many sequences and returns the
# DNA sequences in an array
#######################################################
def convertFileToSequences(filename):
    #all NT of sequences
    seqs = []

    # read file for NT ends
    with open(filename) as f:       
        # read in first line
        line = f.readline()        

        while line:
            tempSeq = ""
            i = line.index('(')
            NT_end = int(line[i+1:len(line)-2])
            line = f.readline()
            if not line:
                break
            while not isHeader(line):
                tempSeq += line.strip()
                line = f.readline()
                if not line:
                    break
        
            seqs.append(tempSeq[:NT_end])

    print(seqs)

    return seqs

convertFileToSequences("test_file_parse.txt")