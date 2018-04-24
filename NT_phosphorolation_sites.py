#S: serine, Y:tyrosine, T:threoine
#Nalani
import random

def hertzStormo(singleton, seqs):
    size = len(singleton)
    profile = createProfile(size, seqs)

def createProfile(numCols, seqs):
     table = []
     row = 0
     # create 2D table initialized with value
     while (row < numRows):
          table.append([])    
          col = 0
          while (col < numCols):
               table[row].append(value)
               col = col + 1
          row = row + 1
     return table

def findPhosSiteMotif(sequences, headers):    
    seqsLen = len(sequences)
    for i in range(seqsLen + 3):
        r = random.randint(0, seqsLen-1)
        s = sequences[r]     
        if len(s) > 3:
            for aminoAcid in s:                
                #phosphorolation sites
                if aminoAcid == 'S' or aminoAcid == 'Y' or aminoAcid == 'T':
                    motif = getMotif(s, aminoAcid)
                    l = len(motif)
                    if l > 4:
                        motifSeqs = []
                        motifHeaders = []
                        for i in range(len(sequences)):
                            if len(sequences[i]) >= l:
                                motifSeqs.append(sequences[i])
                                motifHeaders.append(headers[i])
                        return motif, motifSeqs, motifHeaders    
    return None, None, None

def getMotif(seq, aminoAcid):
    i = seq.index(aminoAcid)
    front = 500
    back = 500
    while ((i-front) < 2) or ((i+back) >= len(seq)):
        front = random.randint(0,5)
        back = random.randint(0,5)
    return seq[i-front:i+back]

def getNumPhosphorolationSites(sequences):
    sAA = []
    yAA = []
    tAA = []
    for s in sequences:
        sI = 0
        yI = 0
        tI = 0
        while sI != -1 and sI+1 < len(s):
            if sI == 0:
                sI = s.find('S', sI)
            else:
                sI = s.find('S', sI+1)
            if sI != -1:
                sAA.append(sI)

        while yI != -1 and yI+1 < len(s):
            if yI == 0:
                yI = s.find('Y', yI)
            else:
                yI = s.find('Y', yI+1)
            if yI != -1:
                yAA.append(yI)

        while tI != -1 and tI+1 < len(s):
            if tI == 0:
                tI = s.find('T', tI)
            else:
                tI = s.find('T', tI+1)
            if tI != -1:
                tAA.append(tI)
    print("SAA: ", sAA)
    print("YAA: ", yAA)
    print("TAA: ", tAA)

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
    headers = []

    # read file for NT ends
    with open(filename) as f:       
        # read in first line
        line = f.readline()        

        while line:
            tempSeq = ""
            i = line.index('(')
            NT_end = int(line[i+1:len(line)-2])
            headers.append(line[:len(line)-1])
            line = f.readline()
            if not line:
                break
            while not isHeader(line):
                tempSeq += line.strip()
                line = f.readline()
                if not line:
                    break
        
            seqs.append(tempSeq[:NT_end])

    return seqs, headers

seqs, headers = convertFileToSequences("new_amino_acids.txt")
#getNumPhosphorolationSites(seqs)
motif, seqs, headers = findPhosSiteMotif(seqs, headers)
print(motif)