#S: serine, Y:tyrosine, T:threoine
#Nalani
import random
import math

"""

"""

def hertzStormo(aminoAcid, singleton, singletonSeq, singletonI, seqs, headers):
    size = len(singleton)
    nMerSet = []
    nMerSet.append(singleton)
    #get initial profile for the chosen singleton
    profile = createProfile(size, nMerSet)

    bestIC = 0.0
    bestSet = [singleton]
    bestHeaders = [headers[singletonSeq]]
    bestLocations = [singletonI]

    doneSeqs = [singletonSeq]

    mini = math.floor(float(len(bestSet[0]))/4.0)+1

    for x in range(len(seqs)-1):
        addIC, addSet, addHeaders, addLocations, addSeq = iteration(seqs, headers, bestSet, doneSeqs, aminoAcid, size)
        same = 0
        for c in range(len(addSet)):
            if addSet[c] == bestSet[0][c]:
                same += 1
        if same >= mini:
            bestHeaders.append(addHeaders)
            bestIC = addIC
            bestLocations.append(addLocations)
            bestSet.append(addSet)
            doneSeqs.append(addSeq)

    return bestIC, bestSet, bestHeaders

def iteration(seqs, headers, currentSet, doneSeqs, aminoAcid, size):    
    bestNMer = ""
    bestHeader = ""
    bestLocation = ""
    bestIC = 0.0
    bestSeq = 0
    for s in range(len(seqs)):
        if s not in doneSeqs:
            i = 0                
            while i != -1:
                i = seqs[s].find(aminoAcid, i+1)
                if i != -1:
                    for x in range(size):
                        if i-x >= 0 and i+size-x < len(seqs[s]):
                            temp = []
                            for c in currentSet:
                                temp.append(c)
                            temp.append(seqs[s][i-x:i+size-x])
                            tempProfile = createProfile(size, temp)
                            val = getIC(tempProfile)
                            if val > bestIC:
                                bestIC = val
                                bestNMer = seqs[s][i-x:i+size-x]
                                bestHeader = headers[s]                            
                                bestLocation = i-x
                                bestSeq = s
    return bestIC, bestNMer, bestHeader, bestLocation, bestSeq

def getIC(profile):
    #background = [.085, .095, .07, .01, .06, .095, .06, .02, .005, .04, .05, .06, .03, .03, .04, .06, .05, .02, .06, .06]    
    background = [.05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]    
    i = 0
    for col in range(len(profile[0])):        
        for row in range(len(profile)):
            check = profile[row][col]
            if check > 0.0:
                i += check * math.log(check/background[row])
    return i
                    
def createProfile(numCols, seqs):
    profile = []  
    for col in range(numCols):            
        avgs = getAvg(seqs, col)
        profile.append(avgs)
    return profile

def getAvg(seqs, col):
    avgs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for s in seqs:
        aminoAcid = s[col]
        if aminoAcid == 'G':
            avgs[0] += 1
        elif aminoAcid == 'A':
            avgs[1] += 1
        elif aminoAcid == 'V':
            avgs[2] += 1
        elif aminoAcid == 'C':
            avgs[3] += 1
        elif aminoAcid == 'P':
            avgs[4] += 1
        elif aminoAcid == 'L':
            avgs[5] += 1
        elif aminoAcid == 'I':
            avgs[6] += 1
        elif aminoAcid == 'M':
            avgs[7] += 1
        elif aminoAcid == 'W':
            avgs[8] += 1
        elif aminoAcid == 'F':
            avgs[9] += 1
        elif aminoAcid == 'S':
            avgs[10] += 1
        elif aminoAcid == 'T':
            avgs[11] += 1
        elif aminoAcid == 'Y':
            avgs[12] += 1
        elif aminoAcid == 'N':
            avgs[13] += 1
        elif aminoAcid == 'Q':
            avgs[14] += 1
        elif aminoAcid == 'K':
            avgs[15] += 1
        elif aminoAcid == 'R':
            avgs[16] += 1
        elif aminoAcid == 'H':
            avgs[17] += 1
        elif aminoAcid == 'D':
            avgs[18] += 1
        elif aminoAcid == 'E':
            avgs[19] += 1
        else:
            print(aminoAcid, " is not accounted for")

    for i in range(20):
        avgs[i] = avgs[i]/len(seqs)
    return avgs

def findPhosSiteMotif(sequences, headers, longest, shortest):    
    seqsLen = len(sequences)
    for i in range(seqsLen + 3):
        r = random.randint(0, seqsLen-1)
        s = sequences[r]     
        if len(s) > 3:
            for aminoAcid in range(len(s)):                
                #phosphorolation sites
                if s[aminoAcid] == 'S' or s[aminoAcid] == 'Y' or s[aminoAcid] == 'T':
                    motif = getMotif(s, s[aminoAcid], longest, shortest)
                    l = len(motif)
                    if l > 4:
                        motifSeqs = []
                        motifHeaders = []
                        motifSeq = 0
                        unaccounted = 0
                        for x in range(seqsLen):
                            if len(sequences[x]) >= l:                                                                                                                                                            
                                motifSeqs.append(sequences[x])
                                motifHeaders.append(headers[x])
                                if s == sequences[x]:                                    
                                    motifSeq = x-unaccounted                                
                            else:
                                unaccounted += 1
                        return s[aminoAcid], motif, motifSeq, aminoAcid, motifSeqs, motifHeaders    
    return None, None, None

def getMotif(seq, aminoAcid, l , s):
    i = seq.index(aminoAcid)
    front = 500
    back = 500
    while ((i-front) < 2) or ((i+back) >= len(seq)) or len(seq[i-front:i+back]) < s:
        front = random.randint(0,l)
        back = random.randint(0,l)
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
            if i == -1:
                print("Sorry, file is in wrong format. Exiting now.")
                return
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

def setHeader(filename):
    with open(filename, "w") as f:
        f.write("Output for NT phosphorolation sites.\n")

def writeData(filename, sets, headers):
    with open(filename, 'a') as f:
        f.write("\n----------- motif "+sets[0]+"-----------\n")
        for i in range(len(sets)):
            f.write(headers[i]+"\n"+sets[i]+"\n")

inputfile = input("What data file would you like to use? ")
outputfile = "output_NT_phos_" + inputfile
longestMotifSize = int(input("What is the longest motif length you want? "))
shortestMotifSize = int(input("What is the shortest motif length you want? "))
numberOfMatches = int(input("What is the least number of occurances you are intested in for each motif? "))
while longestMotifSize < shortestMotifSize or longestMotifSize < 1 or shortestMotifSize < 1 or numberOfMatches < 1:
    print("\nSorry, the last three inputs don't make sense. Try again.\n")
    longestMotifSize = int(input("What is the longest motif length you want? "))
    shortestMotifSize = int(input("What is the shortest motif length you want? "))
    numberOfMatches = int(input("What is the least number of occurances you are intested in for each motif? "))

setHeader(outputfile)
seqs, headers = convertFileToSequences(inputfile)
motifsSoFar = []
for x in range(10):
    aminoAcid, motif, motifSeq, motifI, seqs, headers = findPhosSiteMotif(seqs, headers, longestMotifSize, shortestMotifSize)
    if motif not in motifsSoFar:
        bestIC, bestSet, bestHeaders = hertzStormo(aminoAcid, motif, motifSeq, motifI, seqs, headers)
        if len(bestSet) >= numberOfMatches:
            writeData(outputfile, bestSet, bestHeaders)
        motifsSoFar.append(motif)

print("\nDone finding motifs. To see your results, please open " + outputfile)