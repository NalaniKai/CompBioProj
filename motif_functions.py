#!/usr/bin/python

######################################################
# Expectation-Maximization Algorithm for Finding Motifs
# Customizable # of sequences and motif lengths
# Spring 2018
# Trenton Langer
#######################################################

import random
import math
import time
import sys

##################################################################
# generateSeqsFromFasta - takes a list of FASTA organized sequences
# and store the sequences in an array
##################################################################
def generateSeqsFromFasta(filename):
    # read in file
    file = open(filename)

    # read in entire file
    allLines = file.readlines()
    seqs = []
    count = -1
  
    for i in range(len(allLines)):
        if(allLines[i][0] == '>'):
            seqs.append("")
            count += 1
        else:
            seqs[count] += allLines[i]
    
    for i in range(len(seqs)):
        seqs[i] = seqs[i].replace("\r", "")
        seqs[i] = seqs[i].replace("\n", "")
        
    # close file
    file.close()
    
    return seqs

##########################################################
# Suppose *sequence* is an Amino Acid sequence, and suppose
# profile is organised such as:
# profile[position in motif][amino acid], where
# profile[0][0] would be 'F1', profile[0][1] would be 'L1',
# profile[1][1] would be 'L2', etc...
# Order of Amino Acids is "FLIMVSPTAYHQNKDECWRG"
###########################################################
def getBestMatchInAmino(sequence, profile, motifLength):
    seq = sequence
    bestMatch = 0
    bestMatchIndex = 0
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']

    for i in range(len(sequence)-(motifLength+1)):
        tempMatch = 1
        
        for j in range(motifLength):
            for k in range(len(profile[j])):
                if seq[i+j] == aminoAcids[k]:
                    tempMatch = tempMatch*profile[j][k]
            
        if tempMatch > bestMatch:
            bestMatch = tempMatch
            bestMatchIndex = i
    
    return seq[bestMatchIndex:bestMatchIndex+motifLength]
    
###########################################################
# Calculates and returns the frequency of the given
# amino acid, aa, at given position in list of sequences
# If frequency is 0, return .001 (so info content calculation is
# well-defined)
###########################################################
def getFrequencyOfAminoAtPosition(aa, position, sequences):
    count = 0.0
    numSeqs = len(sequences)
    
    for i in range(numSeqs):
        if sequences[i][position] == aa:
            count += 1
        
    rtnVal = count/numSeqs
    if rtnVal == 0:
        rtnVal += .001
    return rtnVal
    
###########################################################
# Calculates and returns the information content of
# a motif profile, given the values in the profile
# Assumes background frequency of 100/20 = 5% for each nucleotide
###########################################################
def calcInfoContentAmino(profile):
    # To calculate log base 2, use math.log(x, 2)
    sum = 0.0
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(len(profile)):    #uses length of F row of profile table
        for j in range(len(aminoAcids)):
            profileVal = profile[i][j]
            sum += profileVal*math.log(profileVal/0.05, 2)
    
    return sum
    
#########################################################
# Prints the specified motif model to the screen.
#               1       2       3       4       ...
#       F       F1      F2      F3      F4      ...
#       L       L1      L2      L3      L4      ...
#       I       I1      I2      I3      I4      ...
#       ...     ...     ...     ...     ...     ...
#########################################################
def printMotifAmino(profile, motifLength):
    line0 = ""
    lines = []
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(len(aminoAcids)):
        lines.append(aminoAcids[i])

    for i in range(motifLength):
        line0 = line0 + "\t" + str(i)
        for j in range(len(lines)):
            lines[j] += "\t" + str("{0:.3f}".format(profile[i][j]))
    
    print(line0)
    for i in range(len(lines)):
        print(lines[i])
    print()
            
##########################################################
# Uses expectation-maximization algorithm for finding
# the best motif in the any number of aa sequences
##########################################################
def findMotifAmino(sequences, motifLength):

    numSeqs = len(sequences)
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    
    # remember the motif instances from the previous iteration so we know when algorithm converges
    old_instance = []
    randomStart = []
    instance = []
    for i in range(numSeqs):
        old_instance.append("")
        randomStart.append(random.randint(0, len(sequences[i])-motifLength))
        #print("randStart = " + str(randomStart[i]))
        instance.append(sequences[i][randomStart[i]:randomStart[i]+motifLength])
    
    # repeat two steps of EM until convergence.
    convergence = False
    while(not convergence):    
        #temp seq list to test getfrequency
        seqs = []
        for i in range(numSeqs):
            seqs.append(instance[i])
        
        profile = []
        for i in range(motifLength):
            profile.append([])
            for j in range(len(aminoAcids)):
                profile[i].append(getFrequencyOfAminoAtPosition(aminoAcids[j], i, seqs))

        # print the motif model to the screen (comment out when running findBestMotif)
        #printMotifAmino(profile, motifLength)
        
        convergence = True
        for i in range(numSeqs):
            old_instance[i] = instance[i]   # re-assign old instances as the current instances to determine convergence
            instance[i] = getBestMatchInSequence(sequences[i], profile, motifLength)    # find best match in each sequence (step 2 of EM) 
            if(old_instance[i] != instance[i]):
                convergence = False
        
    return profile
    
##################################################################
# Scans through profile and returns list of consensus sequences
##################################################################
def generateConsensusAmino(profile):
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    consensusSeqs = []
    consensusSeqs.append("")
    for i in range(len(profile)):
        maxVal = max(profile[i])
        count = 0
        oldNumSeqs = 1
        for j in range(len(profile[i])):
            if profile[i][j] == maxVal:
                if count == 0:
                    for k in range(len(consensusSeqs)):
                        consensusSeqs[k] += aminoAcids[j]
                    count += 1
                    oldNumSeqs = len(consensusSeqs)
                else:
                    for k in range(oldNumSeqs):
                        cutoff = len(consensusSeqs[k]) - 1
                        consensusSeqs.append(consensusSeqs[k][0:cutoff] + aminoAcids[j])
                        count += 1
                    
        
    #print(consensusSeqs)
    #print(len(consensusSeqs))
    return consensusSeqs      
    
##################################################################
# Scans through profile and returns list of consensus sequences
##################################################################
def consensusMatch(sequences, consensusSeqs):
    matchScores = []
    for i in range(len(sequences)):
        bestScore = -10000
        for j in range(len(consensusSeqs)):
            score = blosumScore(sequences[i], consensusSeqs[j])
            if score > bestScore:
                bestScore = score
        matchScores.append(bestScore)   
    #print(matchScores)
    return matchScores

##################################################################
# Blosum Compare two AA sequences
##################################################################
def blosumScore(seq1, seq2):
    score = 0
    
    return score

##################################################################
# Runs findMotifAmino "repetition" number of times and returns the motif
# with the highest information content
# When calling this function, be sure to comment out
# printing in the findMotif function
##################################################################
def findBestMotifAmino(sequences, motifLength, repetition, seqPrinting):
    bestInfoContent = -10000000
    
    bestProfile = []
    instance = []
    
    totalTime = 0
    last10Time = 0
    for i in range(repetition):
        time1 = time.clock()
        m = findMotifAmino(sequences, motifLength)
        infoContent = calcInfoContentAmino(m)
        #print(infoContent)

        # keep best motif found so far
        if infoContent > bestInfoContent:
                bestProfile = m
                bestInfoContent = infoContent
        time2 = time.clock()
        totalTime += (time2-time1)
        last10Time += (time2-time1)
        if(i % 10 == 0 and True):
            #print(str(i) + "\t" + str(last10Time))
            progress = float(i*100 / repetition)
            etr = last10Time*(repetition/10)*(repetition-i)/(repetition)
            if i != 0:
                sys.stdout.write("Run Progress: %f%%, Estimated Time Remaining: %f seconds \r" % (progress, etr))
                sys.stdout.flush()
            else:
                sys.stdout.write("Run Progress: %f%% \r" % progress)
                sys.stdout.flush()
            last10Time = 0
        
	
    # display information about best motif
    m = bestProfile
    
    # print best motif
    printMotifAmino(m, motifLength)
    
    #generateConsensusAmino sequence
    consensusSeqs = generateConsensusAmino(m)
    
    # find best match in each sequence   
    for i in range(len(sequences)):
        instance.append(getBestMatchInAmino(sequences[i], m, motifLength))
        if(seqPrinting):
            print(sequences[i] + "\t" + instance[i])

    
    #get best consensus match score for each instance
    matchScores = consensusMatch(instance, consensusSeqs)
    
    # print runtime
    print("Total Execution Time: " + "{0:.2f}".format(totalTime) + " seconds")
    
    # print info content
    print ("\nInformation content: " + str(bestInfoContent)+ "\n")
        
    for i in range(len(instance)):
        print(str(i) + "\t" + instance[i] + ", match score: " + str(matchScores[i]))
        
    return instance
    
##########################################################
# Suppose *sequence* is a DNA sequence, and suppose
# profile is organised such as:
# profile[position in motif][nucleotide], where
# profile[0][0] would be 'A1', profile[0][1] would be 'C1',
# profile[1][1] would be 'C2', etc...
# Order of Nucleotides is 'ACGT'
###########################################################
def getBestMatchInSequence(sequence, profile, motifLength):
    seq = sequence
    bestMatch = 0
    bestMatchIndex = 0

    for i in range(len(sequence)-(motifLength+1)):
        tempMatch = 1
        
        for j in range(motifLength):
            if seq[i+j] == 'A':
                tempMatch = tempMatch*profile[j][0] #value for 'A' at position 'j'
            elif seq[i+j] == 'C':
                tempMatch = tempMatch*profile[j][1] #value for 'C' at position 'j'
            elif seq[i+j] == 'G':
                tempMatch = tempMatch*profile[j][2] #value for 'G' at position 'j'
            elif seq[i+j] == 'T' or seq[i+j] == 'U':
                tempMatch = tempMatch*profile[j][3] #value for 'T' at position 'j'
            
        if tempMatch > bestMatch:
            bestMatch = tempMatch
            bestMatchIndex = i
    
    return seq[bestMatchIndex:bestMatchIndex+motifLength]

###########################################################
# Calculates and returns the frequency of the given
# nucleotide, nt, at given position in list of sequences
# If frequency is 0, return .001 (so info content calculation is
# well-defined)
###########################################################
def getFrequencyOfNucleotideAtPosition(nt, position, sequences, amino):
    count = 0.0
    numSeqs = len(sequences)
    if(amino):
        if(nt == 'T'):
            nt = 'U'
    for i in range(numSeqs):
        if sequences[i][position] == nt:
            count += 1
        
    rtnVal = count/numSeqs
    if rtnVal == 0:
        rtnVal += .001
    return rtnVal

###########################################################
# Calculates and returns the information content of
# a motif profile, given the values in the profile
# Assumes background frequency of 25% for each nucleotide
###########################################################
def calcInfoContent(profile):
    # To calculate log base 2, use math.log(x, 2)
    sum = 0.0
    for i in range(len(profile)):    #uses length of A row of profile table
        A = profile[i][0]
        C = profile[i][1]
        G = profile[i][2]
        T = profile[i][3]
        sum += A*math.log(A/0.25, 2) + C*math.log(C/0.25, 2) + G*math.log(G/0.25, 2) + T*math.log(T/0.25, 2)
    
    return sum
    
#########################################################
# Prints the specified motif model to the screen.
#               1       2       3       4       ...
#       A       A1      A2      A3      A4      ...
#       C       C1      C2      C3      C4      ...
#       G       G1      G2      G3      G4      ...
#       T       T1      T2      T3      T4      ...
#########################################################
def printMotif(profile, motifLength):
    line0 = ""
    line1 = "A"
    line2 = "C"
    line3 = "G"
    line4 = "T"

    for i in range(motifLength):
        line0 = line0 + "\t" + str(i)
        line1 = line1 + "\t" + str("{0:.3f}".format(profile[i][0]))
        line2 = line2 + "\t" + str("{0:.3f}".format(profile[i][1]))
        line3 = line3 + "\t" + str("{0:.3f}".format(profile[i][2]))
        line4 = line4 + "\t" + str("{0:.3f}".format(profile[i][3]))
    
    print(line0)
    print(line1)
    print(line2)
    print(line3)
    print(line4)

##########################################################
# Uses expectation-maximization algorithm for finding
# the best motif in the any number of sequences
# amino - boolean for whether giving AA or nucleotide sequences
##########################################################
def findMotif(sequences, motifLength, amino):

    numSeqs = len(sequences)
    # remember the motif instances from the previous iteration so we know when algorithm converges
    old_instance = []
    randomStart = []
    instance = []
    for i in range(numSeqs):
        old_instance.append("")
        randomStart.append(random.randint(0, len(sequences[i])-motifLength))
        if(amino):
            while(randomStart[i] % 3 != 0):
                randomStart[i] = random.randint(0, len(sequences[i])-motifLength)
        #print("randStart = " + str(randomStart[i]))
        instance.append(sequences[i][randomStart[i]:randomStart[i]+motifLength])
    
    # repeat two steps of EM until convergence.
    convergence = False
    while(not convergence):    
        #temp seq list to test getfrequency
        seqs = []
        for i in range(numSeqs):
            seqs.append(instance[i])
        
        profile = []
        for i in range(motifLength):
            profile.append([])
            profile[i].append(getFrequencyOfNucleotideAtPosition("A", i, seqs, amino))
            profile[i].append(getFrequencyOfNucleotideAtPosition("C", i, seqs, amino))
            profile[i].append(getFrequencyOfNucleotideAtPosition("G", i, seqs, amino))
            profile[i].append(getFrequencyOfNucleotideAtPosition("T", i, seqs, amino))

        # print the motif model to the screen (comment out when running findBestMotif)
        #printMotif(profile, motifLength)
        
        convergence = True
        for i in range(numSeqs):
            old_instance[i] = instance[i]   # re-assign old instances as the current instances to determine convergence
            instance[i] = getBestMatchInSequence(sequences[i], profile, motifLength)    # find best match in each sequence (step 2 of EM) 
            if(old_instance[i] != instance[i]):
                convergence = False
        
    return profile

##################################################################
# Runs findMotif "repetition" number of times and returns the motif
# with the highest information content
# When calling this function, be sure to comment out
# printing in the findMotif function
# amino - boolean for whether giving AA or nucleotide sequences
# seqPrinting - whether or not full sequence is printed
##################################################################
def findBestMotif(sequences, motifLength, repetition, amino, seqPrinting):
    bestInfoContent = -10000000
    
    bestProfile = []
    instance = []
    
    totalTime = 0
    last10Time = 0
    for i in range(repetition):
        time1 = time.clock()
        m = findMotif(sequences, motifLength, amino)
        infoContent = calcInfoContent(m)
        #print(infoContent)

        # keep best motif found so far
        if infoContent > bestInfoContent:
                bestProfile = m
                bestInfoContent = infoContent
        time2 = time.clock()
        totalTime += (time2-time1)
        last10Time += (time2-time1)
        if(i % 10 == 0 and True):
            print(str(i) + "\t" + str(last10Time))
            last10Time = 0
	
    # display information about best motif
    m = bestProfile
    
    # print best motif
    if(seqPrinting):
        printMotif(m, motifLength)

    # find best match in each sequence   
    for i in range(len(sequences)):
        instance.append(getBestMatchInSequence(sequences[i], m, motifLength))
        if(seqPrinting):
            print(sequences[i] + "\t" + instance[i])

    # print runtime
    print("Total Execution Time: " + "{0:.2f}".format(totalTime) + " seconds")
    
    # print info content
    print ("\nInformation content: " + str(bestInfoContent)+ "\n")
        
    for i in range(len(instance)):
        print(str(i) + "\t" + instance[i] + "\t" + base2amino(instance[i]))
        
    return instance

##################################################################
# amino2bases - converts amino acids to nucleotide bases using 
# only one of the possible nucleotide sequences, not accounting
# for amino acids made by multiple base sequences
##################################################################
def amino2bases(aaSeq):
    nucSeq = ""
    for i in range(len(aaSeq)):
        if aaSeq[i] == 'F':
            nucSeq += 'UUU'
        elif aaSeq[i] == 'L':
            nucSeq += 'UUA'
        elif aaSeq[i] == 'I':
            nucSeq += 'AUU'
        elif aaSeq[i] == 'M':
            nucSeq += 'AUG'
        elif aaSeq[i] == 'V':
            nucSeq += 'GUU'
        elif aaSeq[i] == 'S':
            nucSeq += 'UCU'
        elif aaSeq[i] == 'P':
            nucSeq += 'CCU'
        elif aaSeq[i] == 'T':
            nucSeq += 'ACU'
        elif aaSeq[i] == 'A':
            nucSeq += 'GCU'
        elif aaSeq[i] == 'Y':
            nucSeq += 'UAU'
        elif aaSeq[i] == 'H':
            nucSeq += 'CAU'
        elif aaSeq[i] == 'Q':
            nucSeq += 'CAA'
        elif aaSeq[i] == 'N':
            nucSeq += 'AAU'
        elif aaSeq[i] == 'K':
            nucSeq += 'AAA'
        elif aaSeq[i] == 'D':
            nucSeq += 'GAU'
        elif aaSeq[i] == 'E':
            nucSeq += 'GAA'
        elif aaSeq[i] == 'C':
            nucSeq += 'UGU'
        elif aaSeq[i] == 'W':
            nucSeq += 'UGG'
        elif aaSeq[i] == 'R':
            nucSeq += 'AGA'
        elif aaSeq[i] == 'G':
            nucSeq += 'GGU'
    return nucSeq
    
##################################################################
# bases2amino - complement of amino2bases to reverse the processing
##################################################################
def base2amino(nucSeq):
    aaSeq = ""
    length = len(nucSeq)
    i = 0
    while (i + 3 <= length):
        if nucSeq[i:i+3] == 'UUU':
            aaSeq += 'F'
        elif nucSeq[i:i+3] == 'UUA':
            aaSeq += 'L'
        elif nucSeq[i:i+3] == 'AUU':
            aaSeq += 'I'
        elif nucSeq[i:i+3] == 'AUG':
            aaSeq += 'M'
        elif nucSeq[i:i+3] == 'GUU':
            aaSeq += 'V'
        elif nucSeq[i:i+3] == 'UCU':
            aaSeq += 'S'
        elif nucSeq[i:i+3] == 'CCU':
            aaSeq += 'P'
        elif nucSeq[i:i+3] == 'ACU':
            aaSeq += 'T'
        elif nucSeq[i:i+3] == 'GCU':
            aaSeq += 'A'
        elif nucSeq[i:i+3] == 'UAU':
            aaSeq += 'Y'
        elif nucSeq[i:i+3] == 'CAU':
            aaSeq += 'H'
        elif nucSeq[i:i+3] == 'CAA':
            aaSeq += 'Q'
        elif nucSeq[i:i+3] == 'AAU':
            aaSeq += 'N'
        elif nucSeq[i:i+3] == 'AAA':
            aaSeq += 'K'
        elif nucSeq[i:i+3] == 'GAU':
            aaSeq += 'D'
        elif nucSeq[i:i+3] == 'GAA':
            aaSeq += 'E'
        elif nucSeq[i:i+3] == 'UGU':
            aaSeq += 'C'
        elif nucSeq[i:i+3] == 'UGG':
            aaSeq += 'W'
        elif nucSeq[i:i+3] == 'AGA':
            aaSeq += 'R'
        elif nucSeq[i:i+3] == 'GGU':
            aaSeq += 'G'
        i += 3
    return aaSeq

######################################
# Testing
######################################
def test1():
    # search for 4-mer (length 4) motif in the following 4 sequences
    seq1 = "GTATACGATGTCTAGTATCAGCGGCATTAG"
    seq2 = "TAGCTGTACGTAGCGGCTTTAGCTGCAT"
    seq3 = "GACAGTCAGCGTTAGCTATATGCT"
    seq4 = "GCAGCAGTTGAGCAGCGATGATTTATCG"
    seq5 = "TCAGCGATTTATTTATTTATTTTA"

    sequences = []
    sequences.append(seq1)
    sequences.append(seq2)
    sequences.append(seq3)
    sequences.append(seq4)
    sequences.append(seq5)

    motifLength = 6

    #printMotif(findMotif(sequences, motifLength), motifLength)
    findBestMotif(sequences, motifLength, 10000, False, True)	# be sure to comment out printing when running findBestMotif
    
#test1()

def test2():
    seqs = generateSeqsFromFasta("ceHaspins_CloseRelatives_ExpressionList.txt")
    for i in range(len(seqs)):
        print(seqs[i])
        
#test2()

#Test Original Code of converting to nucleotides before attempting motif finding
def testAnalyzeHaspins(filename, aminoMotifLength, numberRuns, amino):
    aaSeqs = generateSeqsFromFasta(filename)
    nucSeqs = []
    for i in range(len(aaSeqs)):
        nucSeqs.append(amino2bases(aaSeqs[i]))
    
    instances = findBestMotif(nucSeqs, aminoMotifLength*3, numberRuns, amino, False)

#testAnalyzeHaspins("ceHaspins_CloseRelatives_ExpressionList.txt", 10, 100, True)
    
def testConsensus():
    profile = []
    aminoAcids = ['F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G']
    for i in range(3):
        profile.append([])
        for j in range(len(aminoAcids)):
            if(j < 3):
                profile[i].append(0.25)
            else:
                profile[i].append(0.1)
            
    #profile[2][5] = 0.3
    #profile[1][5] = 0.3

#testConsensus()
        
######################################
# Run Functions Here
######################################
def analyzeHaspins(filename, motifLength, numberRuns):
    aaSeqs = generateSeqsFromFasta(filename)
    printFullSeqs = False
    instances = findBestMotifAmino(aaSeqs, motifLength, numberRuns, printFullSeqs)

analyzeHaspins("ceHaspins_CloseRelatives_ExpressionList.txt", 5, 10)
        

