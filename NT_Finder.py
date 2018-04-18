
########################################################
# convertFileToSequence -- converts FASTA file to string
# from lab 1
# Exercise 0
# convertFileToSequence - takes a FASTA file and returns the
# DNA sequence as a string
#######################################################
def convertFileToSequence(filename):
    # read in file
    file = open(filename)
    
    # read in first line
    header = file.readline()
    if (header[0] == '>'):
        print("in FASTA format")
    else:
        print("invalid format")
        file.close()
        # error code of -1 (so functions that call this can check)
        return -1
    
    # read in rest of file
    sequence = file.read()
    
    # close file
    file.close()

    # remove all return and newline characters
    sequence = sequence.replace("\r", "")
    sequence = sequence.replace("\n", "")
    return sequence


########################################################
# translate -- translate each DNA codon to the corresponding
# single letter amino acid sequence
#######################################################
def translate(inputfile, outputfile):
    # convert the input file to a string
    string = convertFileToSequence(inputfile)
    
    # open the file to write to
    out = open(outputfile, "w")
    out.write("> " + inputfile + " translated to amino acids")
    
    # counter to keep track of how many characters were written to the file
    count = 0
    
    # identifies the codon and finds the corresponding single letter amino acid
    for i in range(0, len(string), 3):
        if count % 60 == 0:
            out.write("\n")
        codon = string[i:i+3]
        if(codon == "TTT" or codon =="TTC" ):
            out.write("F")
            count += 1
        elif(codon == "TTA" or codon =="TTG" or codon =="CTT" or codon =="CTC" or codon == "CTA" or codon =="CTG"):
            out.write( "L")
            count += 1
        elif(codon == "ATT" or codon =="ATC" or codon =="ATA"):
            out.write( "I")
            count += 1
        elif(codon == "GTT" or codon =="GTC" or codon =="GTA"or codon =="GTG"):
            out.write( "V")
            count += 1
        elif(codon == "TTA" or codon =="TTG" or codon =="CTT"or codon =="CTC"or codon =="CTA" or codon =="CTG"):
            out.write( "L")
            count += 1
        elif(codon == "ATG" ):
            out.write("M")
            count += 1
        elif(codon == "TCT" or codon =="TCC" or codon =="TCA" or codon =="TCG"):
            out.write("S")
            count += 1
        elif(codon == "CCT" or codon =="CCC" or codon =="CCA" or codon =="CCG"):
            out.write("P")
            count += 1
        elif(codon == "ACT" or codon =="ACC" or codon =="ACA" or codon =="ACG"):
            out.write( "T" )
            count += 1
        elif(codon == "GCT" or codon =="GCC" or codon =="GCA" or codon =="GCG"):
            out.write( "A")
            count += 1
        elif(codon == "TAT" or codon =="TAC"):
            out.write( "Y")
            count += 1
        elif(codon == "TAA" or codon =="TAG" or codon =="TGA"):
            out.write( "X")
            count += 1
        elif(codon == "CAT" or codon =="CAC"):
            out.write( "H")
            count += 1
        elif(codon == "CAA" or codon =="CAG"):
            out.write( "Q")
            count += 1
        elif(codon == "AAT" or codon =="AAC" ):
            out.write( "N")
            count += 1
        elif(codon == "AAA" or codon =="AAG" ):
            out.write( "K")
            count += 1
        elif(codon == "GAT" or codon =="GAC"):
            out.write("D")
            count += 1
        elif(codon == "GAA" or codon =="GAG"):
            out.write("E")
            count += 1
        elif(codon == "TGT" or codon =="TGC"):
            out.write("C")
            count += 1
        elif(codon == "TGG"):
            out.write( "W" )
            count += 1
        elif(codon == "CGT" or codon =="CGC" or codon =="CGA" or codon =="CGG"):
            out.write( "R")
            count += 1
        elif(codon == "AGT" or codon =="AGC"):
            out.write( "S")
            count += 1
        elif(codon == "AGA" or codon =="AGG" ):
            out.write( "R")
            count += 1
        elif(codon == "GGT" or codon =="GGC" or codon =="GGA" or codon =="GGG"):
            out.write("G")
            count += 1

    # close the file
    out.close()


######################################################################
# localAlignmentScore
# @param s1 the first DNA sequence as a string
# @param s2 the second DNA sequence as a string
# @return the maximum score in the cost table
#
# Determine the score of the optimal local alignment of two strings
######################################################################
def localAlignmentScore(s1, s2):
    
    # Scoring system
    MATCH = 5
    MISMATCH = -4
    GAP = -6
                
    # set table size
    NUM_ROWS = len(s2)+1
    NUM_COLS = len(s1)+1
                
    # Create table and fill it with zeros
    costs = createTable(NUM_ROWS, NUM_COLS, 0)
                
    # Create table for getting back the optimal alignment, fill table with "A"
    # Suggest you use "D", "L", and "T" for diagonal, left, and top
    directions = createTable(NUM_ROWS, NUM_COLS, "A")
                
    # set the first position in the directions table to "F"
    directions[0][0] = "F"
                
    maxValue = 0
    maxRowPos = 0
    maxColPos = 0
                            
    for i in range(1, NUM_ROWS):
        costs[i][0] = 0
        directions[i][0] = "F"
                                    
    for i in range(1, NUM_COLS):
        costs[0][i] = 0
        directions[0][i] = "F"
                                            
    for i in range(1, NUM_ROWS):
        for j in range(1, NUM_COLS):
            if s2[i-1] == s1[j-1]:
                costs[i][j] = max(0, costs[i][j-1] + GAP, costs[i-1][j] + GAP, costs[i-1][j-1] + MATCH)
                if costs[i][j] == 0:
                    directions[i][j] = "F"
                elif costs[i][j] == costs[i][j-1] + GAP:
                    directions[i][j] = "L"
                elif costs[i][j] == costs[i-1][j] + GAP:
                    directions[i][j] = "T"
                else:
                    directions[i][j] = "D"
            else:
                costs[i][j] = max(0, costs[i][j-1] + GAP, costs[i-1][j] + GAP, costs[i-1][j-1] + MISMATCH)
                if costs[i][j] == 0:
                    directions[i][j] = "F"
                elif costs[i][j] == costs[i][j-1] + GAP:
                    directions[i][j] = "L"
                elif costs[i][j] == costs[i-1][j] + GAP:
                    directions[i][j] = "T"
                else:
                    directions[i][j] = "D"
            if costs[i][j] > maxValue:
                maxValue = costs[i][j]
                maxRowPos = i
                maxColPos = j

    # Print out table (only useful for small tables - used for debugging)
    printTable(costs, "costs.txt")
    printTable(directions, "directions.txt")
    
    # find optimal alignment
    align(directions, s1, s2, maxRowPos, maxColPos, "local_alignment.txt")
        
    # return optimal score (lower right-hand cell in table]
    return costs[maxRowPos][maxColPos]


####################################################################
# createTable
# Create a 2D table with the given number of rows and columns
# and fills all entries with value given as a parameter
# (function completed for you)
####################################################################
def createTable(numRows, numCols, value):
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


##################################################################
# printTable
# @param table the 2D table
#
# Print 2D table to file (only useful for small tables for short
# strings)
# Should have tabs between the values on each row
# Useful function for debugging purposes
# Students: Complete this function
##################################################################
def printTable(table, filename):
    file = open(filename, "w")
    for i in range(len(table)):
        for j in range(len(table[0])):
            file.write(str(table[i][j]) + "\t")
        file.write("\n\n")
    file.close()
    return

################################################################
# Reconstruct the optimal alignment and print the alignment
# to a file. Because the sequences can be long, print the
# alignment 50 characters on one line, the other string of 50 characters
# on the next line, and then skip one line, as follows:
# AATT--GGCTATGCT--C-G-TTACGCA-TTACT-AA-TCCGGTC-AGGC
# AAATATGG---TGCTGGCTGCTT---CAGTTA-TGAACTCC---CCAGGC
#
# TATGGGTGCTATGCTCG--T--TACG-CA
# TCAT--TGG---TGCTGGCTGCTT--ACA
#
# align
# @param direction the 2D direction table
# @param s1 the first DNA sequence
# @param s2 the second DNA sequence
# @param maxRow the row of the maximum score in the cost table
# @param maxCol the col of the maximum score in the cost table
# @filename the name of the output file
###############################################################
def align(direction, s1, s2, maxRow, maxCol, filename):
    
    path = direction[maxRow][maxCol]
    
    row = maxRow-1
    col = maxCol-1
    mRow = maxRow
    mCol = maxCol
    index = -1
                    
    optimalAlignments1 = ""
    optimalAlignments2 = ""
    count = 0
                                
    file = open(filename, "w")
        
    while path != "F":
        if path == "D":
            optimalAlignments1 = s1[col] + optimalAlignments1
            optimalAlignments2 = s2[row] + optimalAlignments2
            count += 1
            path = direction[mRow-1][mCol-1]
            if path == "F":
                index = col
            mRow -= 1
            mCol -= 1
            col -= 1
            row -= 1
        elif path == "L":
            optimalAlignments1 = s1[col] + optimalAlignments1
            optimalAlignments2 = "-" + optimalAlignments2
            count += 1
            path = direction[mRow][mCol-1]
            if path == "F":
                index = col
            mCol -= 1
            col -= 1
        else:
            optimalAlignments1 = "-" + optimalAlignments1
            optimalAlignments2 = s2[row] + optimalAlignments2
            count += 1
            path = direction[mRow-1][mCol]
            if path == "F":
                index = col
            mRow -= 1
            row -= 1
        # print the strings to the output file, 50 characters at a time
        if (path == "F") or (count % 50 == 0):
            file.write(optimalAlignments1 + "\n" + optimalAlignments2)
            file.write("\n\n")
            count = 0
            optimalAlignments1 = ""
            optimalAlignments2 = ""
                
    file.close()
    return


###################################################
# main
###################################################
translate("human_haspin_FASTA.txt", "human_haspin_amino_acid.txt")






