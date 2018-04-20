#IC and motifs
#!/usr/bin/python

#Luke Douville's part
#Look for motifs in the N terminus using algorithms that generate motifs.

#Locate the N terminus - it's at the beginning of the amino acid sequence
#Run algorithm to find motifs

import motif_functions

aaFile = open("amino_acids.txt", "r")
aaText = aaFile.read()
d = {}
def isolateAminoAcidSequences():#x, y):#aaString):
    global aaString
    aaString = ""
    for i in range(0, len(aaText)):
        if i + 1 < aaText.find(">", i, len(aaText)):
            aaString += aaText[i]#capture one amino acid sequence
        elif aaText[i] == "\n" and aaText[i + 1] == ">" and i < len(aaText):
            global d
            d[i] = aaString
            aaString = ""
        else:#quit capturing amino acid sequence when the end of the sequence is reached
            aaString = ""

    for i in range(0, len(aaText)):
        if i in d:
            print(i, " ", end = ",")
    print("\n\n", d[5324])
    print("\n\n\n")

isolateAminoAcidSequences()#Each amino acid sequence is now in d
#print("d = ", d)

#Convert amino acid sequences to RNA
j = aaText.find(">", 0, len(aaText))
numAA = 0
n = 1
for i in range(len(aaText)):
    if aaText.find(">", i, len(aaText)) != -1:
        numAA += 1
for i in range(len(aaText)):
#    for j in range(0, len(aaText), j + aaText.find(">", j, len(aaText))):
#    while j < len(d):
    if i in d:
        print("\nDNA sequence from amino acid sequence {}: ".format(n), motif_functions.amino2bases(d[i]))
#            j += aaText.find(">", i, len(aaText))
        n += 1

