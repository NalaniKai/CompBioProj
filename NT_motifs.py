#IC and motifs
#!/usr/bin/python

#Luke Douville's part
#Look for motifs in the N terminus using algorithms that generate motifs.

#Locate the N terminus - it's at the beginning of the amino acid sequence
#Run algorithm to find motifs
#import sys
#
#sys.setrecursionlimit(50000)
aaFile = open("amino_acids.txt", "r")
aaText = aaFile.read()
#global d
d = {}
#global aaString
aaString = ""#For some reason thinks aaString is not initialized
#j = 0
def isolateAminoAcidSequences():#x, y):#aaString):
    global aaString
    aaString = ""
    for i in range(0, len(aaText)):
#        if j > 0:
#            j += 1
        if i + 1 < aaText.find(">", i, len(aaText)):
            aaString += aaText[i]#capture one amino acid sequence
#            for j in range(i, aaText.find(">", i, len(aaText))):
#                d[i] += aaText[j]
#        elif i == aaText.find(">"):
#            break
        elif aaText[i] == "\n" and aaText[i + 1] == ">" and i < len(aaText):
            global d
            d[i] = aaString
            aaString = ""
        else:#quit capturing amino acid sequence when the end of the sequence is reached
#            break
#            d[i] = aaString
#            key = i
#            global d
#            d[i] = aaString
#            d[key] = aaString
#            d[i] = aaString
#            break
            aaString = ""
#    x = i
#    y =
#    isolateAminoAcidSequences(x, y)
#    print("d = ", d)
    for i in range(0, len(aaText)):
        if i in d:
            print(i, " ", end = ",")
    print("\n\n", d[5324])
    print("\n\n\n")

#print(aaString)
isolateAminoAcidSequences()#0, len(aaText))#aaString)#Each amino acid sequence is now in d
#print(d)
#print("\n\n\n", d.keys())
#print(d.items())
#c = aaText.count("\n>")
#for i in range(0, len(aaText)):
#    if
print(d)
