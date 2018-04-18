#IC and motifs
#!/usr/bin/python

#Luke Douville's part
#Look for motifs in the N terminus using algorithms that generate motifs.

#Locate the N terminus - it's at the beginning of the amino acid sequence
#Run algorithm to find motifs
aaFile = open("amino_acids.txt", "r")
aaText = aaFile.read()
d = {}
aaString = ""#For some reason thinks aaString is not initialized
#j = 0
def isolateAminoAcidSequences():
#    aaString = ""
    for i in range(0, len(aaText)):
#        if j > 0:
#            j += 1
        if i + 1 < aaText.find(">", i, len(aaText)):
            aaString += aaText[i]
#            for j in range(i, aaText.find(">", i, len(aaText))):
#                d[i] += aaText[j]
#        elif i == aaText.find(">"):
#            break
#        else:
#            break
#            d[i] = aaString
#            key = i
#            d[key] = aaString
            d[i] = aaString
    print("d = ", d)

#print(aaString)
#print(d)
isolateAminoAcidSequences()
