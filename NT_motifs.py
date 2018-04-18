#IC and motifs
#!/usr/bin/python

#Luke Douville's part
#Look for motifs in the N terminus using algorithms that generate motifs.

#Locate the N terminus - it's at the beginning of the amino acid sequence
#Run algorithm to find motifs
aaFile = open("amino_acids.txt", "r")
aaText = aaFile.read()
aaString = ""
for i in range(0, len(aaText)):
    if i > aaText.find("\n"):
        aaString += aaText[i]
