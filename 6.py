## Targeted DNA string to be translated to it protein confirmation
#print "Using the codon dictionaries, translate  the following nucleotide string:'ATGGCGACTGTCGAACCGGAAACCACCCCTACTCCTAATCCCCCGACTACAGAAGAGGAGAAAACGGAATCTAATCAGGAGGTTGCTAACCCAGAACACTATATTAAACATCCCCTACAGAACAGATGGGCACTCTGGTTTTTTAAAAATGATAAAAGCAAAACTTGGCAAGCAAACCTGCGGCTGATCTCCAAGTTTGATACTGTTGAAGACTTTTGGGCTCTGTACAACCATATCCAGTTGTCTAGTAATTTAATGCCTGGCTGTGACTACTC ACTTTTTAAGGATGGTATTGAGCCTATGTGGGAAGATGAGAAAAACAAACGGGGAGGACGATGGCTAATTACATTGAACAAACAGCAGAGACGAAGTGACCTCGATCGCTTTTGGCTAGAGACACTTCTGTGCCTTATTGGAGAATCTTTTGATGACTACAGTGATGATGTATGTGGCGCTGTTGTTAATGTTAGAGCTAAAGGTGATAAGATAGCAATATGGACTACTGAATGTGAAAACAGAGAAGCTGTTACACATATAGGGAGGGTATACAAGGAAAGGTTAGGACTTCCTCCAAAGATAGTGATTGGTTATCAGTCCCACGCAGACACAGCTACTAAGAGCGGCTCCACCACTAAAAATAGGTTTGTTGTTTAA'into a protein sequence using single-letter amino acid alphabet."
dnasequence = "ATGGCGACTGTCGAACCGGAAACCACCCCTACTCCTAATCCCCCGACTACAGAAGAGGAGAAAACGGAATCTAATCAGGAGGTTGCTAACCCAGAACACTATATTAAACATCCCCTACAGAACAGATGGGCACTCTGGTTTTTTAAAAATGATAAAAGCAAAACTTGGCAAGCAAACCTGCGGCTGATCTCCAAGTTTGATACTGTTGAAGACTTTTGGGCTCTGTACAACCATATCCAGTTGTCTAGTAATTTAATGCCTGGCTGTGACTACTC ACTTTTTAAGGATGGTATTGAGCCTATGTGGGAAGATGAGAAAAACAAACGGGGAGGACGATGGCTAATTACATTGAACAAACAGCAGAGACGAAGTGACCTCGATCGCTTTTGGCTAGAGACACTTCTGTGCCTTATTGGAGAATCTTTTGATGACTACAGTGATGATGTATGTGGCGCTGTTGTTAATGTTAGAGCTAAAGGTGATAAGATAGCAATATGGACTACTGAATGTGAAAACAGAGAAGCTGTTACACATATAGGGAGGGTATACAAGGAAAGGTTAGGACTTCCTCCAAAGATAGTGATTGGTTATCAGTCCCACGCAGACACAGCTACTAAGAGCGGCTCCACCACTAAAAATAGGTTTGTTGTTTAA"
## dictonary table with keys as triplet codons and the values as the proteins of that particular genetic triplet code
#the data in the table is extracted from the online sources wikipedia
genetable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
outputprotein = []  #assinging the blank list so that it can be latter appended to get the results in the lists format
end_dna = len(dnasequence) - (len(dnasequence) %3) - 1 #getting the dna length right
for i in range(0,end_dna,3):                           #for loop to get the dna read in triplets
        codons = dnasequence[i:i+3]
        if codons in genetable:                         #command line to get inside with the gene table dictionary
            aminoacids = genetable[codons]
            outputprotein.append(aminoacids)             #appending the values of the genetable with the output list to generate output if the condition satisfies
        else:
            outputprotein.append("N")                    #this condition satisfies incase the above one doesn't cling
print (dnasequence)
print " The translated protein sequence is listed as below for the above entered dna string"
print outputprotein
print " The translated protein sequence is sourced from the data in the genetable which reads as below"
print (genetable)