## Targeted DNA string to be translated to it protein confirmation
from __future__ import division 
input_file = open("me.fasta", "rU")
from Bio import SeqIO
for cur_record in SeqIO.parse(input_file, "fasta") :
    dnasequence=(cur_record.seq)  
    reverse= dnasequence.reverse_complement()
    print dnasequence
    print reverse

allowed_bases = ["A", "T", "G", "C","U","R","Y","S","W","K","M","B","D","H","V","N",".","-"]
total_dna_bases = 0
for base in allowed_bases:
    total_dna_bases = total_dna_bases + dnasequence.count(base)
    dna_fraction = total_dna_bases / len(dnasequence)
    print(total_dna_bases)
    print dna_fraction
    
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq, translate
from re import finditer    
# Get the starting position for each ORF in the dna sequence and translate.
ORFs = [translate(dnasequence[x.start():], table = 1, stop_symbol = '', to_stop= True) for x in finditer('ATG', str(dnasequence))]
# Get the starting position for each ORF in the reverse complement sequence and translate.
Os = [translate(reverse[x.start():], table = 1, stop_symbol = '', to_stop= True) for x in finditer('ATG', str(reverse))]

# Find the longest ORF.
longest_orf = max(map(str, ORFs), key=len)
longest_rev = max(map(str, Os), key=len)
# Print and save the answer.
print longest_orf
print longest_rev


#dnasequence = "ATGGCGACTGTCGAACCGGAAACCACCCCTACTCCTAATCCCCCGACTACAGAAGAGGAGAAAACGGAATCTAATCAGGAGGTTGCTAACCCAGAACACTATATTAAACATCCCCTACAGAACAGATGGGCACTCTGGTTTTTTAAAAATGATAAAAGCAAAACTTGGCAAGCAAACCTGCGGCTGATCTCCAAGTTTGATACTGTTGAAGACTTTTGGGCTCTGTACAACCATATCCAGTTGTCTAGTAATTTAATGCCTGGCTGTGACTACTC ACTTTTTAAGGATGGTATTGAGCCTATGTGGGAAGATGAGAAAAACAAACGGGGAGGACGATGGCTAATTACATTGAACAAACAGCAGAGACGAAGTGACCTCGATCGCTTTTGGCTAGAGACACTTCTGTGCCTTATTGGAGAATCTTTTGATGACTACAGTGATGATGTATGTGGCGCTGTTGTTAATGTTAGAGCTAAAGGTGATAAGATAGCAATATGGACTACTGAATGTGAAAACAGAGAAGCTGTTACACATATAGGGAGGGTATACAAGGAAAGGTTAGGACTTCCTCCAAAGATAGTGATTGGTTATCAGTCCCACGCAGACACAGCTACTAAGAGCGGCTCCACCACTAAAAATAGGTTTGTTGTTTAA"
## dictonary table with keys as triplet codons and the values as the proteins of that particular genetic triplet code
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
outputprotein = []
end = len(dnasequence) - (len(dnasequence) %3) - 1
for i in range(0,end,3):
        codons = dnasequence[i:i+3]
        if codons in genetable:
            aminoacids = genetable[codons]
            outputprotein.append(aminoacids)
        else:
            outputprotein.append("N")
#print (dnasequence)
#print " The translated protein sequence is listed as below for the above entered dna string"
print outputprotein
#print " The translated protein sequence is sourced from the data in the genetable which reads as below"
#print (genetable)