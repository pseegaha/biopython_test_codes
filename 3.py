#pruthvi swamy
from Bio import SeqIO
nucleotides=['A','C','G','T']      # Nucleotide seq bases which are focussed in this program
motif=[]
infile="sequence.gb"
#infile command for the input of the gen bank sequence
List=[a+b+c+d+f for a in nucleotides for b in nucleotides for c in nucleotides for d in nucleotides for f in nucleotides] # attribuing diff variables for diff nucleotides
utr=[]
#In molecular genetics, the three prime untranslated region (3'-UTR) is the section of messenger RNA (mRNA) #
#that immediately follows the translation termination codon --info source: wikipedia. 
records=SeqIO.parse(infile,'genbank') # parsing command for the acess of the genbank format records
for record in records:   # working on the record of input records
        for feature in record.features: #accessing the features of the selected record
            if feature.type=='CDS':   #using the if syntax for the CDS feature type
                start=feature.location._end.position  # reading the input if the feature locations which start with the termination codon
            else:
                continue                            # otherwise moving to the next step of continuin the if 
        seq=str(record.seq[start:])  #storing the found utr in a variable
        #print len(seq)
        utr.append(seq)                 # using dictionary, appending the variable with list to get output in dictionaries
        #print UTR

for sequence in utr:                    #syntax to find the motifs in the utr
    List1=[0 for m in List]
    motif=dict(zip(List,List1))  
    for m in range(len(sequence)-5):
        for k in motif.keys():
            if sequence[m:m+5]==k:
                motif[k]+=1

print sum(motif.values())              #total motif values
print motif   # The motifs found in the input

    