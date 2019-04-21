#pruthvi swamy
from Bio import SeqIO
sequences = []
handle1=open("homosapiens.gb","r")
#output1=open("homosapiens.txt","w")
handle=SeqIO.parse(handle1,"genbank")
for record in handle:
   for feature in record.features:
     if feature.type=="CDS":
          print(">" + feature.qualifiers["protein_id"][0]+"\n")
          first=(feature.qualifiers["translation"][0]+ "\n")
          print first
          for line in first:
           sequences.append(first)
new = list(set(sequences))
print new
with open('newfile.fasta', 'w') as first:
    for seq in new:
        print seq
        first.write(seq)