from Bio import SeqIO
sequences = []
handle1=open("homosapiens.gb","r")
#output1=open("homosapiens.txt","w")
handle=SeqIO.parse(handle1,"genbank")
for record in handle:
   for feature in record.features:
     if feature.type=="CDS":
          #print(">" + feature.qualifiers["protein_id"][0]+"\n")
          first=feature.qualifiers["translation"][0]+ "\n"
          for line in first:
            sequences.append(first)
new_list = list(set(sequences))
print new_list
with open('new_file.txt', 'w') as first:
    for seq in new_list:
        print seq
        