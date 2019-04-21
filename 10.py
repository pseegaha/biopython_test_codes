
        
sequences = []
with open('homosapiens.gb', 'r') as f:
    for line in f:
        for feature in line:
         if feature.type=="CDS":
            first=(feature.qualifiers["translation"])
            sequences.append(line)
new_list = list(set(sequences))
with open('file.txt', 'w') as f:
    for seq in new_list:
             print(seq)




