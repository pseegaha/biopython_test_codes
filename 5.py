from Bio import SeqIO
bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}    #base pairs
r='TCCTGCTTCAGCCTACAGACCTGGGACTGCCACAGCTCATCACTGTGCCTGCATCCATAATAACTTCTTCAGCATGTTTTGGGCTCAGGCCTCATGGCAGCTGGCCAATGCTTATAAACTACTCTCAATCGCTAGCCCTGTACGTGGCCATTTGCCAAGGGCAGGGTAAAGCAAAGTCCTGGCACGAGAGTAGTTTATAAGCATTGGCCAGCTGCCATGAGGCCAACCCTGCCAAACAAGGACAGGAGACTCCTGAGGCAGGCTCTTCTGTCTTGGGAGGATGGTTCCAGGCCACTGATATTAAGGGTTAGGAGTTCAGTTCTCTGTGAGCTTAAAGGCTGATTATGGGG'
print("The input structure is", r) 
#printing the input structure 


def reversecom(r):              #defining the reverse complement for the input structure
    return ''.join(bases[p] for p in r[::-1])
print ("The reversed complement of the input structure is",reversecom(r)) #printing the reverse complement output


def longstem(r):             #defining for the longest stem 
    m = len(r)    #length of the structure
    #print m
    n = int(m/2) #length of the possible longest possible stem
    #print n
    candidate = ''
    k = 1
    #print k

    while k <= n and len(candidate) == k - 1:
        for j in xrange(m-2*k+1):
            t = r[j:k+j]
            #print t
            if reversecom(t) in r[k+j:]:
                candidate = t
                #print candidate
                break
        k +=1
        #print k
    return candidate
    
print ('The potential longest stem loop sequence is',longstem(r))