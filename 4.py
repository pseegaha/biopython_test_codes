from Bio import SeqIO
from Bio import pairwise2  # This is the version from the pull request
structure='TCCTGCTTCAGCCTACAGACCTGGGACTGCCACAGCTCATCACTGTGCCTGCATCCATAATAACTTCTTCAGCATGTTTTGGGCTCAGGCCTCATGGCAGCTGGCCAATGCTTATAAACTACTCTCAATCGCTAGCCCTGTACGTGGCCATTTGCCAAGGGCAGGGTAAAGCAAAGTCCTGGCACGAGAGTAGTTTATAAGCATTGGCCAGCTGCCATGAGGCCAACCCTGCCAAACAAGGACAGGAGACTCCTGAGGCAGGCTCTTCTGTCTTGGGAGGATGGTTCCAGGCCACTGATATTAAGGGTTAGGAGTTCAGTTCTCTGTGAGCTTAAAGGCTGATTATGGGG'
def complement(structure):     
    replacement1 = structure.replace('A', 'u')     
    replacement2 = replacement1.replace('U', 'a')     
    replacement3 = replacement2.replace('C', 'g')     
    replacement4 = replacement3.replace('G', 'c')     
    return replacement4.upper()  
print structure
RF4 = complement(structure)#[::-1][0:]
print RF4
RF5=reversed(RF4)
print RF5

for a in pairwise2.align.localms(structure, RF4, 1, -2, -2, 0):
        print pairwise2.format_alignment(*a)
        
def is_self_complementary(self):
        """
        Method to check whether a duplex contains self-complementary
        strands.  Always returns False for DNA-RNA duplexes
        """
        if str(self.structure) == str(self.RF4):
            return True
        else:
            return False
print is_self_complementary(structure)
    
def longest_substring(structure, RF4):
    t = [[0]*(1+len(RF4)) for i in range(1+len(structure))]
    l, xl = 0, 0
    for x in range(1,1+len(structure)):
        for y in range(1,1+len(RF4)):
            if structure[x-1] == RF4[y-1]:
                t[x][y] = t[x-1][y-1] 
                if t[x][y]>l:
                    l = t[x][y]
                    xl  = x
            else:
                t[x][y] = 0
    return structure[xl-l: xl]
print longest_substring(structure, RF4)