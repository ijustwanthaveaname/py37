from Bio import AlignIO
a=AlignIO.read("C:\\Users\\62792\\Desktop\\coi.fas","fasta") #用numpy储存：a=np.array([list(rec) for rec in alignment], np.character) 
                                   #若以列形式储存：align_array = np.array([list(rec) for rec in alignment], np.character, order="F")
b=a[:,:-1]
print(a+b)