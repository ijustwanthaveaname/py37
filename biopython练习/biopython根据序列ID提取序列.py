from Bio import SeqIO
IDlist=open(r'D:\test\functionfitered_id.txt','r')
IDlist=[i.strip() for i in IDlist]
file=SeqIO.parse(r"D:\test\protein.faa",'fasta')
filter_protein=[]
for record in file:
    if record.id in IDlist:
        filter_protein.append(record)
SeqIO.write(filter_protein,r"D:\test\functionfitered_id.faa","fasta")