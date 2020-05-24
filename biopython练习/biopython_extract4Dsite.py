####首先需要orthofinder构建物种单拷贝直系同源集，然后利用这个
from Bio import SeqIO
import os
import sys
filepath=sys.argv[1]
codon_table={
'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
file=SeqIO.parse(filepath,'fasta')


wrt = []
for rec in file:
    i=0
    wrt.append('>'+str(rec.id)+'\n')
    while i < len(rec.seq):
        if str(rec.seq[i:i+3]) in codon_table:
            wrt.append(str(rec.seq[i+3]))
        i+=3
    wrt.append('\n')
with open(filepath+'.4Dsite','w') as f:
    f.write(''.join(wrt))


