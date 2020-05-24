from Bio import SeqIO
import sys
import os
from Bio.SeqRecord import SeqRecord
file=sys.argv[1]
a=SeqIO.parse(file, "genbank")
wrt=[]
for rec in a:  #迭代器返回每个SeqRecord对象
    if rec.features:#调用SeqRecord的features对象列表
        for feature in rec.features:
            if feature.type == "CDS":
                if "protein_id" in feature.qualifiers:
                    wrt.append('>'+feature.qualifiers["protein_id"][0]+'\n'+str(feature.location.extract(rec.seq)+'\n'))
                else:
                    wrt.append('>' + feature.qualifiers["gene"][0] + '\n' + str(
                        feature.location.extract(rec.seq) + '\n'))
#feature.location.extract(rec).seq或者feature.extract(rec,seq)都可以
# SeqIO.write(wrt,''.join(file.split(sep='.')[0:-1])+'_cds.fasta',"fasta")
with open(''.join(file.split(sep='.')[0:-1])+'_cds.fasta','w') as f:
    f.write(''.join(wrt))
