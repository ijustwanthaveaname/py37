from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
a=SeqIO.parse("C:\\Users\\62792\\Desktop\\GCF_900518735.1_EBS10Xv2-PRI_genomic.gbff", "genbank")
wrt=[]
for rec in a:  #迭代器返回每个SeqRecord对象
    if rec.features:#调用SeqRecord的features对象列表
        for feature in rec.features:
            if feature.type == "CDS":
                wrt.append(feature.extract(rec))  #feature.location.extract(rec).seq或者feature.extract(rec,seq)都可以
SeqIO.write(wrt,"C:\\Users\\62792\\Desktop\\my.fasta","fasta")