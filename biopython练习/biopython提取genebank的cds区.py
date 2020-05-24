from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
a=SeqIO.parse("C:\\Users\\62792\\Desktop\\GCA_003118565.1_Ppicta_assembly_v1_genomic.gbff", "genbank")
for rec in a:  #迭代器返回每个SeqRecord对象
    if rec.features:#调用SeqRecord的features对象列表
        for feature in rec.features:
            if feature.type == "CDS":
                # print (feature.location)
                print (feature.qualifiers["protein_id"])
                print (str(feature.location.extract(rec.seq)))   #feature.location.extract(rec).seq或者feature.extract(rec,seq)都可以
                