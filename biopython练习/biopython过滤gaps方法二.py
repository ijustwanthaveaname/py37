from Bio import SeqIO
records = (rec for rec in SeqIO.parse("D:\\数据\\Gekko japonicus多疣壁虎\\GCF_001447785.1_Gekko_japonicus_V1.1_genomic.fna", "fasta") if not int(0.4*len(rec.seq))*'N' in rec.seq)
SeqIO.write(records, "filtergekko3.fna", "fasta")
