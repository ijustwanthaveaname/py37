#!/usr/bin/env python3
import gzip
from Bio import SeqIO
long_seq = []
for record in SeqIO.parse("D:\\数据\\Gekko japonicus多疣壁虎\\GCF_001447785.1_Gekko_japonicus_V1.1_genomic.fna","fasta"):
    length=len(record.seq)
    if not int(0.4*length)*'N' in record.seq:
        long_seq.append(record)
SeqIO.write(long_seq, "D:\\数据\\Gekko japonicus多疣壁虎\\filltergecko2.fna", "fasta")