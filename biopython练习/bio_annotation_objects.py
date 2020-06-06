import os
os.chdir('..')
from Bio.Seq import Seq
simple_seq = Seq("GATC")
from Bio.SeqRecord import SeqRecord
simple_seq_r = SeqRecord(simple_seq)
simple_seq_r.id = "AC12345"
simple_seq_r.description = "Made up sequence I wish I could write a paper about"
simple_seq_r.seq
simple_seq_r = SeqRecord(simple_seq, id="AC12345")
simple_seq_r.annotations["evidence"] = "None. I just made it up."
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]

from Bio import SeqIO
# record=SeqIO.read("NC_005816.fna", "fasta")
# record.seq
# record.id
# record.name
# record.description

record = SeqIO.read("NC_005816.gb", "genbank")
record.id#include the version suffix
record.name#from the LOCUS line
record.description#from the DEFINITION line
len(record.features)
record.features#special list which include SeqFeature objects
for feature in record.features:
    print (feature.type)
    print (feature.location)# The location of the SeqFeature on the sequence that you are dealing with
    print (feature.location.strand)#1 for top strand -1 for the bottom strand ,0  if the strand is important but is unknown,None if it doesn’t matter
    print (feature.qualifiers)
    if "protein_id" in feature.qualifiers:
        print (feature.qualifiers["protein_id"])