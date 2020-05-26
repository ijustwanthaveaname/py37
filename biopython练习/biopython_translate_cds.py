from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os
def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_record.seq.translate(cds=True), \
                     id = nuc_record.id, \
                     description = "translation of CDS, using default table")
proteins = (make_protein_record(nuc_rec) for nuc_rec in \
            SeqIO.parse(sys.argv[1], "fasta"))
SeqIO.write(proteins, os.path.dirname(sys.argv[1])+'trans_'+os.path.basename(sys.argv[1]), "fasta")
