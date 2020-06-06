from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import generic_alphabet

protein_seq = Seq("EVRNAK", IUPAC.protein)
dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
protein_seq.alphabet = generic_alphabet
dna_seq.alphabet = generic_alphabet
protein_seq + dna_seq
nuc_seq = Seq("GATCGATGC", generic_nucleotide)#通用核苷酸类。可以和任意核酸类相加
dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
nuc_seq+dna_seq
my_seq=Seq("AGTACACTGGT",IUPAC.unambiguous_dna)
my_seq.alphabet
print(my_seq)
my_seq.complement()
my_seq.reverse_complement()
###类似字符串的操作
for index,letter in enumerate(my_seq):
    print("%i %s"%(index,letter))

my_seq.count('A')
GC(my_seq)
my_seq[4:12]
my_seq[1::3]
str(my_seq)
dna_seq.upper()
dna_seq.lowwer()

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
coding_dna.transcribe()#T→U
coding_dna.reverse_complement().transcribe()#true_transcribe
coding_dna.translate(to_stop=True,cds=True)#RNA和DNA都可以直接翻译,table参数可以选择密码子表

from Bio.Data import CodonTable
standard_table=CodonTable.unambiguous_dna_by_id[1]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_id[2]

from Bio.Seq import MutableSeq
mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
##或者
mutable_seq = my_seq.tomutable()
mutable_seq[5] = "C"
mutable_seq.remove("T")
mutable_seq.reverse()

new_seq=mutable_seq.toseq()

from Bio.Seq import UnknownSeq
unk_dna=UnknownSeq(20,alphabet=IUPAC.ambiguous_dna)

from Bio.Seq import reverse_complement,transcribe,back_transcribe,translate
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
reverse_complement(my_string)
transcribe(my_string)
back_transcribe(my_string)
translate(my_string)

