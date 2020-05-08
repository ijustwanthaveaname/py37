from Bio import SeqIO
file=SeqIO.parse(r"C:\Users\62792\Downloads\rename.phy",'phylip')
SeqIO.write(sequences=file,handle='C:\\Users\\62792\\Downloads\\rename1.fas',format='fasta')

# from Bio import SeqIO
#
# records = SeqIO.parse("THIS_IS_YOUR_INPUT_FILE.fasta", "fasta")
# count = SeqIO.write(records, "THIS_IS_YOUR_OUTPUT_FILE.stockholm", "stockholm")
# print("Converted %i records" % count)