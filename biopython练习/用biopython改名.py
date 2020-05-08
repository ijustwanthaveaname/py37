from Bio import SeqIO
align=SeqIO.parse("C:\\Users\\62792\\Downloads\\.best.fas",'fasta')
seq=[]
idbox=[]#预防重名物种
for record in align:
    record.id=record.id.split('[',1)[1][0:-1]
    record.description = ''
    if record.id not in idbox:
        seq.append(record)
        idbox.append(record.id)
    else:
        record.id=record.id+'2'
        seq.append(record)

SeqIO.write(seq,"C:\\Users\\62792\\Downloads\\rename.fas",'fasta')