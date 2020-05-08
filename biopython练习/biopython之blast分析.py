from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
#result_handle=NCBIWWW.qblast('blastn','nt','8332116')
result_handle =r'D:\test\tblastn.xml'#需要解析的批量blast的xml
blast_records = open(result_handle)
E_VALUE_THRESH=1e-5
i=0
for blast_record in NCBIXML.parse(blast_records):
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            for seq_record in SeqIO.parse('D:\\test\\genomic.fna','fasta'):
                if  alignment.hit_id.split('|')[1].strip() == seq_record.id.strip():
                    if hsp.sbjct_start-1000 < 0:
                        start=0
                    else:
                        start=hsp.sbjct_start-1000
                    if hsp.sbjct_end + 1000 > len(seq_record.seq):
                        end=len(seq_record.seq)
                    else:
                        end=hsp.sbjct_end+1000
                    i+=1
                    subject_protein = open('D:\\protein{}.faa'.format(i), 'w')
                    subject_protein.write('>'+blast_record.query+'\n'+str(hsp.sbjct)+'\n')#以queryid命名
                    subject_protein.close()
                    genome_region = open('D:\\genome_region{}.fna'.format(i),'w')
                    genome_region.write('>'+seq_record.id+'\t'+str(hsp.sbjct_start)+'\t'+str(hsp.sbjct_end)+'\n'+str(seq_record.seq[start:end])+'\n')
                    genome_region.close()
                    # print(seq_record.seq[start:end])



            # if hsp.expect < E_VALUE_THRESH:
            #     print ('****Alignment****')
            #     print ('sequence:', alignment.title)
            #     print ('hitid:',alignment.hit_id.split('|')[1].strip())
            #     print ('length:', alignment.length)
            #     print ('e value:', hsp.expect)
            #     print('query_LOC',hsp.query_start,hsp.query_end)
            #     print (hsp.query[0:75] + '...')
            #     print (hsp.match[0:75] + '...')
            #     print('subject_LOC',hsp.sbjct_start,hsp.sbjct_end)
            #     print (hsp.sbjct[0:75] + '...')

