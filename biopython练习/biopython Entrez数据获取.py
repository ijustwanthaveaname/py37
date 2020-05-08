from Bio import Entrez
Entrez.email='A.N.Other@examole.com'
handle=Entrez.einfo()#查看可用的数据库
result=Entrez.read(handle)
handle = Entrez.einfo(db="pubmed")
record = Entrez.read(handle)
record["DbInfo"]["Description"]
