import os
def get_sequence(inputpath,gene,outputpath):
    inputfasta=os.listdir(inputpath)    #当前文件夹所有的fasta文件名
    outputfile=open(outputpath,"w")     #输出的筛选的fasta文件
    seq=""
    for i in inputfasta:
        inputfile=open(os.path.join(inputpath,i),"r")
        for line in inputfile:
            if line[0]==">" and seq =="":
                geneID=line
            elif line[0]!=">":
                seq=seq+line.strip()+"\n"
            elif line[0]==">"and seq!="":
                if geneID.strip() in gene:
                    outputfile.write(geneID+seq)
                geneID=line
                seq=""
    if geneID.strip() in gene:
        outputfile.write(geneID+seq)
    outputfile.close()


