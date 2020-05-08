import os
inputpath=input("请输入序列文件所在的文件夹:")
genenames=input("请输入需要筛选的geneID文件地址")
outputpath=input("请输入输出文件夹及名称:")
def get_sequence(inputpath,genenames,outputpath):
    inputfasta=os.listdir(inputpath)    #当前文件夹所有的fasta文件名
    outputfile=open(outputpath,"w")     #输出的筛选的fasta文件
    genenamestable=[]
    seq=""
    with open(genenames,"r") as fk:
        for i in fk:
            genenamestable.append(i.strip())
    genenamestable1=[">"+i for i in genenamestable]
    for i in inputfasta:
        inputfile=open(os.path.join(inputpath,i),"r")
        for line in inputfile:
            if line[0]==">" and seq =="":
                geneID=line
            elif line[0]!=">":
                seq=seq+line.strip()+"\n"
            elif line[0]==">"and seq!="":
                if geneID.strip() in genenamestable1:
                    outputfile.write(geneID+seq)
                geneID=line
                seq=""
    if geneID in genenamestable1:
        outputfile.write(geneID+seq)
    outputfile.close()

get_sequence(inputpath,genenames,outputpath)

