import pandas as pd
import numpy as np
import os
inputpath=input("请输入序列文件所在的文件夹:")
orthogroupspath=input("请输入Orthogroups.txt文件路径:")
singlecopypath=input("请输入SingleCopyOrthogroups.txt文件路径:")
outputpath=input("请输入输出文件夹:")
orthogroups=pd.read_csv(orthogroupspath,sep=" ",header=None)
with open(singlecopypath) as fuck:
    singlecopy=[i.strip() for i in fuck]
orthogroups=orthogroups.set_index([0])
scopy=orthogroups.loc[singlecopy]
scopy=scopy.dropna(axis=1,how="all")
scopylist=scopy.values.tolist()
outputname=0
seq=""
inputfasta = os.listdir(inputpath)  # 当前文件夹所有的fasta文件名
for groups in scopylist:
    outputname=outputname+1
    outputpathx=outputpath+"\\"+str(outputname)+".fasta"
    outputfile=open(outputpathx,"w")     #输出的筛选的fasta文件
    for fastaname in inputfasta:
        inputfile=open(os.path.join(inputpath,fastaname),"r")
        for line in inputfile:
            if line[0]==">" and seq =="":
                geneID=line.split()[0].strip()
            elif line[0]!=">":
                seq=seq+line.strip()+"\n"
            elif line[0]==">"and seq !="":
                if geneID[1:] in groups:
                    outputfile.write(geneID+"\n"+seq)
                geneID=line.split()[0].strip()
                seq=""
    if geneID[1:] in groups:
        outputfile.write(geneID+"\n"+seq)
    outputfile.close()


