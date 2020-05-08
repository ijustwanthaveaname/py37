inputpath=input("请输入批量重命名单拷贝基因集序列文件夹:")
outputpath=input("请输入重命名文件输出文件夹:")
def renameid(inputpath,outputpath):
    import os
    taxa=1
    seq=""
    allfasta=os.listdir(inputpath)
    for i in allfasta:
        inputfile= open(os.path.join(inputpath, i), "r")
        outputfile=open(outputpath+"\\"+i,'w')
        for line in inputfile:
            if line[0]==">" and seq=="":
                pass
            elif line[0]!=">":
                seq=seq+line
            elif line[0]==">" and seq!="":
                outputfile.write(">species"+str(taxa)+"\n"+seq)
                seq=""
                taxa=taxa+1
        outputfile.write(">species" + str(taxa) + "\n" + seq)
        seq = ""
        outputfile.close()
        taxa=1
renameid(inputpath,outputpath)

