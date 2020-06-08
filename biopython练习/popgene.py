from Bio import AlignIO
from collections import Counter
from scipy.special import comb, perm
import sys
def _read_alignment(path):
    alignment = AlignIO.read(path, "fasta")
    return alignment
def _sum_pairwise_multiple(data):
    c=0
    for i in range(0,len(data)-1):
        c+=sum([data[i]*sb for sb in data[i+1:]])
    return c

def snp_count(alignmentpath):
    alignment=_read_alignment(alignmentpath)
    snpcount = 0
    for i in range(0,len(alignment[0])):
        if len(set(alignment[:,i]))>1:
            snpcount+=1
    return snpcount

def thetapi_calculate(alignmentpath):
    alignment=_read_alignment(alignmentpath)
    allcounts=[]
    for i in range(0,len(alignment[0])):
        cnt = Counter()
        if len(set(alignment[:,i]))>1:
            for base in list(alignment[:,i]):
                cnt[base]+=1
            mutiply=_sum_pairwise_multiple(list(dict(cnt).values()))
            allcounts.append(mutiply)
        else:
            pass
    combine=comb(len(alignment[:,0]), 2)
    return sum(allcounts)/combine

def thetaw_calculate(alignmentpath):
    alignment=_read_alignment(alignmentpath)
    snpcount=snp_count(alignmentpath)
    thetaw=snpcount/sum([1/i for i in range(1,len(alignment))])
    return thetaw

def calculate_tajimad(alignmentpath):
    alignment=_read_alignment(alignmentpath)
    S=snp_count(alignmentpath)
    n=len(alignment)
    a1=sum([1/i for i in range(1,n)])
    a2=sum([1/i**2 for i in range(1,n)])
    b1=(n+1)/(3*(n-1))
    b2=2*(n**2+n+3)/(9*n*(n-1))
    c1=b1-1/a1
    c2=b2-(n+2)/(a1*n)+a2/a1**2
    e1=c1/a1
    e2=c2/(a1**2+a2)
    thetapi=thetapi_calculate(alignmentpath)
    D=(thetapi-S/a1)/(e1*S+e2*S*(S-1))**(1/2)
    return D

def fst_calculate(allalignmentpath,*alignmentspath):
    allpi=thetapi_calculate(allalignmentpath)
    allcount=0
    for alnpath in alignmentspath:
        allcount+=thetapi_calculate(alnpath)
    return 1-allcount/len(alignmentspath)/allpi

def Nm_calculate(allalignmentpath,*alignmentspath):
    fst=fst_calculate(allalignmentpath,*alignmentspath)
    return 1/(4*fst)-1/4

if sys.argv[1]=="-nsnp":
    print(snp_count(sys.argv[2]))

elif sys.argv[1]=="-pi":
    print(thetapi_calculate(sys.argv[2]))

elif sys.argv[1]=="-w":
    print(thetaw_calculate(sys.argv[2]))

elif sys.argv[1]=="-tajimaD":
    print(calculate_tajimad(sys.argv[2]))

elif sys.argv[1]=="-fst":
    print(fst_calculate(sys.argv[2],*sys.argv[3:]))
    # print (sys.argv[2],sys.argv[3:])
elif sys.argv[1]=="-Nm":
    print(Nm_calculate(sys.argv[2],*sys.argv[3:]))

else:
    raise TypeError("You must input specified arguments!")
# if __name__ == '__main__':
#     allpath=r"C:\Users\62792\Desktop\test.faa"
#     alignmentpath1 = r"C:\Users\62792\Desktop\test1.faa"
#     alignmentpath2 = r"C:\Users\62792\Desktop\test2.faa"
#     # print(thetapi_calculate(allpath))
#     # print(thetapi_calculate(alignmentpath1))
#     # print(thetapi_calculate(alignmentpath2))
#     # print(snp_count(alignmentpath1))
#     print(thetapi_calculate(alignmentpath1))
#     print(thetaw_calculate(alignmentpath1))
#     print(calculate_tajimad(alignmentpath1))
#     print(fst_calculate(allpath,alignmentpath1,alignmentpath2))
#     print(Nm_calculate(allpath,alignmentpath1,alignmentpath2))