from matplotlib import pyplot
filepath=input("Please input address of file :")
file=''.join(open(filepath,'r').readlines())
x=[chr(i) for i in range(65,91)]+[chr(i) for i in range(97,123)]
y=[]
for i in x:
    y.append(file.count(i))
pyplot.bar(x,y)

