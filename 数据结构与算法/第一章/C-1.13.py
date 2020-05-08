def reverse(ls):
    a=[]
    b=0
    for i in ls:
        if not isinstance(i,int):
            raise TypeError
        b-=1
        a.append(ls[b])
    return a
#方法二
# def reverse(ls):
#     return ls[::-1]
