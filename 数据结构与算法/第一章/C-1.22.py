def dot_product(a,b):
    c=[]
    if len(a) != len(b):
        print('dot product must have same length')
    else:
        for i  in range(len(a)):
            c.append(a[i]*b[i])
    return c



