def norm(v,p):
    if p<1:
        raise ValueError
    box=[]
    for i in v:
        box.append(i**p)
    return sum(box)**(1/p)
# def norm(v, p=2):
#     import math
#     return math.sqrt(sum(pow(x, p) for x in v))