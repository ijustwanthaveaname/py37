def is_even(k):
    if not isinstance(k,int):
        raise TypeError
    else:
        k=abs(k)
        while abs(k)>1:
            k-=2
        return k==0