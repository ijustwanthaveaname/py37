def diffnum(data):
    ind=0
    for i in data:
        if not isinstance(i,(int,float)):
            raise TypeError
        ind+=1
        for n in data[ind:]:
            if i==n:
                return False
    return True


