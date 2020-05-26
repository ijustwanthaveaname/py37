def hasoddmulti(data):
    ind=0
    for i in data:
        if not isinstance(i,int):
            raise TypeError
        ind+=1
        for n in data[ind:]:
            if n!=i and (n*i%2)==1:
                return True
            else:
                return False