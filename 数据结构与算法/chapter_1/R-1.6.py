def oddsumsquare(n):
    if not isinstance(n,int):
        raise TypeError
    elif n<=0 or n%2==0:
        raise ValueError
    else:
        i=1
        sum=0
        while i<=n:
            sum+=i**2
            i+=2
        return sum