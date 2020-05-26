def sumsquare(n):
    if not isinstance(n,int):
        raise TypeError
    elif n<=0:
        print('Please input positive integer')
    else:
        i=1
        sum=0
        while i<=n:
            sum+=i**2
            i+=1
        return sum
