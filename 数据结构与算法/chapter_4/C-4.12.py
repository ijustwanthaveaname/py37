def mutiply(m,n):
    if n==0 or m==0:
        return 0
    return n+mutiply(m-1,n)
