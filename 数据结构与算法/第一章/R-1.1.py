def is_multiple(n,m):
    if not isinstance(n,int) and not isinstance(m,int):
        raise TypeError
    elif n%m==0:
        return True
    else:
        return False
