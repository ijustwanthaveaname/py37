def oddsumsquare(n):
    if not isinstance(n, int):
        raise TypeError
    elif n <= 0:
        raise ValueError
    else:
        square=[i**2 for i in range(n+1) if i%2!=0]
        return sum(square)