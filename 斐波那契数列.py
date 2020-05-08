def fibonacci(n,a=1,b=1,i=1):
    print(1,1)
    while i<=n/2-1:
        i=i+1
        a=a+b
        b=a+b
        print(a,b)
    if n%2!=0:
        a=a+b
        print(a)
    else:
        pass
fibonacci(n=9)