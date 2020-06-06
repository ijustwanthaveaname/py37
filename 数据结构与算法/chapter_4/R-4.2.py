def power(x,n):
    if n==0:
        return 1
    else:
        return x*power(x,n-1)

###power(2,5) return 2*2*2*2*2*1
####power(2,4) return 2*2*2*2*1
#####power(2,3) return 2*2*2*1
######power(2,2) return 2*2*1
#######power(2,1) return 2*1
########power(2,0) return 1