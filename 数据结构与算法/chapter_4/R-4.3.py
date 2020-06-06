def power(x,n):
    if n==0:
        return 1
    else:
        partial=power(x,n//2)
        result=partial*partial
        if n % 2==1:
            result*=x
        return result

#power(2,18)partial=power(2,9)=2*16*16 result=2*16*16*16*16*2 return result=2*16*16*16*16*2
##power(2,9) partial=power(2,4)=16 result=16*16 because 9%1==1 result==2*16*16 return 2*16*16
###power(2,4) partial=power(2,2)=4 result=4*4 return 16
####power(2,2) partial=power(2,1)=2 result=2*2 return 4
#####power(2,1) partial=power(2,0)=1 result=1*1 because 1%2==1 result=1*2=2 return 2
######power(2,0) return 1