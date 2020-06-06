def find_max(arr):
    tmp = arr.pop(0)
    if len(arr) == 0:
        return tmp
    max = find_max(arr)
    if max > tmp:
        return max
    else:
        return tmp
"""
find_max([2,1,4]);tmp=2;arr=[1,4];max=find_max([1,4]{tmp=1;arr=[4];max=find_max(4){return 4}=4;return 4}=4
"""

