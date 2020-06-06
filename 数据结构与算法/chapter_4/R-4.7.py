#循环算法
def string2int(string):
    strint=''
    for i in string:
        if i in [str(i) for i in range(10)]:
            strint+=i
    return int(strint)
#递归算法