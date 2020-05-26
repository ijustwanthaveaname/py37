from random import randint
def shuffle(data):
    indbox=[]
    shufflelist=[]
    while len(shufflelist) < len(data):
        ind=randint(0,len(data)-1)
        if ind not in indbox :
            shufflelist.append(data[ind])
            indbox.append(ind)
    return shufflelist



