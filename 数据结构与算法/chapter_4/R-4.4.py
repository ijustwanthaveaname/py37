def reverse(S, start, stop):
    if start < stop:
        S[start], S[stop] = S[stop], S[start]
        reverse(S, start + 1, stop - 1)
#reverse(S,0,5) S[0],S[5]=S[5],S[0]
##reverse(S,1,4) S[1].S[4]=S[4],S[1]
###reverse(S,2,3) S[2],S[3]=S[3],S[2]


