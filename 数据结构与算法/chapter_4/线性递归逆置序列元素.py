def reverse(S,start,stop):
    if start<stop:
        S[start],S[stop]=S[stop],S[start]
        reverse(S,start+1,stop-1)