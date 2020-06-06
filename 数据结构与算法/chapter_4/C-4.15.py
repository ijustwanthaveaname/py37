def subset(S):
    i=0
    while i<=len(S)-1:
        raw = list(S)
        del raw[i]
        print(raw)