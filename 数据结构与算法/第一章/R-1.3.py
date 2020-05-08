def minmax(data):
    minv = data[0]
    maxv = data[0]
    for item in data:
        if minv > item:
            minv = item
        if maxv < item:
            maxv = item
    return minv, maxv
