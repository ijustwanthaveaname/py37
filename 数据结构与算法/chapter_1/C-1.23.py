def getvalue(data, i):
    try:
        return data[i]
    except IndexError:
        print("Don't try buffer overflow attacks in Python!")
