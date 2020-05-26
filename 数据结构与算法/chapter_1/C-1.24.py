def numOfCh(string):
    temp = {'a', 'e', 'i', 'o', 'u', 'A', 'E', 'I', 'O', 'U'}
    count = 0
    for item in string:
        if item in temp:
            count += 1
    return count
