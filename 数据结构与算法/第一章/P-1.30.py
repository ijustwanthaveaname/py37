def function30():
    num = int(input("Please input a number: \n"))
    count = 0
    rest = num
    while rest >= 2:
        rest = int(rest / 2)
        count += 1
    return count