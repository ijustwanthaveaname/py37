def function31():
    coins = {'0.5': 0, '1': 0, '5': 0, '10': 0, '20': 0, '50': 0, '100': 0}
    temp = input("Please input Pay and Total money: \n").split(" ")
    pay, total = int(temp[0]), int(temp[1])
    rest = total - pay

    coins[100] = int(rest / 100)
    rest = int(rest % 100)
    coins[50] = int(rest / 50)
    rest = int(rest % 50)
    coins[20] = int(rest / 20)
    rest = int(rest % 20)
    coins[10] = int(rest / 10)
    rest = int(rest % 10)
    coins[5] = int(rest / 5)
    rest = int(rest % 5)
    coins[1] = int(rest / 1)
    rest = int(rest % 1)
    coins[0.5] = int(rest / 0.5)
    rest = int(rest % 0.5)

    return coins.values()
