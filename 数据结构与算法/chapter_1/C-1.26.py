def function():
    temp = input("Please input three numbers!\n")
    nums = list(temp.split(" "))
    a = int(nums[0])
    b = int(nums[1])
    c = int(nums[2])
    if (a + b == c) or (a == b - c) or ( a * b == c):
        return True
    else:
        return False
