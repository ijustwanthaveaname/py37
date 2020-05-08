from itertools import permutations
def align(n):
    return list(map(''.join,permutations(n)))