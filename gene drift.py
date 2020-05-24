import random
import matplotlib.pyplot as plt
plt.figure(figsize=(10,8),dpi=80)
for a in [i for i in range(1,50)]:
    initial = ['A1'] * 10 + ['A2'] * 10
    frequency = []
    times = [i for i in range(1, 100)]
    for i in times:
        next_generation = []
        while len(next_generation)<20:
            next_generation.append(initial[random.randint(0,19)])
        frequency.append(next_generation.count('A1') / 20)
        initial=next_generation
    plt.plot(times,frequency,label='times{}'.format(a))
plt.legend()
plt.show()