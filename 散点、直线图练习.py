import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
a=range(1,11)
a_label=list("abcdefghij")
b=range(0,10)
fig,ax=plt.subplots(2,2)
ax[0,0].scatter(a,b,s=5)
ax[0,0].plot(b,a,"o",markersize=1)
ax[0,0].set_xticks(np.arange(len(a)))
ax[0,0].set_xticklabels(a_label)
ax[0,0].set_yticks(np.arange(len(b)))
ax[0,0].set_title("fuck you")
fig.suptitle("main title")
plt.show()
