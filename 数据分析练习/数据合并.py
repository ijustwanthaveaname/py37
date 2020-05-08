import pandas as pd
import numpy as np
df_A=pd.DataFrame(data=np.random.rand(3,3),index=[list("ABC")])
df_B=pd.DataFrame(data=np.random.random((3,3)),index=[list('BCD')])
print(df_A.append(df_B))
print(pd.concat(df_A,df_B,join="inner"))
print(pd.merge())




