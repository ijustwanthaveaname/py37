import pandas as pd
import numpy as np
data=np.random.randint(0,50,24).reshape(4,6)
index=pd.MultiIndex.from_product([[2013,2014],[1,2]])
column=pd.MultiIndex.from_product([["Bob","guido","Sue"],["HR","Temp"]])

df=pd.DataFrame(data=data,columns=column,index=index)
df.index.names = ["years", "visit"]
df.columns.names = ["subject", "type"]
idx=pd.IndexSlice
print(df.loc[idx[:,1],idx["Bob","HR"]])
print(df)
#将多级索引返回为列:df.reset_index()
#将列返回为索引df.reset_index().set_index(["years","visit"])
df.mean(axis=1,level="type")