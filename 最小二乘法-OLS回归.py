import numpy as np
def ols(A,row,col,b):
    A=np.matrix(A).reshape(row,col)
    X=((A.transpose()*A).I)*A.transpose()*np.matrix(b).transpose()
    return X

