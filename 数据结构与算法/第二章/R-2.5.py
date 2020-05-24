class Vector:
    def __init__(self,d):
        if isinstance(d,int):
            self._coords=[0]*d
        elif isinstance(d,list):
            self._coords=[i for i in d]
        else:
            raise TypeError
        
    def __len__(self):
        return len(self._coords)

    def __getitem__(self,j):
        return self._coords[j]

    def __setitem__(self,j,val):
        self._coords[j]=val

    # def add(self,other):
    #     if len(self)!=len(other):
    #         raise ValueError('dimensions must agree')
    #     result=Vector(len(self))
    #     for j in range(len(self)):
    #         result[j]=self[j]+other[j]
    #     return result
    def __add__(self,other):
        if len(self) != len(other):
            raise ValueError('dimensions must agree')
        result = Vector(len(self))
        for j in range(len(self)):
            result[j] = self[j] + other[j]
        return result

    def __radd__(self, other):
        if len(self) != len(other):
            raise ValueError('dimensions must agree')
        result = Vector(len(self))
        for j in range(len(self)):
            result[j] = self[j] + other[j]
        return result

    def __sub__(self,other):
        if len(self)!=len(other):
            raise ValueError('dimensions must agree')
        result=Vector(len(self))
        for j in range(len(self)):
            result[j]=self[j]-other[j]
        return result

    # def mul(self,n):
    #     # result = Vector(len(self))
    #     if isinstance(n,Vector) and len(self)==len(n):
    #         result=0
    #         for j in range(len(self)):
    #             result+=self[j]*n[j]
    #     else:
    #         result=Vector(len(self))
    #         for j in range(len(self)):
    #             result[j] = self[j] * n
    #     return result
    def __mul__(self,n):
        if isinstance(n, Vector) and len(self) == len(n):
            result = 0
            for j in range(len(self)):
                result += self[j] * n[j]
        else:
            result = Vector(len(self))
            for j in range(len(self)):
                result[j] = self[j] * n
        return result

    def __rmul__(self, n):
        if isinstance(n, Vector) and len(self) == len(n):
            result = 0
            for j in range(len(self)):
                result += self[j] * n[j]
        else:
            result = Vector(len(self))
            for j in range(len(self)):
                result[j] = self[j] * n
        return result

    def __neg__(self):
        result=Vector(len(self))
        for j in range(len(self)):
            result[j]=-abs(self[j])
        return result

    def __eq__(self,other):
        return self._coords==other._coords

    def __ne__(self,other):
        return not self==other #rely on existing __eq__ definition

    def __str__(self):
        return '<'+str(self._coords)[1:-1]+'>'




