class Progression:
    def __init__(self,start=0):
        self._current=start
    def _adavance(self):
        self._current+=1
    def __next__(self):
        if self._current is None:
            raise StopIteration
        else:
            answer=self._current
            self._adavance()
            return answer
    def __iter__(self):
        return self
    def print_progression(self,n):
        print(' '.join(str(next(self))for j in range(n)))

