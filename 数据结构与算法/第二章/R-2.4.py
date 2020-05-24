class Flower:
    def __init__(self,name="mudan",number=20,price=20.5):
        self._name=name
        self._number=number
        self._price=price

    def get_name(self):
        return self._name

    def set_name(self,names):
        self._name=names

    def get_number(self):
        return self._number

    def set_number(self,numbers):
        self._number=numbers

    def get_price(self):
        return self._price

    def set_price(self,prices):
        self._price=prices