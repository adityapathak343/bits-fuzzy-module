try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib is not installed. Plotting features will not work.")


class FuzzyElement:
    def __init__(self, name, value:float):
        self.name = name
        if 0<=value<=1:
            self.value = value
        else:
            self.value = None


class FuzzySet:
    def __init__(self, elements=[]):
        self.elements = elements
        [self.names, self.values] = [[elements[i].name for i in range(len(elements))], [elements[i].value for i in range(len(elements))]]

    def __str__(self):
        string = ''
        for i in range(len(self.elements)):
            string += ('| ({}, {}) |'.format(self.elements[i].name, self.elements[i].value))
        return string + '\n'

    def __getitem__(self, index):
        if index in self.names:
            return self.values[self.names.index(index)]
        return 0

    def __setitem__(self, item, value):
        if item in self.names:
            i = self.names.index(item)
            self.elements[i].value = value
        else:
            self.elements.append(FuzzyElement(item, value))
        self.__init__(elements=self.elements)

    def __delitem__(self, item):
        if item in self.names:
            i = self.names.index(item)
            del self.elements[i]
            self.__init__(elements=self.elements)

    def __mul__(self, other):
        if type(other) == float and 0<=other<=1:
            return FuzzySet([FuzzyElement(i.name, other*i.value) for i in self.elements])

    def __rmul__(self, other):
        return self*other
    def __and__(self, other):
        w = FuzzySet([])
        for i in self.elements:
            w.elements.append(FuzzyElement(i.name, min(i.value, other[i.name])))
        w.__init__(w.elements)
        return w

    def __rand__(self, other):
        return self&other

    def __or__(self, other):
        w = []
        for i in self.elements:
            w.append(FuzzyElement(i.name, max(i.value, other[i.name])))
        return FuzzySet(w)

    def __ror__(self, other):
        return self|other
    def alphaCut(self, alpha):
        retList = []
        for i in self.elements:
            if i.value < alpha:
                retList.append(FuzzyElement(i.name, 0))
            else:
                retList.append(FuzzyElement(i.name, 1))
        return FuzzySet(retList)


    def plot(self):
        f =  plt.figure()
        plt.scatter(self.names, self.values)
        return f