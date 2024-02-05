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
    def __init__(self, elements):
        self.elements = elements
        [self.names, self.values] = [[elements[i].name for i in range(len(elements))], [elements[i].value for i in range(len(elements))]]

    def __str__(self):
        string = ''
        for i in range(len(self.elements)):
            string += ('| ({}, {}) |'.format(self.elements[i].name, self.elements[i].value))
        return string + '\n'

    def __getitem__(self, index):
        for i in range(len(self.elements)):
            if self.elements[i].name == index:
                return self.elements[i].value

    def __setitem__(self, item, value):
        #TBC

    def __delitem__(self, item):
        if isinstance(item, int):
            del self.value[item]
        else:
            item = sorted(item, reverse=True)
            for elem in item:
                del self.value[elem]

    def plot(self):
        f =  plt.figure()
        plt.plot(self.names, self.values)
        return f