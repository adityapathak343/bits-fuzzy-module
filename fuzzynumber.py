from fuzzysets import *
class intNum:
    def __init__(self, lowLim:float, upLim:float):
        self.lowLim = lowLim
        self.upLim = upLim

    def __str__(self):
        return "INTNUM [{}, {}]".format(str(self.lowLim), str(self.upLim))

    def __neg__(self):
        return intNum(-self.upLim, -self.lowLim)

    def __add__(self, rhs):
        return intNum(self.lowLim+rhs.lowLim, self.upLim+rhs.upLim)

    def __sub__(self, rhs):
        return intNum(self.lowLim-rhs.upLim, self.upLim-rhs.lowLim)

    def __mul__(self, rhs):
        if type(rhs) == float or type(rhs) == int:
            return intNum(self.lowLim*rhs, self.upLim*rhs)
        rho = [self.lowLim*rhs.lowLim, self.lowLim*rhs.upLim, self.upLim*rhs.lowLim, self.upLim*rhs.upLim]
        return intNum(min(rho), max(rho))

    def __rmul__(self, rhs):
        return self*rhs

    def __truediv__(self, other):
        if type(other) == float or type(other) == int:
            return intNum(self.lowLim/other, self.upLim/other)
        return self * intNum(1/other.upLim, 1/other.lowLim)

    def __rtruediv__(self, other):
        return other * intNum(1/self.upLim, 1/self.lowLim)

    def __pow__(self, other):
        try:
            factor = float(other)
            return intNum(self.lowLim**factor, self.upLim**factor)
        except:
            raise NotImplemented
    def __rpow__(self, other):
        try:
            factor = float(other)
            return intNum(factor**self.lowLim, factor**self.upLim)
        except:
            raise NotImplemented

    def __and__(self, rhs):
        if max(self.lowLim, rhs.lowLim) <= min(self.upLim, rhs.upLim):
            return intNum(max(rhs.lowLim, self.lowLim), min(rhs.upLim, self.upLim))

    def __rand__(self, rhs):
        return self&rhs

    def __or__(self, rhs):
        if max(self.lowLim, rhs.lowLim) <= min(self.upLim, rhs.upLim):
            return intNum(min(rhs.lowLim, self.lowLim), max(rhs.upLim, self.upLim))
    def __ror__ (self, rhs):
        return self|rhs

    def inNum(self, x):
        return self.lowLim <= x <= self.upLim

    def toFuzzySet(self, lowLim, upLim, N=1001):
        h = 1 / N
        elements = [FuzzyElement(lowLim + i * h * (upLim - lowLim), float(self.inNum(lowLim + i * h * (upLim - lowLim)))) for i in
                    range(N + 1)]
        return FuzzySet(elements)



# class FuzzyNumber:
#     def __init__(self):


class TFN:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def __str__(self):
        if type(self.a) == type(self):
            return "TFN <{}, {}, {}>".format(str(self.a), str(self.b), str(self.c))
        return "TFN <%f, %f, %f>" % (self.a, self.b, self.c)

    def __neg__(self):
        return TFN(-self.c, -self.b, -self.a)

    def __add__(self, rhs):
        try:
            rhs = float(rhs)
            return TFN(self.a + rhs, self.b + rhs, self.c + rhs)
        except:
            return TFN(self.a + rhs.a, self.b + rhs.b, self.c + rhs.c)

    def __radd__(self, rhs):
        return self+rhs

    def __sub__(self, rhs):
        try:
            rhs = float(rhs)
            return TFN(self.a - rhs, self.b - rhs, self.c - rhs)
        except:
            return TFN(self.a - rhs.c, self.b - rhs.b, self.c - rhs.c)

    def __rsub__(self, other):
        return -(self-other)

    def __eq__(self, rhs):
        return self.a == rhs.a and self.b == rhs.b and self.c == rhs.c

    def __mul__(self, other):
        try:
            factor = float(other)
            return TFN(self.a * factor, self.b * factor, self.c * factor)
        except TypeError:
            if type(self) == type(other):
                return TFN(self.a * other.a, self.b * other.b, self.c * other.c)
            raise NotImplemented


    def __rmul__(self, other):
        return self*other

    def __pow__(self, power, modulo=None):
        if type(power)!=type(5):
            raise NotImplemented
        if power > 0:
            return self*(self**(power-1))
        if power == 0:
            return TFN(1, 1, 1)
        if power < 0:
            return 1/(self**-power)

    def __truediv__(self, other):
        try:
            factor = float(other)
            return TFN(self.a / factor, self.b / factor, self.c / factor)
        except TypeError:
            if type(self) == type(other):
                return TFN(self.a / other.c, self.b / other.b, self.c / other.a)
            raise NotImplemented

    def __rtruediv__(self, other):
        try:
            factor = float(other)
            return TFN(factor / self.c, factor / self.b, factor / self.a)
        except TypeError:
            if type(self) == type(other):
                return other/self
            raise NotImplemented

    def mu(self, x: float):
        if x<self.a:
            return 0
        elif self.a < x < self.b:
            return (x-self.a)/(self.b-self.a)
        elif self.b < x < self.c:
            return (self.c-x)/(self.c-self.b)
        else:
            return 0

    def alphaCut(self, alpha: float):
        if alpha<0:
            return intNum(self.a, self.c)
        elif alpha>1:
            return 0
        else:
            lowLim = self.a + alpha*(self.b-self.a)
            upLim = self.c - alpha*(self.c-self.b)
            return intNum(lowLim, upLim)
    
    def yRankingIndex(self):
        return (self.a + 2*self.b + self.c)/4

    def mean(self):
        return (self.a + self.b + self.c)/3

    def mode(self):
        return 3*self.b -2*self.mean()

    def dof(self):
        return (self.c-self.a)/self.mode()

    def plot(self, xlabel='', ylabel='', N=1000, deev=False):
        N = N if deev else N+1
        h = 1/N
        x = [self.a + (self.c - self.a)*i*h for i in range(N+1)]
        y = [self.mu(i) for i in x]
        f = plt.figure()
        plt.plot(x, y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        return f

    def toFuzzySet(self, lowLim, upLim, N=1001):
        h = 1/N
        elements = [FuzzyElement(lowLim + i*h*(upLim- lowLim), self.mu(lowLim + i*h*(upLim- lowLim))) for i in range(N+1)]
        return FuzzySet(elements)



class ITFN:
    def __init__(self, a1, a2, b, c1, c2):
        self.a1 = a1
        self.a2 = a2
        self.b = b
        self.c1 = c1
        self.c2 = c2

    def __str__(self):
        if type(self.b) == type(self):
            return "ITFN <{}, {}, {}>; <{}, {}, {}>".format(
                str(self.a1), str(self.b), str(self.c1), str(self.a2), str(self.b), str(self.c2))
        return "ITFN <%f, %f, %f>; <%f, %f, %f>" % (self.a1, self.b, self.c1, self.a2, self.b, self.c2)

    def __neg__(self):
        return ITFN(-self.c1, -self.c2, -self.b, -self.a1, -self.a2)

    def __add__(self, rhs):
        try:
            rhs = float(rhs)
            return ITFN(self.a1 + rhs, self.a2 + rhs, self.b + rhs, self.c1 + rhs, self.c2 + rhs)
        except:
            return ITFN(self.a1 + rhs.a1, self.a2 + rhs.a2, self.b + rhs.b, self.c1 + rhs.c1, self.c2 + rhs.c2)

    def __radd__(self, rhs):
        return self + rhs

    def __sub__(self, rhs):
        try:
            rhs = float(rhs)
            return self + (-rhs)
        except TypeError:
            return self + (-rhs)

    def __rsub__(self, other):
        return -(self - other)

    def __eq__(self, rhs):
        return self.a1 == rhs.a1 and self.a2 == rhs.a2 and self.b == rhs.b and self.c1 == rhs.c1 and self.c2 == rhs.c2

    def __mul__(self, other):
        try:
            factor = float(other)
            return ITFN(self.a1 * factor, self.a2 * factor, self.b * factor, self.c1 * factor, self.c2 * factor)
        except TypeError:
            if type(self) == type(other):
                return ITFN(self.a1 * other.a1, self.a2 * other.a2, self.b * other.b, self.c1 * other.c1, self.c2 * other.c2)
            raise NotImplemented

    def __rmul__(self, other):
        return self * other

    def __pow__(self, power, modulo=None):
        if type(power) != type(5):
            raise NotImplemented
        if power > 0:
            return self * (self ** (power - 1))
        if power == 0:
            return ITFN(1, 1, 1, 1, 1)
        if power < 0:
            return 1 / (self ** -power)

    def __truediv__(self, other):
        try:
            factor = float(other)
            return self * ITFN(1 / factor, 1 / factor, 1 / factor, 1 / factor, 1 / factor)
        except TypeError:
            if type(self) == type(other):
                return ITFN(self.a1 / other.c1, self.a2 / other.c2, self.b / other.b, self.c1 / other.a1, self.c2 / other.a2)
            raise NotImplemented

    def __rtruediv__(self, other):
        try:
            factor = float(other)
            return factor * (ITFN(1, 1, 1, 1, 1) / self)
        except TypeError:
            if type(self) == type(other):
                return other / self
            raise NotImplemented

    def mu(self, x: float):
        if x < self.a1:
            return 0
        elif self.a1 < x < self.b:
            return (x - self.a1) / (self.b - self.a1)
        elif self.b < x < self.c1:
            return (self.c1 - x) / (self.c1 - self.b)
        else:
            return 0

    def neta(self, x: float):
        if x < self.a2:
            return 1
        elif self.a2 < x < self.b:
            return (self.b - x) / (self.b - self.a2)
        elif self.b < x < self.c2:
            return (x - self.b) / (self.c2 - self.b)
        else:
            return 1

    def alphaBetaCut(self, alpha: float, beta: float):
        intList = []
        if alpha < 0:
            intList[0] = intNum(self.a1, self.c1)
        elif alpha > 1:
            intList[0] = 0
        else:
            lowLim = self.a1 + alpha * (self.b - self.a1)
            upLim = self.c1 - alpha * (self.c1 - self.b)
            intList[0] = intNum(lowLim, upLim)
        if beta < 0:
            intList[1] = intNum(self.a2, self.c2)
        elif alpha > 1:
            intList[1] = 0
        else:
            lowLim = self.b - beta * (self.b - self.a2)
            upLim = self.b + beta * (self.c2 - self.b)
            intList[1] = intNum(lowLim, upLim)
        return intList

    def defuzz(self):
        return 1/8 * (self.a1 + self.a2 + 4*self.b + self.c1 + self.c2)

    def plot(self, xlabel='', ylabel = '', N=1000, deev=False):
        N = N if deev else N + 1
        h = 1/N
        x = [min(self.a1, self.a2) + (max(self.c1, self.c2) - min(self.a1, self.a2))*i*h for i in range(N+1)]
        y1 = [self.mu(i) for i in x]
        y2 = [self.neta(i) for i in x]
        try:
            f = plt.figure()
            plt.plot(x, y1)
            plt.plot(x, y2)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            return f
        except ModuleNotFoundError:
            print("Matplotlib is needed to use the plot feature")

class TrFN:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d


    def __neg__(self):
        return TrFN(-self.d, -self.c, -self.b, -self.a)
    def __add__(self, rhs):
        try:
            rhs = float(rhs)
            return TrFN(self.a + rhs, self.b + rhs, self.c + rhs, self.d+rhs)
        except:
            return TrFN(self.a + rhs.a, self.b + rhs.b, self.c + rhs.c, self.d+rhs.d)

    def __radd__(self, rhs):
        return self+rhs

    def __sub__(self, rhs):
        try:
            rhs = float(rhs)
            return TrFN(self.a - rhs, self.b - rhs, self.c - rhs)
        except:
            return TrFN(self.a - rhs.d, self.b - rhs.c, self.c - rhs.b, self.d - rhs.a)

    def __rsub__(self, other):
        return -(self-other)

    def __eq__(self, rhs):
        return self.a == rhs.a and self.b == rhs.b and self.c == rhs.c and self.d == rhs.d

    def __mul__(self, other):
        try:
            factor = float(other)
            return TrFN(self.a * factor, self.b * factor, self.c * factor, self.d * factor)
        except TypeError:
            if type(self) == type(other):
                return TrFN(self.a * other.a, self.b * other.b, self.c * other.c, self.d * other.d)
            raise NotImplemented


    def __rmul__(self, other):
        return self*other

    def __pow__(self, power, modulo=None):
        if type(power)!=int:
            raise NotImplemented
        if power > 0:
            return self*(self**(power-1))
        if power == 0:
            return TrFN(1, 1, 1, 1)
        if power < 0:
            return 1/(self**-power)

    def __truediv__(self, other):
        try:
            factor = float(other)
            return TFN(self.a / factor, self.b / factor, self.c / factor)
        except TypeError:
            if type(self) == type(other):
                return TrFN(self.a / other.d, self.b / other.c, self.c / other.b, self.d / other.a)
            raise NotImplemented

    def __rtruediv__(self, other):
        try:
            factor = float(other)
            return TrFN(factor/self.d, factor / self.c, factor / self.b, factor / self.a)
        except TypeError:
            return other/self
    def mu(self, x):
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        if a<=x<b:
            return (x-a)/(b-a)
        elif b<=x<c:
            return 1
        elif c<=x<d:
            return (d-x)/(d-c)
        return 0
    def alphaCut(self, alpha):
        if alpha<=0:
            return intNum(self.a, self.d)
        elif alpha==1:
            return intNum(self.b, self.c)
        lowLim = self.a + alpha*(self.b-self.a)
        upLim = self.d - alpha*(self.d-self.c)
        return intNum(lowLim, upLim)

    def toFuzzySet(self, lowLim, upLim, N=1001):
        h = 1/N
        elements = [FuzzyElement(lowLim + i*h*(upLim- lowLim), self.mu(lowLim + i*h*(upLim- lowLim))) for i in range(N+1)]
        return FuzzySet(elements)
