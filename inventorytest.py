from fuzzynumber import *
import dbpso

def leftFunction(x):
    return 1-x
def rightFunction(x):
    return 1-x

def minimizingFunction(space, c=3, d=500, r=0.8, c3=300, h=1.5, a=100, b=1.22):
    T = space[0]
    return (c*d/r) + (c3/T) + (h*d*T*(a*r+(b*r-1)*d))/(2*(a+b*d)*r)

d = LRFN(leftFunction, rightFunction, 500, 40, 100)
def func(space):
    itff = minimizingFunction(space, d=d)
    return itff.defuzz()

print("RUNNNING!")

solution = dbpso.DBPSO(func, Maxgen=100, N = 100, dim = 1, minx = [0.1], maxx = [3])
sol = solution.solve()
T = sol[0][0]
Z = minimizingFunction([T], d=d)
k = 100 + 1.22*d
t1 = 500*T/(0.8*k)
print(500*T, 0.8*k)
print("Fuzzy t1 is ", t1)
print("Fuzzy Production Rate is ", k)
Q = d*(T-t1)
print("Fuzzy Inventory Lvl is ", Q)
print("Fuzzy Cost is ", Z)
zplt = Z.toFuzzySet(Z.center-Z.alpha, Z.center+Z.beta).plot()
plt.show()
print("Fuzzy Minimizing in T and defuzzifying\nd = %f | t1 = %f | T = %f | k = %f | Q = %f | Z = %f"
      %(d.defuzz(), t1.defuzz(), T, k.defuzz(), Q.defuzz(), Z.defuzz()))