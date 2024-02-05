from fuzzynumber import *
import dbpso

'''a  = float(input("Enter production constant a\n>"))
b = float(input("Enter production constant b\n>"))
r = float(input("Enter reliability\n>"))
c = float(input("Enter production cost per unit\n>"))
c3 = float(input("Enter setup cost per cycle\n>"))
h = float(input("Enter inventory cost per unit quantity per unit time\n>")) '''


def minimizingFunction(space, c=3, d=500, r=0.8, c3=300, h=1.5, a=100, b=1.22):
    T = space[0]
    return (c*d/r) + (c3/T) + (h*d*T*(a*r+(b*r-1)*d))/(2*(a+b*d)*r)

def minimizingFunctionForm2(space, c=3, d=500, r=0.8, c3=300, h=1.5, a=100, b=1.22):
    t1 = space[0]
    T = space[1]
    k = a+b*d
    ans = (c*k*t1 + c3 + h*(0.5*(r*k-d)*t1**2 + d*T*(T-t1) + 0.5*d*(t1**2-T**2)))/T
    return ans

# def minimizingFunctionForm3(space, c=3, d=500, r=0.8, c3=300, h=1.5, a=100, b=1.22):
#     k = a + b * d
#     Q = space[0]
#     T = space[1]
#     ans = (1/T)*((c*k*Q)/(r*k-d) + c3 + (h*r*k*Q**2)/(2*d*(r*k-d)))
#     return ans


solution = dbpso.DBPSO(minimizingFunction, Maxgen=100, N = 100, dim = 1, minx = [0.2], maxx = [2])
sol = solution.solve()
T = sol[0][0]
Z = sol[1]
k = 100 + 1.22*500
t1 = T*500/(0.8*k)
Q = 500*(T-t1)
print("Minimizing in T\nt1 = %f | T = %f | k = %f | Q = %f | Z = %f" %(t1, T, k, Q, Z))

#FUZZY IMPLEMENTATION
print('\n')
d = tfn(460, 500, 600)
print("Fuzzy Demand is ", d)

def func(space):
    tff = minimizingFunction(space, d=d)
    return tff.yRankingIndex()

solution = dbpso.DBPSO(func, Maxgen=100, N = 100, dim = 1, minx = [0.1], maxx = [3])
sol = solution.solve()
T = sol[0][0]
Z = minimizingFunction([T], d=d)
k = 100 + 1.22*d
t1 = 500*T/(0.8*k)
print("Fuzzy t1 is ", t1)
print("Fuzzy Production Rate is ", k)
Q = d*(T-t1)
print("Fuzzy Inventory Lvl is ", Q)
print("Fuzzy Cost is ", Z)
zplt = Z.plot()
plt.show()
print("Fuzzy Minimizing in T and defuzzifying\nd = %f | DOF = %f | t1 = %f | T = %f | k = %f | Q = %f | Z = %f"
      %(d.yRankingIndex(), d.dof(), t1.yRankingIndex(), T, k.yRankingIndex(), Q.yRankingIndex(), Z.yRankingIndex()))


#INTUITIONISTIC FUZZY NUMBERS
d = itfn(460, 400, 500, 600, 660)
print("Intuitionistic Fuzzy Demand is ", d)
def func(space):
    itff = minimizingFunction(space, d=d)
    return itff.defuzz()

solution = dbpso.DBPSO(func, Maxgen=100, N = 100, dim = 1, minx = [0.1], maxx = [3])
sol = solution.solve()
T = sol[0][0]
Z = minimizingFunction([T], d=d)
k = 100 + 1.22*d
t1 = 500*T/(0.8*k)
print("Int. Fuzzy t1 is ", t1)
print("Int. Fuzzy Production Rate is ", k)
Q = d*(T-t1)
print("Int. Fuzzy Inventory Lvl is ", Q)
print("Int. Fuzzy Cost is ", Z)
zplt = Z.plot()
plt.show()
print("Int. Fuzzy Minimizing in T and defuzzifying\nd = %f | t1 = %f | T = %f | k = %f | Q = %f | Z = %f"
      %(d.defuzz(), t1.defuzz(), T, k.defuzz(), Q.defuzz(), Z.defuzz()))

#PENALTY METHOD
