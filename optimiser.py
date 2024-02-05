def derivative(func, x, N=1000000):
    h = 1/N
    return (func(x+h)-func(x))/h

def integral(func, t0, t, N=100000):
    h = 1/N
    line = [(t0+(t-t0)*i*h) for i in range(N+1)]
    sum = 0
    for i in line:
        sum += func(i)*h
    return sum

def xsq(x):
    return x**2

print(derivative(xsq, 2))
print(integral(xsq, 0, 2))
