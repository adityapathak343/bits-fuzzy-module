from fuzzynumber import *
import dbpso
from bees_algorithm import BeesAlgorithm
import numpy as np

def cosh(x):
    return np.cosh(x)

def complex_equation(space):
    T1 = space[0]
    T = space[1]
    if T1>T:
        return -10000000000
    a = 60
    b = 0.1
    p = 100
    c = 4
    m1 = 7
    delta = 0.01
    A = 1000
    omega = 120
    psi = 0.0015
    C0 = 12
    C1 = 2
    C2 = 15
    C3 = 4
    C4 = 10
    Z = 200
    lambda_val = 0.3
    xi = 0.2
    G = 20
    m = 0.8
    k1 = 2
    k2 = 3
    k3 = 1.5
    k4 = 4

    term1 = A

    numerator1 = np.exp(psi * omega) * (m1 - T1 + 1) ** (-np.exp(-psi * omega)) * (
            (m1 + 1) * (a - b * p + 2 * np.exp(2 * psi * omega) * (c * (m1 + 1) ** 2 + 3 * a - 3 * b * p) - np.exp(
        psi * omega) * (5 * a - 5 * b * p)) * (m1 - T1 + 1) ** np.exp(-psi * omega) -
            (m1 + 1) ** np.exp(-psi * omega) * (m1 - T1 + 1) * (c * T1 ** 2 + a - b * p - np.exp(psi * omega) * (
                5 * a - 5 * b * p + c * T1 * (2 * m1 + 3 * T1 + 2)) +
                                                                2 * np.exp(2 * psi * omega) * (
                                                                            3 * a - 3 * b * p + c * (
                                                                                (m1 + 1) ** 2 + T1 * (
                                                                                    m1 + 1) + T1 ** 2)))
    )
    denominator1 = np.exp(2 * psi * omega) * (12 * cosh(psi * omega) - 11) - 1
    term2 = (C0 * numerator1 / denominator1)

    term3 = (
                    -c * (T - T1) * delta * ((3 * T + T1) * delta + 2) +
                    2 * ((a - b * p) * delta ** 2 + c * (T * delta + 1) ** 2) * np.log(T * delta - T1 * delta + 1)
            ) / (2 * delta ** 3)

    term4 = (
                    ((-1 + np.exp(-G * m1)) * xi + 1) * (
                    k1 + k2 * np.exp(psi * omega) * (m1 - T1 + 1) ** (-np.exp(-psi * omega)) * (
                    (m1 + 1) * (a - b * p - 5 * np.exp(psi * omega) * (a - b * p) + 2 * np.exp(2 * psi * omega) * (
                        c * (m1 + 1) ** 2 + 3 * a - 3 * b * p)) * (m1 - T1 + 1) ** np.exp(-psi * omega) -
                    (m1 + 1) ** np.exp(-psi * omega) * (m1 - T1 + 1) * (
                                c * T1 ** 2 + a - b * p - np.exp(psi * omega) * (
                                    5 * a - 5 * b * p + c * T1 * (2 * m1 + 3 * T1 + 2)) +
                                2 * np.exp(2 * psi * omega) * (
                                            3 * a - 3 * b * p + c * (m1 ** 2 + (T1 + 2) * m1 + T1 ** 2 + T1 + 1)))
            )
            )
            ) / denominator1

    numerator2 = np.exp(psi * omega) * k3 * (m1 - T1 + 1) ** (-np.exp(-psi * omega)) * (
            -1 / 12 * T1 * (
            6 * b * (1 - 5 * np.exp(psi * omega) + 6 * np.exp(2 * psi * omega)) * p * (2 * m1 - T1 + 2) +
            6 * a * (1 - 5 * np.exp(psi * omega) + 6 * np.exp(2 * psi * omega)) * (-2 * m1 + T1 - 2) +
            c * ((-4 * m1 + 3 * T1 - 4) * T1 ** 2 - np.exp(psi * omega) * (
                -12 * (m1 + 1) ** 2 - 4 * T1 * (m1 + 1) + 9 * T1 ** 2) * T1 - 6 * np.exp(2 * psi * omega) * (
                             4 * (m1 + 1) ** 3 - T1 ** 3))
    ) * (m1 - T1 + 1) ** np.exp(-psi * omega) -
            np.exp(psi * omega) * (
                        (m1 + 1) ** (1 + np.exp(-psi * omega)) - (m1 - T1 + 1) ** (1 + np.exp(-psi * omega))) * (
                        m1 - T1 + 1) * (c * T1 ** 2 + a - b * p - np.exp(psi * omega) * (
                5 * a - 5 * b * p + c * T1 * (2 * m1 + 3 * T1 + 2)) +
                                        2 * np.exp(2 * psi * omega) * (3 * a - 3 * b * p + c * (
                        m1 ** 2 + (T1 + 2) * m1 + T1 ** 2 + T1 + 1))
                                        )) / (1 + np.exp(psi * omega))

    term5 = numerator2 / denominator1

    term6 = (
                    -3 * np.exp(psi * omega) * (
                        (m1 + 1) ** np.exp(-psi * omega) - (m1 - T1 + 1) ** np.exp(-psi * omega)) * (
                                c * T1 ** 2 + a - b * p - np.exp(psi * omega) * (
                                    5 * a - 5 * b * p + c * T1 * (2 * m1 + 3 * T1 + 2)) +
                                2 * np.exp(2 * psi * omega) * (3 * a - 3 * b * p + c * (
                                    m1 ** 2 + (T1 + 2) * m1 + T1 ** 2 + T1 + 1))) * (m1 - T1 + 1) ** (
                                1 - np.exp(-psi * omega)) +
                    c * T1 ** 3 + 3 * a * T1 - 3 * b * p * T1 - 3 * np.exp(psi * omega) * T1 * (
                                5 * a - 5 * b * p + c * T1 * (m1 + T1 + 1)) + np.exp(2 * psi * omega) * T1 * (
                                18 * a - 18 * b * p + c * (6 * (m1 + 1) ** 2 + 3 * T1 * (m1 + 1) + 2 * T1 ** 2))
            ) / (3 * denominator1)

    term7 = lambda_val * (-T*Z + term4 + term5 + term6)

    numerator3 = C1 * np.exp(psi * omega) * (m1 - T1 + 1) ** (-np.exp(-psi * omega)) * (
            -1 / 12 * T1 * (
            6 * b * (1 - 5 * np.exp(psi * omega) + 6 * np.exp(2 * psi * omega)) * p * (2 * m1 - T1 + 2) +
            6 * a * (1 - 5 * np.exp(psi * omega) + 6 * np.exp(2 * psi * omega)) * (-2 * m1 + T1 - 2) +
            c * ((-4 * m1 + 3 * T1 - 4) * T1 ** 2 - np.exp(psi * omega) * (
                -12 * (m1 + 1) ** 2 - 4 * T1 * (m1 + 1) + 9 * T1 ** 2) * T1 - 6 * np.exp(2 * psi * omega) * (
                             4 * (m1 + 1) ** 3 - T1 ** 3))
    ) * (m1 - T1 + 1) ** np.exp(-psi * omega) -
            np.exp(psi * omega) * (
                        (m1 + 1) ** (1 + np.exp(-psi * omega)) - (m1 - T1 + 1) ** (1 + np.exp(-psi * omega))) * (
                        m1 - T1 + 1) * (c * T1 ** 2 + a - b * p - np.exp(psi * omega) * (
                5 * a - 5 * b * p + c * T1 * (2 * m1 + 3 * T1 + 2)) +
                                        2 * np.exp(2 * psi * omega) * (3 * a - 3 * b * p + c * (
                        m1 ** 2 + (T1 + 2) * m1 + T1 ** 2 + T1 + 1))
                                        )) / (1 + np.exp(psi * omega))

    term8 = numerator3 / denominator1

    numerator4 = C4 * ((T - T1) * delta * (6 * (a - b * p) * delta ** 2 + c * (
                delta * (9 * T + 3 * T1 + 2 * (T ** 2 + T1 * T + T1 ** 2) * delta) + 6)) -
                       6 * ((a - b * p) * delta ** 2 + c * (T * delta + 1) ** 2) * np.log(
                T * delta - T1 * delta + 1))

    term9 = numerator4 / (6 * delta ** 3)

    numerator5 = C3 * ((T - T1) * delta * (6 * (a - b * p) * delta ** 2 + c * (
                delta * (9 * T + 3 * T1 + 2 * (T ** 2 + T1 * T + T1 ** 2) * delta) + 6)) -
                       6 * C3 * ((a - b * p) * delta ** 2 + c * (T * delta + 1) ** 2) * np.log(
                T * delta - T1 * delta + 1))

    term10 = numerator5 / (6 * delta ** 4)

    term11 = C2 * (
            -3 * np.exp(psi * omega) * (
                (m1 + 1) ** np.exp(-psi * omega) - (m1 - T1 + 1) ** np.exp(-psi * omega)) * (
                        c * T1 ** 2 + a - b * p - np.exp(psi * omega) * (
                            5 * a - 5 * b * p + c * T1 * (2 * m1 + 3 * T1 + 2)) +
                        2 * np.exp(2 * psi * omega) * (
                                    3 * a - 3 * b * p + c * ((m1 + 1) ** 2 + T1 * (m1 + 1) + T1 ** 2)) * (
                                    m1 - T1 + 1) ** (1 - np.exp(-psi * omega)) + c * T1 ** 3 +
                        3 * a * T1 - 3 * b * p * T1 - 3 * np.exp(psi * omega) * T1 * (
                                    5 * a - 5 * b * p + c * T1 * (m1 + T1 + 1)) + np.exp(2 * psi * omega) * T1 * (
                                    18 * a - 18 * b * p + c * (6 * (m1 + 1) ** 2 + 3 * T1 * (m1 + 1) + 2 * T1 ** 2))
                        ) / (3 * denominator1))

    final_result = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + term11 + T * omega)/T
    return -1000000000 if final_result <0 else -final_result



search_boundaries = [[0.01, 0.01], [3, 3]]
alg = BeesAlgorithm(complex_equation, search_boundaries[0], search_boundaries[1])
alg.performFullOptimisation(max_iteration=5000)
best = alg.best_solution
print(best.score, best.values)

#     # Example usage:
# result = complex_equation(space, A, C0, m1, a, b, p, psi, omega, c, delta, G, k1, k2, k3, k4, xi, lambda_val, C1, C2, C3, C4)


    #PENALTY METHOD
