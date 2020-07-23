from scipy.optimize import fsolve
from math import tan, tanh, sqrt, pi


def mode_constants(vars):
    k, kappa, omega = vars
    distensibility, length, period_s, rho = 1.29e-4, 10, 6, 1030
    eq1 = (k * length / 2) * tanh(kappa * length / 2) - (kappa * length / 2) * tan(k * length / 2)
    eq2 = k ** 2 - ((2 * pi) / (distensibility * period_s)) * (
                sqrt(1 + ((period_s / pi) * rho * distensibility ** 2 * omega ** 2)) - 1)
    eq3 = kappa ** 2 - ((2 * pi) / (distensibility * period_s)) * (
                sqrt(1 + ((period_s / pi) * rho * distensibility ** 2 * omega ** 2)) + 1)
    return [eq1, eq2, eq3]

k, kappa, omega = fsolve(mode_constants, (5, 5, 5))
print(k, kappa, omega)