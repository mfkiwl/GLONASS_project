import math
from consts import *
from ephemeris import Ephemeris


def func(arg, sat: Ephemeris):
    f = [0 for i in range(6)]
    f[0] = arg[3]
    f[1] = arg[4]
    f[2] = arg[5]

    r = math.sqrt(arg[0] ** 2 + arg[1] ** 2 + arg[2] ** 2)
    A = mu / r ** 3

    f[3] = (omega_e ** 2 - A) * arg[0] + 2 * omega_e * arg[4] + 3 / 2 * C20 * mu * R ** 2 / r ** 5 * arg[0] * (
            1 - 5 * arg[2] ** 2 / r ** 2) + sat.w[0]

    f[4] = (omega_e ** 2 - A) * arg[1] - 2 * omega_e * arg[3] + 3 / 2 * C20 * mu * R ** 2 / r ** 5 * arg[1] * (
            1 - 5 * arg[2] ** 2 / r ** 2) + sat.w[1]

    f[5] = -A * arg[2] + 3 / 2 * C20 * mu * R ** 2 / r ** 5 * arg[2] * (3 - 5 * arg[2] ** 2 / r ** 2) + sat.w[2]

    return f


def rk4(sat: Ephemeris, h, T):
    s = [0 for i in range(6)]  # массив "состояния" спутника
    for i in range(3):
        s[i] = sat.r[i]
        s[i + 3] = sat.v[i]
    time = sat.tb
    flag = True

    while flag:
        if h > 0:
            if time + h < T:
                time += h
            else:
                flag = False
                h = T - time
        else:
            if time + h > T:
                time += h
            else:
                flag = False
                h = T - time

        arg = [i for i in s]
        k1 = [h * i for i in func(arg, sat)]

        arg = [i + 1 / 2 * j for i, j in zip(s, k1)]
        k2 = [h * i for i in func(arg, sat)]

        arg = [i + 1 / 2 * j for i, j in zip(s, k2)]
        k3 = [h * i for i in func(arg, sat)]

        arg = [i + j for i, j in zip(s, k3)]
        k4 = [h * i for i in func(arg, sat)]

        s = [si + (k1i + 2 * k2i + 2 * k3i + k4i) / 6 for si, k1i, k2i, k3i, k4i in zip(s, k1, k2, k3, k4)]

    return s