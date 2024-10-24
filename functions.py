import math
import numpy as np
from ephemeris import Ephemeris
from consts import *


def for_rk4(arg: list[float], sat: Ephemeris) -> list[float]:
    """
    Вспомогательная функция для функции rk4
    """

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


def rk4(sat: Ephemeris, h: int, T) -> list[float]:
    """
    Рассчитывает состояние спутника (координаты и компоненты вектора скорости) 
    с помощью численного интегрирования методом Рунге-Кутты 4-го порядка

    :param sat: Спутник (объект класса Ephemeris)
    :type sat: Ephemeris
    :param h: Шаг интегрирования 
    :type h: int
    :param T: Показания часов МДВ на момент предшествия  
    :type T: float
    
    :return: Состояние спутника - его координаты и компоненты вектора скорости 
    :rtype: list[float] 
    """

    s = [0 for i in range(6)]  # Массив "состояния" спутника
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
        k1 = [h * i for i in for_rk4(arg, sat)]

        arg = [i + 1 / 2 * j for i, j in zip(s, k1)]
        k2 = [h * i for i in for_rk4(arg, sat)]

        arg = [i + 1 / 2 * j for i, j in zip(s, k2)]
        k3 = [h * i for i in for_rk4(arg, sat)]

        arg = [i + j for i, j in zip(s, k3)]
        k4 = [h * i for i in for_rk4(arg, sat)]

        s = [si + (k1i + 2 * k2i + 2 * k3i + k4i) / 6 for si, k1i, k2i, k3i, k4i in zip(s, k1, k2, k3, k4)]

    return s


def dist(a: list[float], b: list[float]) -> float:
    """
    Рассчитывает расстояние между двумя точками

    :param a: Первая точка
    :type a: list[float]
    :param b: Вторая точка 
    :type b: list[float]
    
    :return: Расстояние между двумя точками
    :rtype: float 
    """

    res = ((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)**0.5
    return res


def func_sat_time(rec_time, pseudo_time):
    """
    Рассчитывает показания часов спутника на момент предшествия

    :param rec_time: Показания часов приемника на момент измерения 
    :type rec_time: float
    :param pseudo_time: Псевдозадержка 
    :type pseudo_time: float
    
    :return: Показания часов спутника на момент предшествия
    :rtype: float 
    """

    res = rec_time - pseudo_time
    return res


def func_mos_time(sat_time, tb, tau, gamma, tauSys):
    """
    Рассчитывает показания часов МДВ на момент предшествия

    :param sat_time: Показания часов спутника на момент предшествия 
    :type sat_time: float
    :param tb: Эфемеридная информация  
    :type tb: float
    :param tau: Эфемеридная информация  
    :type tau: float
    :param gamma: Эфемеридная информация  
    :type gamma: float
    :param tauSys: Эфемеридная информация  
    :type tauSys: float
    
    :return: Показания часов МДВ на момент предшествия
    :rtype: float 
    """

    res = sat_time + tau - gamma * (sat_time - tb) + tauSys
    while res < 0:
        res += 86400
    return res


def func_rec_coord(c: float,
                   sat_num: int,
                   pseudo_time: list[float],
                   dTsys: list[float],
                   sat_state: list[float], eps=1e-3) -> tuple[list[float], list[list[float]]]:
    """
    Расcчитывает координаты приемника и матрицу направляющих косинусов H 

    :param c: Скорость света в м/с
    :type c: float
    :param sat_num: Количество спутников
    :type sat_num: int
    :param pseudo_time: Псевдозадержки для каждого спутника
    :type pseudo_time: list[float]
    :param dTsys: Разница показаний часов спутников и часов МДВ на моменты предшествия
    :type dTsys: list[float]
    :param sat_state: "Состояния" всех спутников
    :type sat_state: list[float]
    :param eps: Ошибка
    :type eps: float
    :return: Координаты приемника ([0]) и матрицу H ([1])
    :rtype gamma: tuple[list[float], list[list[float]]]
    """
    
    cond = [0.0 for i in range(4)]  # x, y, z, dD - координаты и смещение показаний часов приемника
    hx = [0.0 for i in range(sat_num)]
    hy = [0.0 for i in range(sat_num)]
    hz = [0.0 for i in range(sat_num)]
    D = [0.0 for i in range(sat_num)]
    xi = [0.0 for i in range(sat_num)]
    H = [[1 for _ in range(4)] for i in range(sat_num)]
    for _ in range(100):
        for i in range(sat_num):
            D[i] = dist(cond, sat_state[i])
            hx[i] = (cond[0] - sat_state[i][0]) / D[i]
            hy[i] = (cond[1] - sat_state[i][1]) / D[i]
            hz[i] = (cond[2] - sat_state[i][2]) / D[i]
            xi[i] = c * pseudo_time[i] - c * dTsys[i] - D[i] - cond[3]
            H[i][0] = hx[i]
            H[i][1] = hy[i]
            H[i][2] = hz[i]

        H = np.array(H)
        xi = np.array(xi)
        HT = H.T
        HTdotH = HT.dot(H)
        HTdotH_inv = np.linalg.inv(HTdotH)
        dtheta = HTdotH_inv.dot(HT)
        dtheta = dtheta.dot(xi)

        for i in range(4):
            cond[i] += dtheta[i]

        if abs(np.max(dtheta)) < eps:
            break

    # Учитываем вращение Земли

    # Времена распространения сигналов от всех спутников:
    tau = [0.0 for i in range(sat_num)]
    for i in range(sat_num):
        tau[i] = dist(cond, sat_state[i]) / c

    # Вносим поправки в определение координат спутников:

    for i in range(sat_num):
        alpha = omega_e * tau[i]
        xTemp = sat_state[i][0]
        yTemp = sat_state[i][1]
        sat_state[i][0] = xTemp + yTemp * alpha
        sat_state[i][1] = -xTemp * alpha + yTemp

    # Пересчитываем координаты приемника:
    for i in range(sat_num):
        D[i] = dist(cond, sat_state[i])
        hx[i] = (cond[0] - sat_state[i][0]) / D[i]
        hy[i] = (cond[1] - sat_state[i][1]) / D[i]
        hz[i] = (cond[2] - sat_state[i][2]) / D[i]
        xi[i] = c * pseudo_time[i] - c * dTsys[i] - D[i] - cond[3]
        H[i][0] = hx[i]
        H[i][1] = hy[i]
        H[i][2] = hz[i]

    H = np.array(H)
    xi = np.array(xi)
    HT = H.T
    HTdotH = HT.dot(H)
    HTdotH_inv = np.linalg.inv(HTdotH)
    dtheta = HTdotH_inv.dot(HT)
    dtheta = dtheta.dot(xi)

    for i in range(4):
        cond[i] += dtheta[i]

    rec_coord = [i for i in cond[:3]]
    return rec_coord, H


def get_nhe_coord(coord: list[float], coord0: list[float]) -> list[float]:
    """
    Переводит координаты точки из геодезической в местную систему координат

    :param coord: Координаты точки, которые будем переводить в местную СК
    :type coord: list[float]
    :param coord0: Координаты начальной точки
    :type coord0: list[float]
    :return: Координаты точки в новой СК
    :rtype dTsys: list[float]
    """

    x, y, z = coord[0], coord[1], coord[2]
    x0, y0, z0 = coord0[0], coord0[1], coord0[2]

    N = (-x0*z0*(x - x0) - y0*z0*(y - y0) + (x0**2 + y0**2)*(z - z0)) / ((x0*z0)**2 + (y0*z0)**2 + (x0**2 + y0**2)**2)**0.5
    H = (x0*(x - x0) + y0*(y - y0) + z0*(z - z0)) / (x0**2 + y0**2 + z0**2)**0.5
    E = (-y0*(x - x0) + x0*(y - y0)) / (x0**2 + y0**2)**0.5

    return [N, H, E]


def get_pseudo_dop(rec_coord: list[float],
                       rec_speed: list[float],
                       sat: Ephemeris,
                       sat_state: list[float],
                       params: dict[str: float]) -> tuple[float, float]:
    """
    Рассчитывает пседодоплеровское смещение и длину волны спутникого сигнала

    :param rec_coord: Координаты приемника
    :type rec_coord: list[float]
    :param rec_speed: Скорость приемника
    :type rec_speed: list[float]
    :param sat: Спутник
    :type sat: Ephemeris
    :param sat_state: Состояние спутника
    :type sat_state: list[float]
    :param params: Константы 
    :return: Пседодоплеровское смещение и длину волны спутникого сигнала
    :rtype: tuple[float, float]
    """

    r = [sat_state[i] - rec_coord[i] for i in range(3)]  # Вектор, направленный от приемника к спутнику
    mod_r = (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 0.5
    D = sum(r[i] * (rec_speed[i] - sat_state[i + 3]) for i in range(3)) / mod_r

    sat_freq = params['f0'] + sat.frqNum * params['df']
    lam = params['c'] / sat_freq
    k = sat_freq / params['f_mo']

    pseudo_dop = - D / lam - k * params['df_mo']
    return pseudo_dop, lam


def func_rec_speed(sat_num: int, 
                   sat_state: list[list[float]], 
                   H: list[list[float]], 
                   pseudo_dop: list[float], 
                   lam: list[float], 
                   params: dict[str: float]) -> list[float]:
    """
    Рассчитывает компоненты вектора скорости приемника

    :param sat_num: Количество спутников 
    :type sat_num: int
    :param sat_state: Состояния спутников 
    :type sat_state: list[list[float]]
    :param H: Матрица направляющих косинусов 
    :type H: list[list[float]]
    :param pseudo_dop: Псевдодоплеровское смещение для каждого спутника
    :type pseudo_dop: list[float]
    :param lam: Длины волн спутниковых сигналов 
    :type lam: list[float]
    :param params: Константы 
    :type params: dict[str: float]
    :return: Компоненты вектора скорости приемника
    :rtype: list[float]
    """

    lam_mo = params['c'] / params['f_mo']  # Длина волны несущего колебания задающего генератора
    for i in range(sat_num):
        H[i][3] = -lam_mo

    hx, hy, hz, _ = H.T

    ksi = [0 for i in range(sat_num)]
    for i in range(sat_num):
        ksi[i] = lam[i] * pseudo_dop[i] + hx[i] * sat_state[i][3] + hy[i] * sat_state[i][4] + hz[i] * sat_state[i][5]
    ksi = np.array(ksi)

    HT = H.T
    HTdotH = HT.dot(H)
    HTdotH_inv = np.linalg.inv(HTdotH)
    res = HTdotH_inv.dot(HT)
    rec_speed = res.dot(ksi)
    return rec_speed[:3]


def get_nhe_speed(speed: list[float], coord0: list[float]) -> list[float]:
    """
    Переводит вектор скорости из геодезической в местную систему координат

    :param speed: Вектор скорости, который будем переводить в местную СК
    :param coord0: Координаты начальной точки
    :return: Вектор скорости в местной СК
    """

    point = [i + j for i, j in zip(coord0, speed)]
    res = get_nhe_coord(point, coord0)
    return res

