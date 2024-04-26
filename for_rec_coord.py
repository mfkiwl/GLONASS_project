from consts import c, omega_e
from for_modeling import dist
import numpy as np
from main1 import sat_num, pseudo_time, dTsys, sat_state


def func_rec_coord(c: float,
                   sat_num: int,
                   pseudo_time: list[float],
                   dTsys: list[float],
                   sat_state: list[float], eps=1e-3):
    """
    :param c: Скорость света в м/с
    :param sat_num: Количество спутников
    :param pseudo_time: Псевдозадержки для каждого спутника
    :param dTsys: Разница показаний часов спутников и часов МДВ на моменты предшествия
    :param sat_state: "Состояния" всех спутников
    :param eps: Ошибка
    :return: Координаты приемника ([0]) и матрицу H ([1])
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


