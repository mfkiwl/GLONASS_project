import math
import numpy as np
from consts import *
from ephemeris import *
from for_time import *
from rk4 import *


sat_num = 9     # количество спутников
sat_arr = set_sat()  # массив спутнков
rec_time = 79185.000  # Показания часов приемника на момент измерения, [с]
pseudo_time = [0.073360802463432, 0.063914220797371, 0.069785759291336, 0.076138991917530, 0.078236013447872,
               0.077133859571741, 0.065296928503503, 0.065833067985743, 0.077911542130084]  # Псевдозадержки, [с]
sat_state = [[0] * 6 for i in range(sat_num)]  # двумерный массив с "состоянием" каждого спутника
dTsys = [0.0 for i in range(sat_num)]  # разница показаний часов спутника и часов МДВ на момент предшествия

# Координаты и скорости спутников на момент предшествия

for i in range(sat_num):
    temp1 = func_sat_time(rec_time, pseudo_time[i])
    # Показания часов МДВ на момент предшествия:
    temp = func_mos_time(temp1, sat_arr[i].tb, sat_arr[i].tau, sat_arr[i].gamma, sat_arr[i].tauSys)
    # print('t_пр =', temp)
    if sat_arr[i].tb < temp:
        h = 10
    else:
        h = -10
    res = [1e3 * i for i in rk4(sat_arr[i], h, temp)]
    sat_state[i] = res
    dTsys[i] = temp - temp1

# координаты приемника

sat_coord = [[0] * 3 for i in range(9)]
for i in range(9):
    for j in range(3):
        sat_coord[i][j] = sat_state[i][j]

sat_coord = [[1e-3 * i for i in j] for j in sat_coord]  # перевели координаты спутников в км
cond = [0.0 for i in range(4)]  # x, y, z, dD - координаты и смещение показаний часов приемника
hx = [0.0 for i in range(9)]
hy = [0.0 for i in range(9)]
hz = [0.0 for i in range(9)]
D = [0.0 for i in range(9)]
xi = [0.0 for i in range(9)]
H = [[1, 1, 1, 1] for i in range(9)]
eps = 1e-3  # ошибка
for counter in range(100):
    for j in range(9):
        D[j] = math.sqrt(
            (cond[0] - sat_coord[j][0]) ** 2 + (cond[1] - sat_coord[j][1]) ** 2 + (cond[2] - sat_coord[j][2]) ** 2)
        hx[j] = (cond[0] - sat_coord[j][0]) / D[j]
        hy[j] = (cond[1] - sat_coord[j][1]) / D[j]
        hz[j] = (cond[2] - sat_coord[j][2]) / D[j]
        xi[j] = c * pseudo_time[j] - c * dTsys[j] - D[j] - cond[3]
        H[j][0] = hx[j]
        H[j][1] = hy[j]
        H[j][2] = hz[j]

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
tau = [0.0 for i in range(9)]
for j in range(9):
    tau[j] = math.sqrt(
        (cond[0] - sat_coord[j][0]) ** 2 + (cond[1] - sat_coord[j][1]) ** 2 + (cond[2] - sat_coord[j][2]) ** 2) / c

# Вносим поправки в определение координат спутников:

for j in range(9):
    alpha = omega_e * tau[j]
    xTemp = sat_coord[j][0]
    yTemp = sat_coord[j][1]
    sat_coord[j][0] = xTemp + yTemp * alpha
    sat_coord[j][1] = -xTemp * alpha + yTemp

# Пересчитываем координаты приемника:
for j in range(9):
    D[j] = math.sqrt(
        (cond[0] - sat_coord[j][0]) ** 2 + (cond[1] - sat_coord[j][1]) ** 2 + (cond[2] - sat_coord[j][2]) ** 2)
    hx[j] = (cond[0] - sat_coord[j][0]) / D[j]
    hy[j] = (cond[1] - sat_coord[j][1]) / D[j]
    hz[j] = (cond[2] - sat_coord[j][2]) / D[j]
    xi[j] = c * pseudo_time[j] - c * dTsys[j] - D[j] - cond[3]
    H[j][0] = hx[j]
    H[j][1] = hy[j]
    H[j][2] = hz[j]

H = np.array(H)
xi = np.array(xi)
HT = H.T
HTdotH = HT.dot(H)
HTdotH_inv = np.linalg.inv(HTdotH)
dtheta = HTdotH_inv.dot(HT)
dtheta = dtheta.dot(xi)

for i in range(4):
    cond[i] += dtheta[i]

rec_coord = [1e3 * i for i in cond[:3]]
# print('r_пр =', rec_coord)
# print(sat_state)

