import numpy as np
from main1 import sat_num, sat_arr, sat_state, rec_coord, H, hx, hy, hz
from for_rec_speed import *
from consts import *


c = c * 1e3  # Перевели скорость света в м/с
# Моделирование псевдодоплеровских смещений частот

f_mo = 5*1e6  # Master oscillator - задающий генератор, тут это его номинальная частота
df_mo = 5  # Нестабильность задающего генератора
pseudo_dop = [0 for i in range(sat_num)]  # Псевдодоплеровские смещения несущих частот спутников
lam = [0 for i in range(sat_num)]  # Длины волн несущих колебаний спутниковых сигналов
for i in range(sat_num):
    D = dist_change(rec_coord, [0, 0, 0], sat_state[i])
    sat_freq = get_freq(sat_arr[i])
    lam[i] = c / sat_freq
    k = sat_freq / f_mo  # Коэффициент - отношение частоты сутникого сигнала к частоте задающего генератора
    pseudo_dop[i] = get_pseudo_dop_not_mod(D, lam[i], k, df_mo)

# print(pseudo_dop)
params = {'c': 299792458.0,  # Скорость света
          'f0': 1602 * 1e6,  # Для расчета номинальной частоты спутникова сигнала
          'df': 562.5 * 1e3,  # Для расчета номинальной частоты спутникова сигнала
          'f_mo': 5 * 1e6,  # Номинальная частота задающего генератора
          'df_mo': 5  # Нестабильность задающего генератора
          }
# for i in range(sat_num):
#     print(get_pseudo_dop(rec_coord, [0, 0, 0], sat_arr[i], sat_state[i], params)[0])


# Расчет скорости приемника по псевдодоплеровским смещениям
lam_mo = c / f_mo  # Длина волны несущего колебания задающего генератора
for i in range(sat_num):
    H[i][3] = -lam_mo

ksi = [0 for i in range(sat_num)]
for i in range(sat_num):
    ksi[i] = lam[i] * pseudo_dop[i] + hx[i] * sat_state[i][3] + hy[i] * sat_state[i][4] + hz[i] * sat_state[i][5]
ksi = np.array(ksi)

HT = H.T
HTdotH = HT.dot(H)
HTdotH_inv = np.linalg.inv(HTdotH)
res = HTdotH_inv.dot(HT)
res = res.dot(ksi)
print('Скорость приемника =', res[:3])

# print(func_rec_speed(sat_num, sat_state, H, pseudo_dop, lam, params))
