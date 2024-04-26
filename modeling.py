from pprint import pprint
from ephemeris import *
from main1 import sat_num, sat_arr
from consts import R
import math
from rk4 import *
from consts import *
from for_modeling import *
from for_time import *
import numpy as np
from for_rec_coord import func_rec_coord

R = R * 1e3  # Перевели радиус Земли в м
c = c * 1e3  # Перевели скорость света в м/с

# Для простоты будем считать, что показания часов на приемнике равны показаниям часов МДВ

points_num = 3  # Количество точек траектории
rec_coord = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num + 1)}
rec_speed = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num)}
rec_time = {f'point_{i}': 0 for i in range(1, points_num + 1)}

# 3 времени измерения, [c] (максимум +- 10 мин от tb = 78300 с)
time_step = 5  # Шаг по времени, [с]
rec_time['point_1'] = 78305
rec_time['point_2'] = 78310
rec_time['point_3'] = 78315

# 3 точки траектории приемника в географической системе координат (широта, долгота)
# Средняя  скорость бега человека - 8 км/ч (2,2 м/с), автомобиля - 40 км/ч (11,1 м/с), перемещение не сильно больше 100м
# Здесь будет примерно человек, двигающийся по широте (долгота = const)

rec_coord_geographic = {f'point_{i}': [0 for _ in range(2)] for i in range(1, points_num + 1)}
# Москва:
rec_coord_geographic['point_1'] = [0.95993108859, 0.65839]
rec_coord_geographic['point_2'] = [0.95993308859, 0.65839]
rec_coord_geographic['point_3'] = [0.95993508859, 0.65839]
# Что-то рандомное (первый варинт, вроде, тоже рабочий)
# rec_coord_geographic['point_1'] = [0.972970, 0.65839]
# rec_coord_geographic['point_2'] = [0.972972, 0.65839]
# rec_coord_geographic['point_3'] = [0.972974, 0.65839]

# Переводим в геодезическую систему координат
for point in rec_coord.keys():
    rec_coord[point][0] = R * math.cos(rec_coord_geographic[point][0]) * math.cos(rec_coord_geographic[point][1])
    rec_coord[point][1] = R * math.cos(rec_coord_geographic[point][0]) * math.sin(rec_coord_geographic[point][1])
    rec_coord[point][2] = R * math.sin(rec_coord_geographic[point][0])

# 3 вектора скорости приемника в м/с   УТОЧНИТЬ, КАК СЧИТАТЬ!!  Пока считаем просто как приращение коорд на время
# Только для первых 2-х точек:
for i in range(1, points_num):
    for j in range(3):
        rec_speed[f'point_{i}'][j] = (rec_coord[f'point_{i + 1}'][j] - rec_coord[f'point_{i}'][j]) / time_step


# Определение координат и скоростей спутников на момент предшествия для каждой точки траектории приемника

# Словарик из состояний всех спутников для каждой точки
sat_state = {f'point_{i}': [[0 for i in range(6)] for i in range(sat_num)] for i in range(1, points_num + 1)}

# Показания часов МДВ (и приемника) на момент предшествия спутников для всех точек:
tpr = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}

# Задержка на время распространения сигнала (разница показаний часов МДВ на момент изм и момент предшествия)
tz = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}

# Начальное приближение: tпр = tизм, далее корректируем
# 10 раз будем приближать tпр к "реальному" от начального приближения - tизм

for point in sat_state.keys():
    for i in range(sat_num):
        for _ in range(10):
            tpr_temp = rec_time[point] - tz[point][i]  # Показ. часов МДВ на момент предш. на данном этапе приближения
            if sat_arr[i].tb < tpr_temp:
                h = 10
            else:
                h = -10
            sat_state[point][i] = [1e3 * _ for _ in rk4(sat_arr[i], h, tpr_temp)]
            tz[point][i] = dist(sat_state[point][i][:3], rec_coord[point]) / c


for point in tpr.keys():
    tpr[point] = [(rec_time[point] - tz[point][i]) for i in range(sat_num)]


# По tпр моделируем показания часов спутника на момент предшествия и псевдозадержку

# Показания часов спутников на момент предшествия
T_sat = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}
for point in T_sat.keys():
    for i in range(sat_num):
        T_sat[point][i] = tpr[point][i] - sat_arr[i].tau + sat_arr[i].gamma * (tpr[point][i] - sat_arr[i].tb) - sat_arr[i].tauSys

# Псевдозадержка
pseudo_time = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}
for point in T_sat.keys():
    for i in range(sat_num):
        pseudo_time[point][i] = rec_time[point] - T_sat[point][i]

# Находим координаты приемника по смоделированным псевдозадержкам

# Считаем координаты и скорости спутников по смоделированной псевдозедержке
sat_state_calc = {f'point_{i}': [[0 for i in range(6)] for i in range(sat_num)] for i in range(1, points_num + 1)}

# Разница показаний часов спутника и часов МДВ на момент предшествия
dTsys = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}

for point in sat_state_calc.keys():
    for i in range(sat_num):
        T_sat_calc = func_sat_time(rec_time[point], pseudo_time[point][i])
        tpr_calc = func_mos_time(T_sat_calc, sat_arr[i].tb, sat_arr[i].tau, sat_arr[i].gamma, sat_arr[i].tauSys)
        if sat_arr[i].tb < tpr_calc:
            h = 10
        else:
            h = -10
        sat_state_calc[point][i] = [1e3 * i for i in rk4(sat_arr[i], h, tpr_calc)]
        dTsys[point][i] = tpr_calc - T_sat_calc

# Считаем координаты приемника по смоделированной псевдозадержке,
# точнее - по посчитанным по ней координатам и скоростям спутников, а также dTsys

rec_coord_calc = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num + 1)}

for point in rec_coord_calc.keys():
    rec_coord_calc[point] = func_rec_coord(c, sat_num, pseudo_time[point], dTsys[point], sat_state_calc[point])[0]

# for point in rec_coord.keys():
#     print(point, (rec_coord[point][0] - rec_coord_calc[point][0]) / rec_coord[point][0], (rec_coord[point][1] - rec_coord_calc[point][1]) / rec_coord[point][1], (rec_coord[point][2] - rec_coord_calc[point][2]) / rec_coord[point][2])

# На данном этапе получилось смоделировать псевдозадержки (идеальное измерение)
# и рассчитать координаты приемника (3 точки траектории)
# Далее пробуем смоделировать большее количество точек и зашумленные измерения псевдозадержек

# pprint(rec_coord_calc)
# pprint(rec_coord)
