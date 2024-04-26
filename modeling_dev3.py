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
import statistics as stat
import matplotlib.pyplot as plt
from for_rec_speed import *
import matplotlib.pyplot as plt

R = R * 1e3  # Перевели радиус Земли в м
c = c * 1e3  # Перевели скорость света в м/с

# Для простоты будем считать, что показания часов на приемнике равны показаниям часов МДВ

points_num = 10  # Количество точек траектории
rec_coord = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num + 1)}
rec_speed = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num)}

# points_num времен измерения, [c] (максимум +- 10 мин от tb = 78300 с)
time_step = 5  # Шаг по времени, [с]
rec_time = {f'point_{i}': (78300 + time_step * i) for i in range(1, points_num + 1)}

# points_num точек траектории приемника в географической системе координат (широта, долгота)
# Средняя  скорость бега человека - 8 км/ч (2,2 м/с), автомобиля - 40 км/ч (11,1 м/с), перемещение не сильно больше 100м
# Здесь будет примерно человек, двигающийся по широте (долгота = const)

# Москва:
lat_step = 0.000002  # Шаг по широте [рад]
rec_coord_geographic = {f'point_{i}': [0.95993108859 + lat_step * (i - 1), 0.65839] for i in range(1, points_num + 1)}

# Переводим в геодезическую систему координат
for point in rec_coord.keys():
    rec_coord[point][0] = R * math.cos(rec_coord_geographic[point][0]) * math.cos(rec_coord_geographic[point][1])
    rec_coord[point][1] = R * math.cos(rec_coord_geographic[point][0]) * math.sin(rec_coord_geographic[point][1])
    rec_coord[point][2] = R * math.sin(rec_coord_geographic[point][0])

# 3 вектора скорости приемника в м/с   УТОЧНИТЬ, КАК СЧИТАТЬ!!  Пока считаем просто как приращение коорд на время
# Не считаем скорость в последней точке траектории:
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
        T_sat[point][i] = tpr[point][i] - sat_arr[i].tau + sat_arr[i].gamma * (tpr[point][i] - sat_arr[i].tb) - sat_arr[
            i].tauSys

# Псевдозадержка (без шумов)
pseudo_time = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}
for point in T_sat.keys():
    for i in range(sat_num):
        pseudo_time[point][i] = rec_time[point] - T_sat[point][i]

# Добавляем нормальный шум к псевдозадержкам
noise_num_coord = 50  # Количество различных шумов
rep_num_coord = 50  # Количество расчетов с шумом одной дисперсии

# Дисперсии шумов:
noise_arr_coord = [(i+1)*3e-7 for i in range(noise_num_coord)]
# noise_arr_coord = [0.4e-8, 1e-8]  # можно задавать вручую

pseudo_time_noisy = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}': [
            [
                (
                        pseudo_time[f'point_{_}'][j]
                        + np.random.normal(loc=0.0, scale=noise_arr_coord[i])
                )
                for j in range(sat_num)
            ]
            for t in range(rep_num_coord)
        ]
        for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

# Находим координаты приемника по смоделированным псевдозадержкам

# Считаем координаты и скорости спутников по смоделированной псевдозедержке

# По незашумленным псевдозадержкам:
sat_state_calc = {f'point_{i}': [[0 for i in range(6)] for i in range(sat_num)] for i in range(1, points_num + 1)}
# По зашумленнаям псевдозаржкам:
sat_state_calc_noisy = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}': [
            [
                [0 for i in range(6)]
                for i in range(sat_num)
            ]
            for j in range(rep_num_coord)
        ]
        for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

# Разница показаний часов спутника и часов МДВ на момент предшествия

# По незашумленным псевдозадержкам:
dTsys = {f'point_{i}': [0 for i in range(sat_num)] for i in range(1, points_num + 1)}
# По зашумленным псевдозадержкам:
dTsys_noisy = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}': [
            [0 for i in range(sat_num)] for j in range(rep_num_coord)
        ]
        for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

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

for point in sat_state_calc_noisy.keys():
    for noise in sat_state_calc_noisy[point].keys():
        for j in range(rep_num_coord):
            for i in range(sat_num):
                T_sat_calc = func_sat_time(rec_time[point], pseudo_time_noisy[point][noise][j][i])
                tpr_calc = func_mos_time(T_sat_calc, sat_arr[i].tb, sat_arr[i].tau, sat_arr[i].gamma, sat_arr[i].tauSys)
                if sat_arr[i].tb < tpr_calc:
                    h = 10
                else:
                    h = -10
                sat_state_calc_noisy[point][noise][j][i] = [1e3 * i for i in rk4(sat_arr[i], h, tpr_calc)]
                dTsys_noisy[point][noise][j][i] = tpr_calc - T_sat_calc

# Считаем координаты приемника по смоделированным псевдозадержкам (зашумленным и незашумленным),
# точнее - по посчитанным по ней координатам и скоростям спутников, а также dTsys

# По незашумленным псевдозадержкам:
rec_coord_calc = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num + 1)}

# По зашумленным псевдозадержкам:
rec_coord_calc_noisy = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}': [
            [0 for i in range(3)] for j in range(rep_num_coord)
        ]
        for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

# Считаем матрицы H:
# По незашумленным псевдозадержкам:
H = {f'point_{_}': np.array([[0 for i in range(4)]for j in range(sat_num)]) for _ in range(1, points_num + 1)}

# По зашумленным псевдозадержкам:
# H_noisy = {
#     f'point_{_}': {
#         f'noise_{noise_arr_coord[__]}': [
#             [
#                 [0 for i in range(4)] for j in range(sat_num)
#             ]
#             for ___ in range(rep_num_coord)
#         ]
#         for __ in range(noise_num_coord)
#     }
#     for _ in range(1, points_num + 1)
# }
# H_noisy = np.array(H_noisy)

for point in rec_coord_calc.keys():
    rec_coord_calc[point], H[point] = func_rec_coord(c,
                                                     sat_num,
                                                     pseudo_time[point],
                                                     dTsys[point],
                                                     sat_state_calc[point])

# Для пересчета в местную СК:
rec_coord_calc_noisy_loc = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}': [
            [0 for i in range(3)] for j in range(rep_num_coord)
        ]
        for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

for point in rec_coord_calc_noisy.keys():
    for noise in rec_coord_calc_noisy[point].keys():
        for i in range(rep_num_coord):
            rec_coord_calc_noisy[point][noise][i] = func_rec_coord(
                c,
                sat_num,
                pseudo_time_noisy[point][noise][i],
                dTsys_noisy[point][noise][i],
                sat_state_calc_noisy[point][noise][i]
            )[0]
            rec_coord_calc_noisy_loc[point][noise][i] = get_nhe_coord(rec_coord_calc_noisy[point][noise][i],
                                                                      rec_coord_calc['point_1'])
# y = [0 for i in range(points_num)]
# y0 = [0 for i in range(points_num)]
# for i in range(points_num):
#     y[i] = rec_coord_calc_noisy_loc[f'point_{i+1}'][f'noise_{noise_arr_coord[0]}'][0][0]
#     y0[i] = rec_coord_calc_noisy[f'point_{i+1}'][f'noise_{noise_arr_coord[0]}'][0][0]
# t = [sat_arr[0].tb + time_step*(i + 1) for i in range(points_num)]
# plt.plot(t, y)
# plt.show()
# plt.plot(t, y0)
# plt.show()
# Посчитали rep_num раз координаты приемника для одной точки с шумом (определение псевдозадержек) одной дисперсии
# Теперь найдем дисперсию ошибок определения этих координат в зависимости от дисперсии зашумленных псевдозадержек


# Средение (по шуму с одной дисперсией) значения координат приемника для каждой точки и для каждого шума:
rec_coord_mean = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}':
            [0 for j in range(3)] for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

for point in rec_coord_mean.keys():
    for noise in rec_coord_mean[point].keys():
        for i in range(3):
            temp = np.array(rec_coord_calc_noisy[point][noise]).T[i]
            rec_coord_mean[point][noise][i] = stat.mean(temp)

# Среднеквадратичное отклонение значений координат приемника (population standard deviation of data):
rec_coord_sigma = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}':
            [0 for j in range(3)] for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}
rec_coord_sigma_loc = {
    f'point_{_}': {
        f'noise_{noise_arr_coord[i]}':
            [0 for j in range(3)] for i in range(noise_num_coord)
    }
    for _ in range(1, points_num + 1)
}

for point in rec_coord_mean.keys():
    for noise in rec_coord_mean[point].keys():
        for i in range(3):
            temp = np.array(rec_coord_calc_noisy[point][noise]).T[i]
            rec_coord_sigma[point][noise][i] = stat.pstdev(temp)

            if point == 'point_1':
                temp1 = np.array(rec_coord_calc_noisy_loc['point_1'][noise]).T[i]
                rec_coord_sigma_loc['point_1'][noise][i] = stat.pstdev(temp1)

# Работаем со скоростью (последнюю точку не рассматриваем!!)

# Все числа, которые нужны для работы с псевдодоплеровскими смещениями:
params = {'c': 299792458.0,  # Скорость света
          'f0': 1602 * 1e6,  # Для расчета номинальной частоты спутникова сигнала
          'df': 562.5 * 1e3,  # Для расчета номинальной частоты спутникова сигнала
          'f_mo': 5 * 1e6,  # Номинальная частота задающего генератора
          'df_mo': 5  # Нестабильность задающего генератора
          }
# Моделируем псевдодоплеровские смещения:
# Без зашумления:
pseudo_dop = {
    f'point_{_}': [0 for i in range(sat_num)]
    for _ in range(1, points_num)
}
# Длины волн спутниковых сигналов:
lam = [0 for i in range(sat_num)]
for point in pseudo_dop.keys():
    for i in range(sat_num):
        pseudo_dop[point][i], lam[i] = get_pseudo_dop(rec_coord[point],
                                                      rec_speed[point],
                                                      sat_arr[i],
                                                      sat_state[point][i],
                                                      params)

# С шумом:
noise_num_speed = 50  # Количество различных шумов
rep_num_speed = 50  # Количество расчетов с шумом одной дисперсии

# Дисперсии шумов:
noise_arr_speed = [(i+1)*1e-4 for i in range(noise_num_speed)]
# noise_arr_speed = [5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 1e-1, 1.1e-1, 1.2e-1]  # можно задавать вручую

pseudo_dop_noisy = {
    f'point_{_}': {
        f'noise_{noise_arr_speed[i]}': [
            [
                (
                        pseudo_dop[f'point_{_}'][j]
                        + np.random.normal(loc=0.0, scale=noise_arr_speed[i])
                )
                for j in range(sat_num)
            ]
            for t in range(rep_num_speed)
        ]
        for i in range(noise_num_speed)
    }
    for _ in range(1, points_num)
}

# Считаем по смоделированным псевдодоплерам (незашумл. и зашумл.) скорости приемника

rec_speed_calc = {f'point_{i}': [0 for _ in range(3)] for i in range(1, points_num)}
for point in rec_speed_calc.keys():
    rec_speed_calc[point] = func_rec_speed(sat_num, sat_state[point], H[point], pseudo_dop[point], lam, params)

rec_speed_calc_noisy = {
    f'point_{_}': {
        f'noise_{noise_arr_speed[__]}': [
            [0 for i in range(3)] for j in range(rep_num_speed)
        ]
        for __ in range(noise_num_speed)
    }
    for _ in range(1, points_num)
}

rec_speed_calc_noisy_loc = {
    f'point_{_}': {
        f'noise_{noise_arr_speed[__]}': [
            [0 for i in range(3)] for j in range(rep_num_speed)
        ]
        for __ in range(noise_num_speed)
    }
    for _ in range(1, points_num)
}

for point in rec_speed_calc_noisy.keys():
    for noise in rec_speed_calc_noisy[point].keys():
        for i in range(rep_num_speed):
            rec_speed_calc_noisy[point][noise][i] = func_rec_speed(sat_num,
                                                                   sat_state[point],
                                                                   H[point],
                                                                   pseudo_dop_noisy[point][noise][i],
                                                                   lam,
                                                                   params)

            rec_speed_calc_noisy_loc[point][noise][i] = get_nhe_speed(rec_speed_calc_noisy[point][noise][i],
                                                                      rec_coord_calc['point_1'])


# Средение (по шуму с одной дисперсией) значения компонент скорости приемника для каждой точки и для каждого шума:
rec_speed_mean = {
    f'point_{_}': {
        f'noise_{noise_arr_speed[i]}':
            [0 for j in range(3)] for i in range(noise_num_speed)
    }
    for _ in range(1, points_num)
}

for point in rec_speed_mean.keys():
    for noise in rec_speed_mean[point].keys():
        for i in range(3):
            temp = np.array(rec_speed_calc_noisy[point][noise]).T[i]
            rec_speed_mean[point][noise][i] = stat.mean(temp)

# Среднеквадратичное отклонение значений компонент скорости приемника (population standard deviation of data):
rec_speed_sigma = {
    f'point_{_}': {
        f'noise_{noise_arr_speed[i]}':
            [0 for j in range(3)] for i in range(noise_num_speed)
    }
    for _ in range(1, points_num)
}
rec_speed_sigma_loc = {
    f'point_{_}': {
        f'noise_{noise_arr_speed[i]}':
            [0 for j in range(3)] for i in range(noise_num_speed)
    }
    for _ in range(1, points_num)
}

for point in rec_speed_mean.keys():
    for noise in rec_speed_mean[point].keys():
        for i in range(3):
            temp = np.array(rec_speed_calc_noisy[point][noise]).T[i]
            rec_speed_sigma[point][noise][i] = stat.pstdev(temp)

            if point == 'point_1':
                temp1 = np.array(rec_speed_calc_noisy_loc[point][noise]).T[i]
                rec_speed_sigma_loc['point_1'][noise][i] = stat.pstdev(temp1)
            else:
                pass


# Строим графики


#
# for point in rec_speed_calc.keys():
#     print(point)
#     print((rec_speed_calc[point][0] - rec_speed_calc_noisy[point][f'noise_{noise_arr_speed[0]}'][0][0])*100)
#     print((rec_speed_calc[point][0] - rec_speed_calc_noisy[point][f'noise_{noise_arr_speed[1]}'][0][0])*100)
#
#     print((rec_speed_calc[point][0], rec_speed_calc_noisy[point][f'noise_{noise_arr_speed[0]}'][0][0])*100)
#     print((rec_speed_calc[point][0], rec_speed_calc_noisy[point][f'noise_{noise_arr_speed[1]}'][0][0])*100)

# Строим координатный сравнительный график (компонента x):
# plot_comp(rec_coord_calc, rec_coord_calc_noisy, noise_arr_coord, 0, time_step, sat_arr[0].tb, points_num, 0, True)
# plot_comp(rec_speed_calc, rec_speed_calc_noisy, noise_arr_speed, 0, time_step, sat_arr[0].tb, points_num - 1, 0, False)

# # Строим скоростной сравнительный график (компонента x):
# plot_comp(rec_speed_calc, rec_speed_calc_noisy, noise_arr_speed, 0, time_step, sat_arr[0].tb, points_num - 1, 0, False)
#
# # Строим график оценки сигмы ошибок координат от сигмы шума (геодезическая СК):
# pprint(rec_speed_sigma)
# print(noise_arr_speed)
plot_std_dev(rec_coord_sigma_loc, noise_arr_coord, noise_num_coord, True, False)
plot_std_dev(rec_speed_sigma_loc, noise_arr_speed, noise_num_speed, False, False)

# pprint(rec_coord_sigma)
# plot_std_dev(rec_coord_sigma, noise_arr_coord, noise_num_coord, True, True)
#
# # Строим график оценки сигмы ошибок скорости от сигмы шума (геодезическая СК):
# plot_std_dev(rec_speed_sigma, noise_arr_speed, noise_num_speed, False)


# # Нужно построить графики в местной системе координат
# rec_coord_sigma_loc = {'point_1': rec_coord_sigma_loc}
# rec_speed_sigma_loc = {'point_1': rec_speed_sigma_calc}
# plot_std_dev(rec_coord_sigma_loc, noise_arr_coord, noise_num_coord, True, False)
# plot_std_dev(rec_speed_sigma_loc, noise_arr_speed, noise_num_speed, False, False)

# print('Сравниваем координаты')
# print(rec_coord['point_1'])
# print(rec_coord_calc['point_1'])
# print('Сигма ошибки = ', rec_coord_sigma['point_1'][f'noise_{noise_arr_coord[0]}'])
#
# print('Сравниваем скорость')
# print(rec_speed['point_1'])
# print(rec_speed_calc['point_1'])
# print('Сигма ошибки = ', rec_speed_sigma['point_1'][f'noise_{noise_arr_speed[0]}'])
#
# print('Сравниваем calc и mean координаты')
# print(rec_coord_calc['point_1'])
# print(rec_coord_mean['point_1'][f'noise_{noise_arr_coord[0]}'])
#
# print('Сравниваем calc и mean скорости')
# print(rec_speed_calc['point_1'])
# print(rec_speed_mean['point_1'][f'noise_{noise_arr_speed[0]}'])
#
# print('Сравниваем сптниковые состояния - "реальные " и посчитанные')
# print(rec_speed['point_1'])
# print(rec_speed_calc)


# Еще для спутниковых координат сравнение (calc и оценка методом посл. приближ.)
# Спросить какие вообще еще могут быть ошибки и как мы их НЕ учичтываем
# Спросить про размерности псевдозадержки и псевдодоплеровского смещения
