import numpy as np
import statistics as stat
from class_Model import Model
from ephemeris import sat_num, sat_arr
from consts import *
from functions import *


class Point:
    def __init__(self, num: int, time_step: int, lat_step: float) -> None:  
        self.num = num
        self.rec_time = self.init_rec_time(time_step)
        self.state = Model(num, lat_step)   
        self.sat = [[0 for _ in range(6)] for _ in range(sat_num)] 
        self.set_sat_state()
        self.set_tpr()

    def init_rec_time(self, time_step: int) -> int:
        """
        Рассчитывает показания часов приемника на момент измерения 

        :param time_step: Шаг моделирования по времени 
        :type time_step: int
        :return: показания часов приемника на момент измерения, [c]
        :rtype: int
        """
        return 78300 + time_step * self.num

    def set_sat_state(self) -> tuple[list[float], list[list[float]]]:
        """
        Рассчитывает состояния спутников (их координаты и компоненты вектора скорости) на моменты предшествия 
        и задержки на время распространения сигналов со спутников до приемника

        :return: Кортеж, содержащий:
            Список состояний каждого спутника (координаты и компоненты вектора скорости)
            Список задержек сигналов для каждого спутника
        :rtype: tuple[list[float], list[list[float]]]
        """
        tz: list[float] = [0 for _ in range(sat_num)]
        sat_state: list[list[float]] = [[0 for i in range(6)] for i in range(sat_num)]
        for i in range(sat_num):
            for _ in range(10):
                tpr_temp = (
                    self.rec_time - tz[i]
                )  # Показания часов МДВ на момент предшествия на данном этапе приближения
                if sat_arr[i].tb < tpr_temp:
                    h = 10
                else:
                    h = -10
                sat_state[i] = [1e3 * _ for _ in rk4(sat_arr[i], h, tpr_temp)]
                tz[i] = dist(sat_state[i][:3], self.state.coord) / c

        self.sat_state = sat_state
        self.tz = tz
        return sat_state, tz

    def set_tpr(self) -> list[float]:
        """
        Вычисляет показания часов МДВ (и приемника) на моменты предшествия для каждого спутника

        :return: Список показаний часов МДВ (и приемника) на моменты предшествия
        :rtype: list[float]
        """
        self.tpr = [self.rec_time - self.tz[i] for i in range(sat_num)]
        return self.tpr

    def _T_sat(self) -> list[float]:
        """
        Вычисляет показания часов спутников на моменты предшествия

        :return: Список показаний часов спутников на моменты предшествия
        :rtype: list[float]
        """
        T_sat = [0 for i in range(sat_num)]
        for i in range(sat_num):
            T_sat[i] = (
                self.tpr[i]
                - sat_arr[i].tau
                + sat_arr[i].gamma * (self.tpr[i] - sat_arr[i].tb)
                - sat_arr[i].tauSys
            )
        self.T_sat = T_sat

        return self.T_sat

    def _pseudo_time(self) -> list[float]:
        """
        Вычисляет истинные значения псевдозадержек для каждого спутника

        :return: Список псевдозадержек для каждого спутника
        :rtype: list[float]
        """
        pseudo_time = [0 for i in range(sat_num)]
        for i in range(sat_num):
            pseudo_time[i] = self.rec_time - self.T_sat[i]
        self.pseudo_time = pseudo_time
        return self.pseudo_time

    def _pseudo_time_noisy(self, noise_num_coord, rep_num_coord, noise_arr_coord):
        """
        Вычисляет зашумленные значения псевдозадержек для каждого спутника

        :param noise_num_coord: Количество различных шумов  
        :type noise_num_coord: int
        :param rep_num_coord: Количество повторений расчетов с одним шумом 
        :type rep_num_coord: int
        :param noise_arr_coord: Дисперсии шумов  
        :type noise_arr_coord: list[float]
        :return: Трехмерный массив зашумленных значений псевдозадержек \n
            - Первый уровень: различные дисперсии шума \n
            - Второй уровень: повторения для каждого шума \n
            - Третий уровень: зашумленные значения псевдозадержек для каждого спутника \n
        :rtype: list[list[list[float]]]
        """
        self.pseudo_time_noisy = []
        for i in range(noise_num_coord):
            rep_arr = []
            for _ in range(rep_num_coord):
                sat_noise = []
                for j in range(sat_num):
                    sat_noise.append(
                        self.pseudo_time[j]
                        + np.random.normal(loc=0.0, scale=noise_arr_coord[i])
                    )
                rep_arr.append(sat_noise)
            self.pseudo_time_noisy.append(rep_arr)
        return self.pseudo_time_noisy

    def _sat_state_calc(self) -> tuple[list[list[float]], list[int]]:
        """
        Рассчитывает состояния спутников (их координаты и компоненты вектора скорости) 
        и разницу показаний часов спутника и часов МДВ на моменты предшествия по истинным значениям псевдозадержек 
        
        :return: Кортеж, содержащий:\n
            - Состояние каждого спутника на момент предшествия \n
            - Разница показаний часов спутника и часов МДВ на момент предшествия \n
        :rtype: tuple[list[list[float]], list[int]]
        """
        sat_state_calc = [[0 for _ in range(6)] for _ in range(sat_num)]
        dTsys = [0 for i in range(sat_num)]
        for i in range(sat_num):
            T_sat_calc = func_sat_time(self.rec_time, self.pseudo_time[i])
            tpr_calc = func_mos_time(
                T_sat_calc,
                sat_arr[i].tb,
                sat_arr[i].tau,
                sat_arr[i].gamma,
                sat_arr[i].tauSys,
            )
            if sat_arr[i].tb < tpr_calc:
                h = 10
            else:
                h = -10
            sat_state_calc[i] = [1e3 * i for i in rk4(sat_arr[i], h, tpr_calc)]
            dTsys[i] = tpr_calc - T_sat_calc
        self.sat_state_calc = sat_state_calc
        self.dTsys = dTsys
        return self.sat_state_calc, self.dTsys

    def _sat_state_calc_noisy(
        self, noise_num_coord, rep_num_coord, noise_arr_coord
    ) -> None:
        """
        Рассчитывает состояния спутников (их координаты и компоненты вектора скорости) 
        и разницу показаний часов спутника и часов МДВ на моменты предшествия по зашумленным значениям псевдозадержек

        :param noise_num_coord: Количество различных шумов  
        :type noise_num_coord: int
        :param rep_num_coord: Количество повторений расчетов с одним шумом 
        :type rep_num_coord: int
        :param noise_arr_coord: Дисперсии шумов  
        :type noise_arr_coord: list[float]
        """
        sat_state_calc_noisy = {
            f"noise_{noise_arr_coord[_]}": [
                [[0 for _ in range(6)] for _ in range(sat_num)]
                for _ in range(rep_num_coord)
            ]
            for _ in range(noise_num_coord)
        }
        dTsys_noisy = {
            f"noise_{noise_arr_coord[i]}": [
                [0 for i in range(sat_num)] for j in range(rep_num_coord)
            ]
            for i in range(noise_num_coord)
        }
        for noise, noise_i in zip(
            sat_state_calc_noisy.keys(), range(len(sat_state_calc_noisy.keys()))
        ):
            for j in range(rep_num_coord):
                for i in range(sat_num):
                    T_sat_calc = func_sat_time(
                        self.rec_time, self.pseudo_time_noisy[noise_i][j][i]
                    )
                    tpr_calc = func_mos_time(
                        T_sat_calc,
                        sat_arr[i].tb,
                        sat_arr[i].tau,
                        sat_arr[i].gamma,
                        sat_arr[i].tauSys,
                    )
                    if sat_arr[i].tb < tpr_calc:
                        h = 10
                    else:
                        h = -10
                    sat_state_calc_noisy[noise][j][i] = [
                        1e3 * i for i in rk4(sat_arr[i], h, tpr_calc)
                    ]
                    dTsys_noisy[noise][j][i] = tpr_calc - T_sat_calc
        self.sat_state_calc_noisy = sat_state_calc_noisy
        self.dTsys_noisy = dTsys_noisy

    def _rec_coord_calc(self) -> None:
        """
        Вычисляет координаты приемника и матрицу направляющих косинусов по истинным значениям псевдозадержек 
        """

        self.rec_coord_calc, self.H = func_rec_coord(
            c, sat_num, self.pseudo_time, self.dTsys, self.sat_state_calc
        )

    def _rec_coord_noisy(self, rep_num_coord, noise_arr_coord, noise_num_coord) -> None:
        """
        Вычисляет координаты приемника в геодезической и местной системах координат 
        по зашумленным значениям псевдозадержек

        :param rep_num_coord: Количество повторений расчетов с одним шумом 
        :type rep_num_coord: int
        :param noise_arr_coord: Дисперсии шумов  
        :type noise_arr_coord: list[float]
        :param noise_num_coord: Количество различных шумов  
        :type noise_num_coord: int
        """
        rec_coord_calc_noisy = {
            f"noise_{noise_arr_coord[i]}": [
                [0 for _ in range(3)] for _ in range(rep_num_coord)
            ]
            for i in range(noise_num_coord)
        }
        rec_coord_calc_noisy_loc = {
            f"noise_{noise_arr_coord[i]}": [
                [0 for i in range(3)] for j in range(rep_num_coord)
            ]
            for i in range(noise_num_coord)
        }
        for noise_i, noise in enumerate(list(rec_coord_calc_noisy.keys())):
            for i in range(rep_num_coord):
                rec_coord_calc_noisy[noise][i] = func_rec_coord(
                    c,
                    sat_num,
                    self.pseudo_time_noisy[noise_i][i],
                    self.dTsys_noisy[noise][i],
                    self.sat_state_calc_noisy[noise][i],
                )[0]
                if self.num == 1:
                    rec_coord_calc_noisy_loc[noise][i] = get_nhe_coord(
                        rec_coord_calc_noisy[noise][i], self.rec_coord_calc
                    )

        self.rec_coord_calc_noisy = rec_coord_calc_noisy
        self.rec_coord_calc_noisy_loc = rec_coord_calc_noisy_loc

    def _rec_coord_mean(self, noise_arr_coord, noise_num_coord) -> None:
        """
        Вычисляет среднее значение координат приемника для различных шумов  

        :param noise_arr_coord: Дисперсии шумов  
        :type noise_arr_coord: list[float]
        :param noise_num_coord: Количество различных шумов  
        :type noise_num_coord: int
        """
        rec_coord_mean = {
            f"noise_{noise_arr_coord[i]}": [0 for j in range(3)]
            for i in range(noise_num_coord)
        }

        for noise in rec_coord_mean.keys():
            for i in range(3):
                temp = np.array(self.rec_coord_calc_noisy[noise]).T[i]
                rec_coord_mean[noise][i] = stat.mean(temp)
        self.rec_coord_mean = rec_coord_mean

    def _pseudo_dop(self, params, points_num: int) -> None:
        """
        Вычисляет истинные значения псевдодоплеровского смещения для каждого спутника
        и длины волн спутниковых сигналов  

        :param params: Константы   
        :param points_num: Количество точек траектории  
        :type noise_num_coord: int
        """

        pseudo_dop = [0 for i in range(sat_num)]
        # Длины волн спутниковых сигналов:
        lam = [0 for i in range(sat_num)]
        if self.num != points_num:
            for i in range(sat_num):
                pseudo_dop[i], lam[i] = get_pseudo_dop(
                    self.state.coord,
                    self.state.speed,
                    sat_arr[i],
                    self.sat_state[i],
                    params,
                )
        self.pseudo_dop, self.lam = pseudo_dop, lam

    def _pseudo_dop_noisy(self, noise_num_speed: int, noise_arr_speed: list[float], rep_num_speed: int) -> None:
        """
        Добавляет шум к истинным значениям псевдодоплеровского смещения   

        :param noise_num_speed: Количество различных шумов  
        :type noise_num_speed: int
        :param noise_arr_speed: Дисперсии шумов  
        :type noise_arr_speed: list[float]
        :param rep_num_speed: Количество повторений расчетов с одним шумом  
        :type rep_num_speed: list[float]        
        """

        pseudo_dop_noisy = {
            f"noise_{noise_arr_speed[i]}": [
                [
                    (
                        self.pseudo_dop[j]
                        + np.random.normal(loc=0.0, scale=noise_arr_speed[i])
                    )
                    for j in range(sat_num)
                ]
                for t in range(rep_num_speed)
            ]
            for i in range(noise_num_speed)
        }
        self.pseudo_dop_noisy = pseudo_dop_noisy

    def _rec_speed_calc(self, params) -> None:
        """
        Рассчитывает компоненты вектора скорости приемника по истинным значениям псевдодоплеровского смещения    

        :param params: Константы   
        """

        self.rec_speed_calc = func_rec_speed(
            sat_num, self.sat_state, self.H, self.pseudo_dop, self.lam, params
        )

    def _rec_speed_calc_noisy(self, 
                              noise_arr_speed: list[float], 
                              rep_num_speed: int, 
                              noise_num_speed: int, 
                              params, 
                              points_num: int) -> None:
        """
        Рассчитывает компоненты вектора скорости приемника в геодезической и местной системах координат
        по зашумленным значениям псевдодоплеровского смещения    

        :param noise_arr_speed: Дисперсии различных шумов  
        :type noise_arr_speed: list[float]
        :param rep_num_speed: Количество повторений расчетов с одним шумом   
        :type rep_num_speed: int
        :param noise_num_speed: Количество различных шумов  
        :type noise_num_speed: int   
        :param params: Константы   
        :param points_num: Количество точек траектории  
        :type points_num: int
        """

        rec_speed_calc_noisy = {
            f"noise_{noise_arr_speed[_]}": [
                [0 for i in range(3)] for j in range(rep_num_speed)
            ]
            for _ in range(noise_num_speed)
        }

        rec_speed_calc_noisy_loc = {
            f"noise_{noise_arr_speed[_]}": [
                [0 for i in range(3)] for j in range(rep_num_speed)
            ]
            for _ in range(noise_num_speed)
        }

        if self.num != points_num:

            for noise in rec_speed_calc_noisy.keys():
                for i in range(rep_num_speed):
                    rec_speed_calc_noisy[noise][i] = func_rec_speed(
                        sat_num,
                        self.sat_state,
                        self.H,
                        self.pseudo_dop_noisy[noise][i],
                        self.lam,
                        params,
                    )

                    if self.num == 1:
                        rec_speed_calc_noisy_loc[noise][i] = get_nhe_speed(
                            rec_speed_calc_noisy[noise][i], self.rec_coord_calc
                        )

        self.rec_speed_calc_noisy = rec_speed_calc_noisy
        self.rec_speed_calc_noisy_loc = rec_speed_calc_noisy_loc

    def _rec_speed_mean(self, noise_arr_speed: list[float], noise_num_speed: int) -> None:
        """
        Вычисляет среднее значение компонент вектора скорости приемника для различных шумов

        :param noise_arr_speed: Дисперсии различных шумов  
        :type noise_arr_speed: list[float]
        :param noise_num_speed: Количество различных шумов  
        :type noise_num_speed: int      
        """
        rec_speed_mean = {
            f"noise_{noise_arr_speed[i]}": [0 for j in range(3)]
            for i in range(noise_num_speed)
        }

        for noise in rec_speed_mean.keys():
            for i in range(3):
                temp = np.array(self.rec_speed_calc_noisy[noise]).T[i]
                rec_speed_mean[noise][i] = stat.mean(temp)

        self.rec_speed_mean = rec_speed_mean

    def _rec_sigma(self, 
                   noise_arr: list[float], 
                   noise_num: int, 
                   speed: bool = False) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
        """
        Вычисляет стандартное отклонение координат и компонент вектора скорости приемника 
        в геодезической и местной системах координат

        :param noise_arr: Дисперсии различных шумов  
        :type noise_arr: list[float]
        :param noise_num: Количество различных шумов  
        :type noise_num: int
        :param speed: Флаг: False - рассчитывает стандартное отклонение координат приемника, 
        True - рассчитывает стандартное отклонение компонент вектора скорости приемника
        :type speed: bool
        """

        rec_sigma = {
            f"noise_{noise_arr[i]}": [0 for j in range(3)] for i in range(noise_num)
        }

        rec_sigma_loc = {
            f"noise_{noise_arr[i]}": [0 for j in range(3)] for i in range(noise_num)
        }
        if not speed:
            rec_mean = self.rec_coord_mean
            rec_calc_noisy = self.rec_coord_calc_noisy
        else:
            rec_mean = self.rec_speed_mean
            rec_calc_noisy = self.rec_speed_calc_noisy
        for noise in rec_mean.keys():
            for i in range(3):
                temp = np.array(rec_calc_noisy[noise]).T[i]
                rec_sigma[noise][i] = stat.pstdev(temp)

                if self.num == 1:
                    temp1 = np.array(rec_calc_noisy[noise]).T[i]
                    rec_sigma_loc[noise][i] = stat.pstdev(temp1)

        return rec_sigma, rec_sigma_loc
