from consts import *
from class_Point import Point


def modeling(points_num: int, 
             noise_pseudo_time: list[float], rep_num_pseudo_time: int, 
             noise_pseudo_dop: list[float], rep_num_pseudo_dop: int) -> list[Point]:
    """
    Моделирует движение приемника (координаты и компоненты вектора скорости), 
    а также значения производимых в нем измерений псевдозадержек и псевдодоплеровских смещений. 
    Производит расчет параметров движения приемника по смоделированным измерениям.  
    Все расчеты для каждой точки хранятся в классе Point.

    :param points_num: Количество точек траектории
    :type points_num: int
    :param noise_pseudo_time: Массив дисперсий шума, который будет добавлен к смоделированным значениям псевдозадержек 
    :type noise_pseudo_time: list[float]
    :param rep_num_pseudo_time: Количество повторений расчетов координат приемника 
    :type rep_num_pseudo_time: int
    :param noise_pseudo_dop: Массив дисперсий шума, который будет добавлен к смоделированным значениям псевдодоплеровских смещений  
    :type noise_pseudo_dop: list[float]
    :param rep_num_pseudo_dop: Количество повторений расчетов компонент вектора скорости приемника
    :type rep_num_pseudo_dop: int
    :return: Массив элементов класса Point с произведенными рассчетами  
    :rtype: list[Point] 
    """
    
    time_step: int = 5  # Шаг по времени, [c]
    lat_step: float = 0.000002 # Шаг по широте, [рад]

    points = [Point(num=i + 1, time_step=time_step, lat_step=lat_step) for i in range(points_num)]

    for point_prev, point_next in zip(points[:-1], points[1:]):
        point_prev.state.init_speed(point_next, time_step=time_step)

    for point in points:
        point.set_sat_state()  
        point.set_tpr()  
        point._T_sat()  
        point._pseudo_time()  
        point._pseudo_time_noisy(
            len(noise_pseudo_time), rep_num_pseudo_time, noise_pseudo_time
        )  
        point._sat_state_calc()  
        point._sat_state_calc_noisy(
            len(noise_pseudo_time), rep_num_pseudo_time, noise_pseudo_time
        )  
        point._rec_coord_calc()
        point._rec_coord_noisy(rep_num_pseudo_time, noise_pseudo_time, len(noise_pseudo_time))
        point._rec_coord_mean(noise_pseudo_time, len(noise_pseudo_time))
        if point.num != points_num:
            point._pseudo_dop(params, points_num)
            point._pseudo_dop_noisy(len(noise_pseudo_dop), noise_pseudo_dop, rep_num_pseudo_dop)
            point._rec_speed_calc(params)
            point._rec_speed_calc_noisy(
                noise_pseudo_dop, rep_num_pseudo_dop, len(noise_pseudo_dop), params, points_num
            )
            point._rec_speed_mean(noise_pseudo_dop, len(noise_pseudo_dop))

    return points


def sigma(points: list[Point], 
          points_num: int, 
          noise_arr_coord: list[float], 
          noise_num_coord: int, 
          noise_arr_speed: list[float], 
          noise_num_speed: int) -> tuple[list[dict[str, list[float]]]]:
    
    points_coord_sigma = []
    points_coord_sigma_loc = []
    points_speed_sigma = []
    points_speed_sigma_loc = []
    for point in points:
        coord_sigma, coord_sigma_loc = point._rec_sigma(
            noise_arr_coord, noise_num_coord
        )
        points_coord_sigma.append(coord_sigma)
        points_coord_sigma_loc.append(coord_sigma_loc)
        if point.num != points_num:
            speed_sigma, speed_sigma_loc = point._rec_sigma(
                noise_arr_speed, noise_num_speed, True
            )
            points_speed_sigma.append(speed_sigma)
            points_speed_sigma_loc.append(speed_sigma_loc)

    return points_coord_sigma, points_speed_sigma, points_coord_sigma_loc, points_speed_sigma_loc
