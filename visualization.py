import matplotlib.pyplot as plt
from class_Point import Point
from ephemeris import sat_arr


def plot_comp(
        points: list[Point],
        noise_arr: list[float],
        noise_num: int,
        points_num: int,
        index: int,
        flag: bool,
        save: bool = False
) -> None:
    """
    Строит сравнительный график истинных и рассчитанных значений координат или компонент вектора скорости от времени 
    
    :param points: Траектория движения приемника
    :type points: list[Point]
    :param noise_arr: Дисперсии шумов 
    :type noise_num: list[float]
    :param noise_num: Индекс шума в массиве шумов
    :type noise_num: int
    :param points_num: Количество точек траектории
    :type points_num: int
    :param index: Компонента, которую хотим отобразить на графике
    :type index: int
    :param flag: True - координаты, False - скорость
    :type flag: bool
    :param save: True - сохранит график в файл, False - не сохранит
    :type save: bool
    """

    time_step = 5
    tb = sat_arr[0].tb 
    
    comp = ''
    if index == 0:
        comp = 'x'
    else:
        pass
    if index == 1:
        comp = 'y'
    else:
        pass
    if index == 2:
        comp = 'z'
    else:
        pass
 
    if flag:
        t = [tb + time_step*(i + 1) for i in range(points_num)]

        file_name = 'Сравнительный график для координат.png'

        true_y = [points[i].rec_coord_calc[index] for i in range(points_num)]
        noisy_y = [points[i].rec_coord_calc_noisy[f'noise_{noise_arr[noise_num]}'][0][index] for i in range(points_num)]

        plt.plot(t, true_y, color='red')
        plt.plot(t, noisy_y, color='blue')

        noisy_y = [points[i].rec_coord_calc_noisy[f'noise_{noise_arr[noise_num + 1]}'][0][index] for i in range(points_num)]

        plt.plot(t, noisy_y, color='green')

        plt.legend([f'Real coordinate (component {comp})',
                    f'Calculated coordinate (noise sigma = {noise_arr[noise_num]})',
                    f'Calculated coordinate (noise sigma = {noise_arr[noise_num + 1]})'])
        plt.ylabel(f'{comp}, m')
        plt.title(f'Comparative plot \n Coordinate (component {comp})' )

    else:
        t = [tb + time_step*(i + 1) for i in range(points_num - 1)]

        file_name = 'Сравнительный график для скоростей.png'

        true_y = [points[i].rec_speed_calc[index] for i in range(points_num - 1)]
        noisy_y = [points[i].rec_speed_calc_noisy[f'noise_{noise_arr[noise_num]}'][0][index] for i in range(points_num - 1)]
        print(true_y)
        print(noisy_y)

        plt.plot(t, true_y, color='red')
        plt.plot(t, noisy_y, color='blue')

        noisy_y = [points[i].rec_speed_calc_noisy[f'noise_{noise_arr[noise_num + 1]}'][0][index] for i in range(points_num - 1)]
        print(noisy_y)

        plt.plot(t, noisy_y, color='green')

        plt.legend([f'Real speed (component {comp})',
                    f'Calculated speed (noise sigma = {noise_arr[noise_num]})',
                    f'Calculated speed (noise sigma = {noise_arr[noise_num + 1]})'])
        plt.ylabel(f'V{comp}, m/c')
        plt.title(f'Comparative plot \n Speed (component {comp})')

    plt.xlabel('time, c')

    if save:
        plt.savefig(f'GLONASS_project/images/{file_name}')

    plt.show()


def plot_std_dev(a: dict[str, list[float]], noise_arr: list[float], noise_num: int, flag: bool, flag_cs: bool, save: bool = False):

    """
    Строит график зависимости стандартного отклонения координат или компонент вектора скорости приемника 
    от стандартного оклонения добаленного в смоделированные измерения псевдозадержек и псевдодоплеровских смещений шума 
    
    :param a: Стандартное отклонение координат или скорости приемника
    :type a: list[Point]
    :param noise_arr: Дисперсии шумов 
    :type noise_arr: list[float]
    :param noise_num: Количество шумов
    :type noise_num: int
    :param flag: True - координаты, False - скорость
    :param flag_cs: True - xyz, False - NHE
    :param save: True - сохранит график в файл, False - не сохранит
    """

    y = [a[f'noise_{noise_arr[i]}'][0] for i in range(noise_num)]
    plt.plot(noise_arr, y, color='green')

    y = [a[f'noise_{noise_arr[i]}'][1] for i in range(noise_num)]
    plt.plot(noise_arr, y, color='blue')

    y = [a[f'noise_{noise_arr[i]}'][2] for i in range(noise_num)]
    plt.plot(noise_arr, y, color='red')

    file_name = ''

    if flag_cs:
        plt.legend(['Component x', 'Component y', 'Component z'])
        file_name += 'xyz_'
    else:
        plt.legend(['Component N', 'Component H', 'Component E'])
        file_name += 'NHE_'

    if flag:
        plt.ylabel('Standard deviation of the coordinate, m')
        plt.xlabel('Standard deviation of the noise, c')
        file_name += 'coordinate_noise.png'
    else:
        plt.ylabel('Standard deviation of the speed, m/c')
        plt.xlabel('Standard deviation of the noise, Hz')
        file_name += 'speed_noise.png'

    if save:
        plt.savefig(f'GLONASS_project/images/{file_name}')

    plt.show()
    