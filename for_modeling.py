import matplotlib.pyplot as plt


def dist(a: list[float], b: list[float]) -> float:
    res = ((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)**0.5
    return res


def plot_comp(
        true: dict[str: list[float]],
        noisy: dict,
        noise_arr: list[float],
        noise_num: int,
        time_step: float,
        tb: float,
        points_num: int,
        index: int,
        flag: bool
):
    """
    :param true: Истинные значения
    :param noisy: Зашумленные значения
    :param noise_arr: Массив шумов
    :param noise_num: Индекс шума (в массиве шумов)
    :param time_step: Шаг по времени
    :param tb: Узловой момент
    :param points_num: Количество точек
    :param index: Компонента, которую хотим отобразить на графике
    :param flag: true - координаты, false - скорость
    :return: строит сравнительный график
    """

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

    t = [tb + time_step*(i + 1) for i in range(points_num)]
    true_y = [true[f'point_{i}'][index] for i in range(1, points_num + 1)]
    noisy_y = [noisy[f'point_{i}'][f'noise_{noise_arr[noise_num]}'][0][index] for i in range(1, points_num + 1)]
    plt.plot(t, true_y, color='red')
    plt.plot(t, noisy_y, color='blue')

    # Строим следующий в списке шум:
    noisy_y = [noisy[f'point_{i}'][f'noise_{noise_arr[noise_num + 1]}'][0][index] for i in range(1, points_num + 1)]
    plt.plot(t, noisy_y, color='green')
    if flag:
        plt.legend([f'Real coordinate (component {comp})',
                    f'Calculated coordinate (noise sigma = {noise_arr[noise_num]})',
                    f'Calculated coordinate (noise sigma = {noise_arr[noise_num + 1]})'])
        plt.ylabel(f'{comp}, m')
        plt.title(f'Comparative plot \n Coordinate (component {comp})' )

    else:
        plt.legend([f'Real speed (component {comp})',
                    f'Calculated speed (noise sigma = {noise_arr[noise_num]})',
                    f'Calculated speed (noise sigma = {noise_arr[noise_num + 1]})'])
        plt.ylabel(f'V{comp}, m/c')
        plt.title(f'Comparative plot \n Speed (component {comp})')

    plt.xlabel('time, c')
    plt.show()


def plot_std_dev(a: dict, noise_arr: list[float], noise_num: int, flag: bool, flag_cs: bool):

    """
    :param a: Сигма ошибок определения координат / скорости
    :param noise_arr: Сигма шумов
    :param noise_num: Количество шумов
    :param flag: true - координаты, false - скорость
    :param flag_cs: true - xyz, false - NHE
    :return:
    """

    y = [a['point_1'][f'noise_{noise_arr[i]}'][0] for i in range(noise_num)]
    plt.plot(noise_arr, y, color='green')

    y = [a['point_1'][f'noise_{noise_arr[i]}'][1] for i in range(noise_num)]
    plt.plot(noise_arr, y, color='blue')

    y = [a['point_1'][f'noise_{noise_arr[i]}'][2] for i in range(noise_num)]
    plt.plot(noise_arr, y, color='red')

    if flag_cs:
        plt.legend(['Component x', 'Component y', 'Component z'])
    else:
        plt.legend(['Component N', 'Component H', 'Component E'])

    if flag:
        plt.ylabel('Standard deviation of the coordinate, m')
        plt.xlabel('Standard deviation of the noise, c')
    else:
        plt.ylabel('Standard deviation of the speed, m/c')
        plt.xlabel('Standard deviation of the noise, Hz')

    plt.show()


def get_nhe_coord(coord: list[float], coord0: list[float]) -> list[float]:
    """
    :param coord: Координаты точки, которые будем переводить в местную СК
    :param coord0: Координаты начальной точки
    :return: Координаты точки в новой СК
    """

    x, y, z = coord[0], coord[1], coord[2]
    x0, y0, z0 = coord0[0], coord0[1], coord0[2]

    N = (-x0*z0*(x - x0) - y0*z0*(y - y0) + (x0**2 + y0**2)*(z - z0)) / ((x0*z0)**2 + (y0*z0)**2 + (x0**2 + y0**2)**2)**0.5
    H = (x0*(x - x0) + y0*(y - y0) + z0*(z - z0)) / (x0**2 + y0**2 + z0**2)**0.5
    E = (-y0*(x - x0) + x0*(y - y0)) / (x0**2 + y0**2)**0.5

    return [N, H, E]


def get_nhe_speed(speed: list[float], coord0: list[float]) -> list[float]:
    """
    :param speed: Вектор скорости, который будем переводить в местную СК
    :param coord0: Координаты начальной точки
    :return: Вектор скорости в местной СК
    """

    point = [i + j for i, j in zip(coord0, speed)]
    res = get_nhe_coord(point, coord0)
    return res

