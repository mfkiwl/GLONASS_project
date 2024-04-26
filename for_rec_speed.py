from ephemeris import Ephemeris
import numpy as np


def dist_change(rec_coord: list[float], rec_speed: list[float], sat_state: list[float]) -> float:
    """
    :param rec_coord: координаты приемника
    :param rec_speed: скорость приемника
    :param sat_state: "состояние" спутника
    :return: скорость изменения расстояния между спутником и приемником
    """
    r = [sat_state[i] - rec_coord[i] for i in range(3)]  # вектор из приемника в спутник
    mod_r = (r[0]**2 + r[1]**2 + r[2]**2) ** 0.5
    res = sum(r[i] * (rec_speed[i] - sat_state[i + 3]) for i in range(3)) / mod_r
    return res


def get_freq(sat: Ephemeris):
    f0 = 1602*1e6
    df = 562.5*1e3
    freq = f0 + sat.frqNum*df
    return freq


def get_pseudo_dop_not_mod(D, lam, k, df):
    '''
    :param D: скорость изменения расстояния между спутником и приемником
    :param lam: длина волны несущей частоты спцтникого сигнала
    :param k: коэффициент - отношение частоты сутникого сигнала к частоте задающего генератора
    :param df: нестабильность задающего генератора
    '''
    res = -D/lam - k * df
    return res


def get_pseudo_dop(rec_coord: list[float],
                       rec_speed: list[float],
                       sat: Ephemeris,
                       sat_state: list[float],
                       params: dict[str: float]):
    """
    :param rec_coord: Координаты приемника
    :param rec_speed: Скорость приемника
    :param sat: Спутник
    :param sat_state: Состояние спутника
    :param params: Все числа, нужные для расчетов
    :return: Пседодоплеровское смещение для данной скорости приемника и данного спутника, Длину волны спутникого сигнала
    """

    r = [sat_state[i] - rec_coord[i] for i in range(3)]  # вектор из приемника в спутник
    mod_r = (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** 0.5
    D = sum(r[i] * (rec_speed[i] - sat_state[i + 3]) for i in range(3)) / mod_r

    sat_freq = params['f0'] + sat.frqNum * params['df']
    lam = params['c'] / sat_freq
    k = sat_freq / params['f_mo']

    pseudo_dop = - D / lam - k * params['df_mo']
    return pseudo_dop, lam


def func_rec_speed(sat_num: int, sat_state, H, pseudo_dop: list[float], lam: list[float], params: dict[str: float]):
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

