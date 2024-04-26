# расчет показаний часов спутника на момент предшествия
def func_sat_time(rec_time, pseudo_time):
    res = rec_time - pseudo_time
    return res


# расчет показаний часов МДВ на момент предшествия
def func_mos_time(sat_time, tb, tau, gamma, tauSys):
    res = sat_time + tau - gamma * (sat_time - tb) + tauSys
    while res < 0:
        res += 86400
    return res
