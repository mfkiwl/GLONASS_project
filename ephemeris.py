import json

class Ephemeris:
    def __init__(self, health=None, age=None, flags=None, sv=None, frqNum=None, dne=None, tk=None, tb=None,
                 r=None, v=None, w=None, tauSys=None, tau=None, gamma=None, name=None):
        self.health = health  # Состояние спутника
        self.age = age  # Возраст эфемеридной информации
        self.flags = flags  # Флаги
        self.sv = sv  # Номер слота орбиты спутника от 1 до 24
        self.frqNum = frqNum  # Номер частотного канала спутника от –7 до +6

        self.dne = dne  # Номер суток ГЛОНАСС внутри четырехлетнего периода
        self.tk = tk  # Время начала кадра в текущих сутках ГЛОНАСС [с]
        self.tb = tb  # Узловой момент времени эфемерид по шкале МДВ [с] внутри суток dne

        self.r = r  # Координаты спутника в системе PЗ-90 на момент tb
        self.v = v  # Скорости спутника в системе PЗ-90 на момент tb
        self.w = w  # Ускорения спутника в системе PЗ-90 на момент tb

        self.tauSys = tauSys  # Поправка к шкале времени ГЛОНАСС для перехода на шкалу МДВ (tauSys = Тмдв - Тглн [с])
        self.tau = tau  # Поправка к шкале времени спутника для перехода на шкалу ГЛОНАСС (tau = Тглн - Тка [с])
        self.gamma = gamma  # Относительное отклонение прогнозируемого значения несущей частоты на момент tb

        self.name = name  # "Название" спутника

    def __str__(self) -> str:
        return self.name
    
    def __repr__(self) -> str:
        return self.name


def read_sat_info():
    result = []
    with open('GLONASS_project/ephemeris.json', 'r') as f:
        data = json.load(f)
        for sat in data.keys():
            sat_obj = Ephemeris(**data[sat])
            result.append(sat_obj)
    return result

sat_arr =read_sat_info()
sat_num = len(sat_arr)
