import math
from consts import *


class Model:
    R = R * 1e3
    c = c * 1e3

    def __init__(self, num: int, lat_step:float) -> None:
        """
        Моделирует истинные координаты и компоненты вектора скорости приемника 

        :param num: Номер точки траектории движения приемника
        :type num: int
        :param num: Шаг по широте, [рад]
        :type num: float
        :return: Объект класса Model, содержащий географические и геодезические координаты приемника, а также компоненты его вектора скорости
        :rtype: State
        """
        self.num = num
        self.x_speed, self.y_speed, self.z_speed, self.speed = 0, 0, 0, []
        self.geographic_coord = self._geography(lat_step)
        self._geodezic()
        
    def __getitem__(self, i: int):
        return self.coord[i]

    def _geodezic(self) -> tuple[float]:
        """
        Преобразует координаты приемника из географической системы координат в геодезические

        :return: Координаты приемника в геодезической системе координат 
        :rtype: tuple[float]
        """
        self.x = (
            self.R * math.cos(self.geographic_coord[0]) * math.cos(self.geographic_coord[1])
        )
        self.y = (
            self.R * math.cos(self.geographic_coord[0]) * math.sin(self.geographic_coord[1])
        )
        self.z = self.R * math.sin(self.geographic_coord[0])
        self.coord = [self.x, self.y, self.z]
        return (self.x, self.y, self.z)

    def _geography(self, lat_step):
        """
        Рассчитывает координаты точки траектории движения приемника в географической системе координат

        :param num: Номер точки траектории движения приемника
        :type num: int
        :return: Координаты точки в геогрфической системе координат (широта, долгота)
        :rtype: list[float]
        """
        return [0.95993108859 + lat_step * (self.num - 1), 0.65839]

    def init_speed(self, next_point, time_step: int) -> tuple[float]:
        """
        Рассчитывает компоненты вектора скорости приемника

        :param next_point: Следующая точка траектории
        :type next_point: Point
        :param next_point: Шаг моделирования по времени
        :type next_point: int
        :return: Компоненты вектора скорости приемника
        :rtype: tuple[float]
        """
        self.x_speed = (next_point.state[0] - self[0]) / time_step
        self.y_speed = (next_point.state[1] - self[1]) / time_step
        self.z_speed = (next_point.state[2] - self[2]) / time_step
        self.speed = [self.x_speed, self.y_speed, self.z_speed]
        return (self.x_speed, self.y_speed, self.z_speed, self.speed)
    