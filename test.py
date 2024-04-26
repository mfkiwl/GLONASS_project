import math
import numpy as np
import matplotlib.pyplot as plt
from for_modeling import get_nhe


r = 6378.136 * 1e3  # Радиус Земли в м

# Определение констант для преобразования градусов в радианы
deg_to_rad = math.pi / 180.0
deg_to_rad = 1

# Ввод координат широты и долготы
# latitude = float(input("Введите широту в градусах: "))
# longitude = float(input("Введите долготу в градусах: "))
lat1 = 0.97297
lat2 = 0.972972
long1 = 0.65839
long2 = 0.65839

# Радиус Земли в метрах
R = 6378136.0

# Преобразование координат из географической в геодезическую систему
x1 = R * math.cos(lat1 * deg_to_rad) * math.cos(long1 * deg_to_rad)
y1 = R * math.cos(lat1 * deg_to_rad) * math.sin(long1 * deg_to_rad)
z1 = R * math.sin(lat1 * deg_to_rad)
x2 = R * math.cos(lat2 * deg_to_rad) * math.cos(long2 * deg_to_rad)
y2 = R * math.cos(lat2 * deg_to_rad) * math.sin(long2 * deg_to_rad)
z2 = R * math.sin(lat2 * deg_to_rad)

# print(f"Координаты в геодезической системе: X={x1}, Y={y1}, Z={z1}")
# print(f"Координаты в геодезической системе: X={x2}, Y={y2}, Z={z2}")
# print (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5)


# noise_sigma = 0.001  # сигма нормального распределения
# noise_sigma = 0.002
# noise = np.random.normal(loc=0.0, scale=noise_sigma, size=10000) # Передается именно сигма, не сигма**2
# print(noise)
# plt.hist(noise, color='blue', edgecolor='black')#, bins=int(180/5))
# plt.show()

# plt.plot([1, 2, 3], [4, 5, 6])
# plt.show()

# print(get_nhe([0, 0, r], [r, 0, 0]))
