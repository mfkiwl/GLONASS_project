from modeling import modeling, sigma
from visualization import plot_comp, plot_std_dev


# Задаем количество дисперсий различных шумов,
# которые будут добавлены к измерениям псевдозадержек и псевдодоплеровских смещений соответственно:
noise_num_pseudo_time = 50
noise_num_pseudo_dop = 50

# Задаем дисперсии этих шумов:
noises_pseudo_time = [(i + 1) * 3e-5 for i in range(noise_num_pseudo_time)]
noises_pseudo_dop = [(i + 1) * 1e-4 for i in range(noise_num_pseudo_dop)]

points_num = 10  # Задаем количество точек траектории движения приемника
trajectory = modeling(
    points_num=points_num,
    noise_pseudo_time=noises_pseudo_time,
    rep_num_pseudo_time=50,
    noise_pseudo_dop=noises_pseudo_dop,
    rep_num_pseudo_dop=50,
)

coord_sigma, speed_sigma, coord_sigma_loc, speed_sigma_loc = sigma(
    trajectory,
    points_num,
    noises_pseudo_time,
    len(noises_pseudo_time),
    noises_pseudo_dop,
    len(noises_pseudo_dop),
)

plot_comp(trajectory, noises_pseudo_time, 0, points_num, 0, True, True)
plot_std_dev(coord_sigma[0], noises_pseudo_time, len(noises_pseudo_time), True, True, True)
