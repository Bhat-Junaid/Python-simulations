import numpy as np
import time
from rk5 import rk5
from adamspc10 import adpc11
from force_func import ncorps
from Conversion import cartoelt
from numba import njit
from bodies import Init

start_time = time.perf_counter()

@njit
def compute_trajectory():
    h = 8640.                # 1.15 percent of smallest planet
    hp = h / 100.
    t_initial = 0.
    NN = 365250
    yr = 21915
    lyr = NN - yr
    err = 4.e-6
    t0 = 86400 * 365.25


    GM = Init()
    y0 = [
        7.78479e11, 0.0, 0.0,0.0, 13.07e3, 0.0,   # Jupiter
        2.57846487e+11, 3.28055458e+11, -2.87016173e+10, -14389.41688278, 9884.18537063,3105.78091161     # Ceres
    ]

    trajectory_data = []
    t = t_initial
    dydx0 = ncorps(t, y0)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx1 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx2 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx3 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx4 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx5 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx6 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx7 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx8 = ncorps(t, y1)
    for i in range(10):
        y1 = rk5(y0, t, hp, ncorps)
        t += hp
        y0 = y1

    dydx9 = ncorps(t, y1)
    y9 = y1
    j = 0
    k = 0

    for i in range(NN):
        y10 = adpc11(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, y9, err, t, h, ncorps)
        y9 = y10
        dydx0 = dydx1
        dydx1 = dydx2
        dydx2 = dydx3
        dydx3 = dydx4
        dydx4 = dydx5
        dydx5 = dydx6
        dydx6 = dydx7
        dydx7 = dydx8
        dydx8 = dydx9
        t += h
        dydx9 = ncorps(t, y9)

        ellplanetb = cartoelt(y10[:3], y10[3:6], GM[0] + GM[1])
        ellplanetc = cartoelt(y10[6:9], y10[9:12], GM[0] + GM[2])

        if i <=yr:
            if j == 9:
                trajectory_data.append((t/t0, ellplanetb))
                trajectory_data.append((t/t0, ellplanetc))
                j = 0
            else:
                j += 1
        elif yr < i <= lyr:
            if j == 999:
                trajectory_data.append((t / t0, ellplanetb))
                trajectory_data.append((t / t0, ellplanetc))
                j = 0
            else:
                j += 1
        elif i > lyr:
            if k == 9:
                trajectory_data.append((t / t0, ellplanetb))
                trajectory_data.append((t / t0, ellplanetc))
                k = 0
            else:
                k += 1

    return trajectory_data

def write_to_files(trajectory_data):
    filenames = ["jupiter.dat", "ceres.dat"]
    for i, filename in enumerate(filenames):
        with open(filename, "w") as file:
            for t, ell in trajectory_data[i::len(filenames)]:
                file.write(f"{t} {' '.join(map(str, ell))}\n")

trajectory_data = compute_trajectory()
write_to_files(trajectory_data)

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print("Elapsed Time: ", elapsed_time, " sec")

