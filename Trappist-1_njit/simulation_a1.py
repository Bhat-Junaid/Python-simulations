# Importing all the neccessary dependencies 
import time
from rk5 import rk5                       # Runge-kutta 5th order for initial values
from adamspc10 import adpc11              # ABM-11, integrator
from force_func import ncorps             # Gravitational forces in the system
from Conversion import cartoelt, eletocar # Cartesian to elemental coordinates conversion
from bodies import Init                   # Bodies in the system
from numba import njit


start_time = time.perf_counter()        # Optional, time counter, just to check how fast it is

@njit
def compute_trajectory():
    h = 1500.               # Main step size, 1.15 percent of smallest planet
    hp = h / 10.            # smaller time step for RK5 implementation
    t_initial = 0.          # Initial time 
    NN = 21038400           # integration Steps
    yr = 21038              #  1 earth year in terms of number of steps
    lyr = NN - yr           #  Last year of simulation
    err = 4.e-6             # Error tolerance for Predictor-Corrector

    # Gravitational parameter of all the bodies
    GM = Init()
 
    # Initial Conditions for positions and velocities of planets
    # [a,h,k,q,p,\lambda]
    elem = [
        1741719667.5325758, -0.00215, 0.00217, 0.0000005, 0.0, 1.5834738519202087,
        2352198529.2521772, 0.00055, 0.00001, 0.0000005, 0.0, -1.7011231035507217,
        3331462499.7237353, -0.00496, 0.00267, 0.0000005, 0.0, 1.323972877958588,
        4395196370.392811, 0.00433, -0.00461, 0.0000005, 0.0, 0.10674692636961337,
        5512261142.634852, -0.00840, -0.00051, 0.0000005, 0.0, 0.5837847930221738,
        6966704751.692981, 0.00380, 0.00128, 0.0000005, 0.0, 0.11047396920875246,
        9077291856.855492, -0.00365, -0.00002, 0.0000005, 0.0, 3.0423518820022775
    ]

    def convert(elem, GM):
        y0 = []

        # process it in chunks of 6 elements
        for i in range(0, len(elem), 6):
            ell_chunk = elem[i:i + 6]
            pos, vit = eletocar(ell_chunk, GM)

            # Append position and velocity to y0 in the specified format
            y0.extend([pos[0], pos[1], pos[2], vit[0], vit[1], vit[2]])

        return y0

    # Convert initial orbital elements to Cartesian coordinates
    y0 = convert(elem, GM[0])

    # List to store the trajectory data
    trajectory_data = []
    t = t_initial
    dydx0 = ncorps(t, y0)

     # Use RK5 to initialize the first 10 steps
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
    k = 0 # Parameter to record points in file


    # Main integration loop using Adams-Bashforth-Moulton predictor-corrector method
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

        #Cartesian coordinates to orbital elements 
        ellplanetb = cartoelt(y10[:3], y10[3:6], GM[0] + GM[1])
        ellplanetc = cartoelt(y10[6:9], y10[9:12], GM[0] + GM[2])
        ellplanetd = cartoelt(y10[12:15], y10[15:18], GM[0] + GM[3])
        ellplanete = cartoelt(y10[18:21], y10[21:24], GM[0] + GM[4])
        ellplanetf = cartoelt(y10[24:27], y10[27:30], GM[0] + GM[5])
        ellplanetg = cartoelt(y10[30:33], y10[33:36], GM[0] + GM[6])
        ellplaneth = cartoelt(y10[36:39],y10[39:42], GM[0] + GM[7])

        # Record points for first year with a small sampling rate
        if i <= yr:
            if j == 9:
                trajectory_data.append((t, ellplanetb))
                trajectory_data.append((t, ellplanetc))
                trajectory_data.append((t, ellplanetd))
                trajectory_data.append((t, ellplanete))
                trajectory_data.append((t, ellplanetf))
                trajectory_data.append((t, ellplanetg))
                trajectory_data.append((t, ellplaneth))
                j = 0
            else:
                j += 1
        # Record points between first and last year with different sampling rate
        elif yr < i <= lyr:
            if j == 21038:
                trajectory_data.append((t, ellplanetb))
                trajectory_data.append((t, ellplanetc))
                trajectory_data.append((t, ellplanetd))
                trajectory_data.append((t, ellplanete))
                trajectory_data.append((t, ellplanetf))
                trajectory_data.append((t, ellplanetg))
                trajectory_data.append((t, ellplaneth))
                j = 0
            else:
                j += 1

        # Record points for last year with a small sampling rate
        elif i > lyr:
            if k == 9:
                trajectory_data.append((t, ellplanetb))
                trajectory_data.append((t, ellplanetc))
                trajectory_data.append((t, ellplanetd))
                trajectory_data.append((t, ellplanete))
                trajectory_data.append((t, ellplanetf))
                trajectory_data.append((t, ellplanetg))
                trajectory_data.append((t, ellplaneth))
                k = 0
            else:
                k += 1


    return trajectory_data


# Stores the trajectory data into files 
def write_to_files(trajectory_data):
    filenames = ["planetb_a1.dat", "planetc_a1.dat", "planetd_a1.dat", "planete_a1.dat", "planetf_a1.dat", "planetg_a1.dat", "planeth_a1.dat"]
    for i, filename in enumerate(filenames):
        with open(filename, "w") as file:
            for t, ell in trajectory_data[i::len(filenames)]:
                file.write(f"{t} {' '.join(map(str, ell))}\n")

trajectory_data = compute_trajectory()
write_to_files(trajectory_data)


# Optional,Record the end time and print the elapsed time
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print("Elapsed Time: ", elapsed_time, " sec")

