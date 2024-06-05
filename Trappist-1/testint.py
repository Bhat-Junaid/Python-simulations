##################################################################################################################
################################            Simulation                  ##########################################
##################################################################################################################

# Importing all the neccessary dependencies 
import time
#import cProfile
from rk5 import rk5               # Runge-kutta 5th order for initial values
from adamspc10 import adpc11      # ABM-11, integrator
from force_func import ncorps     # Gravitational forces in the system
from Conversion import cartoelt   # Cartesian to elemental coordinates conversion
from bodies import Init           # Bodies in the system


start_time = time.perf_counter()       # Optional, time counter, just to check how fast it is

def main():
    h = 1500.               # Main step size, 1.15 percent of smallest planet
    hp = h / 10.            # smaller time step for RK5 implementation
    t_initial = 0.          # Initial time
    NN = 100000             # Integration steps
    err = 4.e-6             # Error tolerance for Predictor-Corrector

    # Gravitational parameter of all the bodies
    GM = Init()

    # Initial Conditions for positions and velocities of planets
    # x, y, z, vx, vy, vz for each planet
    y0 = [
        1.727e9, 0., 0., 0., 83384.07, 0.,
        2.365e9, 0., 0., 0., 71072.18, 0.,
        3.356e9, 0., 0., 0., 59747.91, 0.,
        4.380e9, 0., 0., 0., 52427.03, 0.,
        5.763e9, 0., 0., 0., 45201.49, 0.,
        7.010e9, 0., 0., 0., 41447.00, 0.,
        9.010e9, 0., 0., 0., 36544.36, 0.
    ]

    # Open files to save the orbital elements for each planet
    f1 = open("planetb.dat", "w")
    f2 = open("planetc.dat", "w")
    f3 = open("planetd.dat", "w")
    f4 = open("planete.dat", "w")
    f5 = open("planetf.dat", "w")
    f6 = open("planetg.dat", "w")
    f7 = open("planeth.dat", "w")

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

        # Convert Cartesian coordinates to orbital elements for each planet
        ellplanetb = cartoelt(y10[:3], y10[3:6], GM[0] + GM[1])
        ellplanetc = cartoelt(y10[6:9], y10[9:12], GM[0] + GM[2])
        ellplanetd = cartoelt(y10[12:15], y10[15:18], GM[0] + GM[3])
        ellplanete = cartoelt(y10[18:21], y10[21:24], GM[0] + GM[4])
        ellplanetf = cartoelt(y10[24:27], y10[27:30], GM[0] + GM[5])
        ellplanetg = cartoelt(y10[30:33], y10[33:36], GM[0] + GM[6])
        ellplaneth = cartoelt(y10[36:39],y10[39:42], GM[0] + GM[7])
        
        # Write orbital elements to files every jth step
        if j == 9:
            f1.write(
                f"{t} {ellplanetb[0]} {ellplanetb[1]} {ellplanetb[2]} {ellplanetb[3]} {ellplanetb[4]} {ellplanetb[5]}\n")
            f2.write(
                f"{t} {ellplanetc[0]} {ellplanetc[1]} {ellplanetc[2]} {ellplanetc[3]} {ellplanetc[4]} {ellplanetc[5]}\n")
            f3.write(
                f"{t} {ellplanetd[0]} {ellplanetd[1]} {ellplanetd[2]} {ellplanetd[3]} {ellplanetd[4]} {ellplanetd[5]}\n")
            f4.write(
                f"{t} {ellplanete[0]} {ellplanete[1]} {ellplanete[2]} {ellplanete[3]} {ellplanete[4]} {ellplanete[5]}\n")
            f5.write(
                f"{t} {ellplanetf[0]} {ellplanetf[1]} {ellplanetf[2]} {ellplanetf[3]} {ellplanetf[4]} {ellplanetf[5]}\n")
            f6.write(
                f"{t} {ellplanetg[0]} {ellplanetg[1]} {ellplanetg[2]} {ellplanetg[3]} {ellplanetg[4]} {ellplanetg[5]}\n")
            f7.write(
                f"{t} {ellplaneth[0]} {ellplaneth[1]} {ellplaneth[2]} {ellplaneth[3]} {ellplaneth[4]} {ellplaneth[5]}\n")
            j = 0
        else:
            j += 1

    # Close all files
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()
    f7.close()

if __name__ == "__main__":
    main()                # Run the main simulation function


# Optional, Record the end time and print the elapsed time
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print("Elapsed Time: ", elapsed_time, " sec")
#cProfile.run('main()')