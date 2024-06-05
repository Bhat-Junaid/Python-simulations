#######################################################################################################################
###### Force function implementation with the central body at the origin    ################################
#######################################################################################################################

import numpy as np
from bodies import Init


def ncorps(t, rrp):                                     # rrp: array containing initial conditions
    global J2, C22, J3, J4, J6, J8, R

    # Perturbation Coefficients (set to zero)
    J2 = 0.0
    C22 = 0.0
    J3 = 0.0
    J4 = 0.0
    J6 = 0.0
    J8 = 0.0
    R = 0.0
    GM = Init()                                        # Gravitational Parameter

    if len(rrp) != 6 * (len(GM) - 1):
        print("Error in nbodies, the vectors have wrong sizes")
        print(len(rrp), len(GM))
        exit(1)

    nbcorps = len(GM) - 1                              # number of planets
    vecsol = [0.0] * len(rrp)
    unsr3 = []

    # Inverse of cube of distance from star
    for i in range(nbcorps):
        unsr3.append((rrp[i * 6] ** 2 + rrp[i * 6 + 1] ** 2 + rrp[i * 6 + 2] ** 2) ** -1.5)

    # Force due to the Star
    for i in range(nbcorps):
        x, y, z = rrp[i * 6], rrp[i * 6 + 1], rrp[i * 6 + 2]
        vecsol[i * 6] += rrp[i * 6 + 3]
        vecsol[i * 6 + 1] += rrp[i * 6 + 4]
        vecsol[i * 6 + 2] += rrp[i * 6 + 5]
        vecsol[i * 6 + 3] -= (GM[0] + GM[i + 1]) * x * unsr3[i]
        vecsol[i * 6 + 4] -= (GM[0] + GM[i + 1]) * y * unsr3[i]
        vecsol[i * 6 + 5] -= (GM[0] + GM[i + 1]) * z * unsr3[i]

        # Force due to all other bodies
        for j in range(nbcorps):
            if i != j:
                unsdist3 = ((rrp[i * 6] - rrp[j * 6]) ** 2 + (rrp[i * 6 + 1] - rrp[j * 6 + 1]) ** 2 + (
                            rrp[i * 6 + 2] - rrp[j * 6 + 2]) ** 2) ** -1.5
                vecsol[i * 6 + 3] += GM[j + 1] * ((rrp[j * 6] - rrp[i * 6]) * unsdist3 - rrp[j * 6] * unsr3[j])
                vecsol[i * 6 + 4] += GM[j + 1] * (
                            (rrp[j * 6 + 1] - rrp[i * 6 + 1]) * unsdist3 - rrp[j * 6 + 1] * unsr3[j])
                vecsol[i * 6 + 5] += GM[j + 1] * (
                            (rrp[j * 6 + 2] - rrp[i * 6 + 2]) * unsdist3 - rrp[j * 6 + 2] * unsr3[j])

        truc = 1.0 / np.sqrt(x * x + y * y + z * z)
        xn, yn, zn = x * truc, y * truc, z * truc

        # Irrevelant for the current computations
        # Takes in account the relativistic factors
        # Can be easily commented off and the code still works
        
        xJ2 = -1.5 * R ** 2 * xn * (J2 * (xn * xn + yn * yn - 4 * zn * zn) + C22 * (
                    6 * xn * xn - 14 * yn * yn - 4 * zn * zn)) * truc ** 4
        yJ2 = -1.5 * R ** 2 * yn * (J2 * (xn * xn + yn * yn - 4 * zn * zn) + C22 * (
                    14 * xn * xn - 6 * yn * yn + 4 * zn * zn)) * truc ** 4
        zJ2 = -1.5 * R ** 2 * zn * (
                    J2 * (3 * xn * xn + 3 * yn * yn - 2 * zn * zn) + 10 * C22 * (xn * xn - yn * yn)) * truc ** 4

        xJ3 = -2.5 * R ** 3 * J3 * xn * zn * (3 * (xn * xn + yn * yn) - 4 * zn * zn) * truc ** 5
        yJ3 = -2.5 * R ** 3 * J3 * yn * zn * (3 * (xn * xn + yn * yn) - 4 * zn * zn) * truc ** 5
        zJ3 = 0.5 * R ** 3 * J3 * (3 * (xn ** 4 + yn ** 4) - 24 * zn * zn * (
                    xn * xn + yn * yn) + 6 * xn * xn * yn * yn + 8 * zn ** 4) * truc ** 5

        xJ4 = 15 / 8 * R ** 4 * J4 * xn * (8 * zn ** 4 - 12 * zn * zn * (
                    xn * xn + yn * yn) + 2 * xn * xn * yn * yn + xn ** 4 + yn ** 4) * truc ** 6
        yJ4 = 15 / 8 * R ** 4 * J4 * yn * (8 * zn ** 4 - 12 * zn * zn * (
                    xn * xn + yn * yn) + 2 * xn * xn * yn * yn + xn ** 4 + yn ** 4) * truc ** 6
        zJ4 = 5 / 8 * R ** 4 * J4 * zn * (
                    8 * zn ** 4 - 40 * zn * zn * (xn * xn + yn * yn) + 30 * xn * xn * yn * yn + 15 * (
                        xn ** 4 + yn ** 4)) * truc ** 6

        xJ6 = -7 / 16 * R ** 6 * J6 * xn * (
                    -64 * zn ** 6 + 240 * zn ** 4 * (xn * xn + yn * yn) - 120 * zn * zn * (xn ** 4 + yn ** 4) - 240 * (
                        xn * yn * zn) ** 2 + 5 * xn ** 6 + 15 * xn * xn * yn * yn * (
                                xn * xn + yn * yn) + 5 * yn ** 6) * truc ** 8
        yJ6 = -7 / 16 * R ** 6 * J6 * yn * (
                    -64 * zn ** 6 + 240 * zn ** 4 * (xn * xn + yn * yn) - 120 * zn * zn * (xn ** 4 + yn ** 4) - 240 * (
                        xn * yn * zn) ** 2 + 5 * xn ** 6 + 15 * xn * xn * yn * yn * (
                                xn * xn + yn * yn) + 5 * yn ** 6) * truc ** 8
        zJ6 = -7 / 16 * R ** 6 * J6 * zn * (-16 * zn ** 6 + 168 * zn ** 4 * (xn * xn + yn * yn) - 210 * zn * zn * (
                    xn * xn + yn * yn) ** 2 + 35 * (xn ** 6 + yn ** 6) + 105 * xn * xn * yn * yn * (
                                                        xn * xn + yn * yn)) * truc ** 8

        xJ8 = -45 / 128 * R ** 8 * J8 * xn * (1120 * zn ** 4 * (xn * xn + yn * yn) ** 2 - 840 * (xn * yn * zn) ** 2 * (
                    xn * xn + yn * yn) + 128 * zn ** 8 + 7 * (xn ** 8 + yn ** 8) - 896 * zn ** 6 * (
                                                          xn * xn + yn * yn) - 280 * zn * zn * (
                                                          xn ** 6 + yn ** 6) + 28 * (xn * xn + yn * yn) * (
                                                          xn ** 4 + yn ** 4) + 42 * (xn * yn) ** 4) * truc ** 10
        yJ8 = -45 / 128 * R ** 8 * J8 * yn * (1120 * zn ** 4 * (xn * xn + yn * yn) ** 2 - 840 * (xn * yn * zn) ** 2 * (
                    xn * xn + yn * yn) + 128 * zn ** 8 + 7 * (xn ** 8 + yn ** 8) - 896 * zn ** 6 * (
                                                          xn * xn + yn * yn) - 280 * zn * zn * (
                                                          xn ** 6 + yn ** 6) + 28 * (xn * xn + yn * yn) * (
                                                          xn ** 4 + yn ** 4) + 42 * (xn * yn) ** 4) * truc ** 10
        zJ8 = -9 / 128 * R ** 8 * J8 * zn * (
                    (xn * yn * zn) ** 2 * (12096 * zn ** 2 - 10080 * (xn * xn + yn * yn)) + 128 * zn ** 8 + 315 * (
                        xn ** 8 + yn ** 8) - 2304 * zn ** 6 * (xn * xn + yn * yn) + 6048 * zn ** 4 * (
                                xn ** 4 + yn ** 4) - 3360 * zn * zn * (xn ** 6 + yn ** 6) + 1260 * (xn * yn) ** 2 * (
                                xn ** 4 + yn ** 4) + 1890 * (xn * yn) ** 4) * truc ** 10

        vecsol[i * 6 + 3] += (GM[0] + GM[i + 1]) * (xJ2 + xJ3 + xJ4 + xJ6 + xJ8)
        vecsol[i * 6 + 4] += (GM[0] + GM[i + 1]) * (yJ2 + yJ3 + yJ4 + yJ6 + yJ8)
        vecsol[i * 6 + 5] += (GM[0] + GM[i + 1]) * (zJ2 + zJ3 + zJ4 + zJ6 + zJ8)

    return vecsol
