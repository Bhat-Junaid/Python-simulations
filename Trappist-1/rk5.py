########################################################################################################################
#####################                  Runge Kutta 5th Order                            ################################
########################################################################################################################


def rk5(y, x, h, deriv):
    n = len(y)

    # Initialize vectors for intermediate calculations (k1,k2,k3,k4,k5,k6)
    dydx1 = deriv(x, y)
    dydx2 = deriv(x + h / 4.0, [y[i] + h / 4.0 * dydx1[i] for i in range(n)])
    dydx3 = deriv(x + 0.375 * h, [y[i] + 0.09375 * h * dydx1[i] + 0.28125 * h * dydx2[i] for i in range(n)])
    dydx4 = deriv(x + 12.0 / 13.0 * h, [y[i] + 1932.0 / 2197.0 * h * dydx1[i] - 7200.0 / 2197.0 * h * dydx2[i] +
                                        7296.0 / 2197.0 * h * dydx3[i] for i in range(n)])
    dydx5 = deriv(x + h, [y[i] + h * (439.0 / 216.0 * dydx1[i] - 8.0 * dydx2[i] + 3680.0 / 513.0 * dydx3[i] -
                                      845.0 / 4104.0 * dydx4[i]) for i in range(n)])
    dydx6 = deriv(x + h / 2.0, [y[i] + h * (-8.0 / 27.0 * dydx1[i] + 2.0 * dydx2[i] - 3544.0 / 2565.0 * dydx3[i] +
                                            1859.0 / 4104.0 * dydx4[i] - 11.0 / 40.0 * dydx5[i]) for i in range(n)])

    # Compute the next step using the Runge-Kutta 5th order formula
    yout = [y[i] + h * (16.0 / 135.0 * dydx1[i] + 6656.0 / 12825.0 * dydx3[i] + 28561.0 / 56430.0 * dydx4[i] -
                        9.0 / 50.0 * dydx5[i] + 2.0 / 55.0 * dydx6[i]) for i in range(n)]

    return yout
