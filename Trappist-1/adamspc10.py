########################################################################################################################
#####################             Adams Bashforth Moulton Predictor Corrector           ################################
########################################################################################################################

from adamsb10 import ab10
from adamsm10 import am11


def adpc11(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, y9, err, x, h, deriv):
    N = len(dydx0)
     # Check if all input arrays have the same length
    if any(len(arr) != N for arr in [dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, y9]):
        print("Error in adpc11: wrong size")
        exit(1)
    
    # Predictor step: Adams-Bashforth 10th-order method
    y10b = ab10(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, y9, h)
    dydx10 = deriv(x + h, y10b)

    # Corrector step: Adams-Moulton 11th-order method
    y10a = am11(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, dydx10, y9, h)
    
    # Calculate the difference between predictor and corrector
    diff = 0.0
    for i in range(N):
        if abs(y10a[i] - y10b[i]) > diff:
            diff = abs(y10a[i] - y10b[i])
 
    # Iterate until the difference is within the acceptable error tolerance
    while diff > err:
        y10b = y10a
        dydx10 = deriv(x + h, y10b)
        y10a = am11(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, dydx10, y9, h)
        diff = 0.0
        for i in range(N):
            if abs(y10a[i] - y10b[i]) > diff:
                diff = abs(y10a[i] - y10b[i])

    return y10b    # Return the corrected solution
