########################################################################################################################
#####################                  Adams Bashforth 10th order                       ################################
########################################################################################################################
def ab10(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, y9, h):
    N = len(dydx0)

    # Check if all input arrays have the same length
    if any(len(arr) != N for arr in [dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, y9]):
        print("error in ab10: wrong size")
        exit(1)
 
    # Initialize the next step solution array
    y10 = []

    # Calculate the next step solution using the Adams-Bashforth 10th-order formula
    for i in range(N):
        truc = y9[i] + h / 7257600. * (30277247. * dydx9[i] - 104995189. * dydx8[i] \
                                       + 265932680. * dydx7[i] - 454661776. * dydx6[i] + 538363838. * dydx5[i] \
                                       - 444772162. * dydx4[i] + 252618224. * dydx3[i] - 94307320. * dydx2[i] \
                                       + 20884811. * dydx1[i] - 2082753. * dydx0[i])
        y10.append(truc)
    return y10                       # Return the solution at the next step
