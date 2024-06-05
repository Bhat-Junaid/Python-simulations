########################################################################################################################
#####################                  Adams Moulton 11 th order                       ################################
########################################################################################################################

def am11(dydx0, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, dydx10, y9, h):
    
    # Number of variables 
    N = len(dydx0)

    # Check if all input arrays have the same length
    if any(len(arr) != N for arr in [dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8, dydx9, dydx10, y9]):
        print("Error in am11: wrong size")
        exit(1)
    
    # Initialize the next step solution array
    y10 = []
    pouet = 2.08767569878680989792e-9

    # Calculate the next step solution using the Adams-Moulton 11th-order formula
    for i in range(N):
        truc = y9[i] + h * pouet * (
                    134211265 * dydx10[i] + 656185652 * dydx9[i] - 890175549 * dydx8[i] + 1446205080 * dydx7[
                i] - 1823311566 * dydx6[i] + 1710774528 * dydx5[i] - 1170597042 * dydx4[i] + 567450984 * dydx3[
                        i] - 184776195 * dydx2[i] + 36284876 * dydx1[i] - 3250433 * dydx0[i])
        y10.append(truc)

    return y10      # Return the solution at the next step
