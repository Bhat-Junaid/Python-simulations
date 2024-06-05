############################# Planetary System ################################

# Gravitational constant in m^3 kg^-1 s^-2, G = 6.6743e-11

def Init():
    GM = [
        1.989e30 * 0.09 * 6.6743e-11,  # star
        1.374 * 5.972e24 * 6.6743e-11,    # planet b
        1.308 * 5.972e24 * 6.6743e-11,    # planet c
        0.388 * 5.972e24 * 6.6743e-11,    # planet d
        0.692 * 5.972e24 * 6.6743e-11,    # planet e
        1.039 * 5.972e24 * 6.6743e-11,    # planet f
        1.321 * 5.972e24 * 6.6743e-11,    # planet g
        0.326 * 5.972e24 * 6.6743e-11     # planet h
    ]                                     # G * mass of each body
    return GM       # Return the list of gravitational parameters 
