############################# Conversion of Cartesian to elemental cordinates #####################################
import numpy as np

def cartoelt(pos, vit, GM):
    pos = np.array(pos)
    vit = np.array(vit)

    ell = np.zeros(6)                            # Initialize the orbital elements array
    rayon = np.sqrt(np.sum(pos ** 2))            # distance from the star
    v2 = np.sum(vit ** 2)                        # Squared velocity
    a = GM * rayon / (2 * GM - rayon * v2)       # semi major axis
    
    # Angular momentum vector components
    gx = pos[1] * vit[2] - pos[2] * vit[1]
    gy = pos[2] * vit[0] - pos[0] * vit[2]
    gz = pos[0] * vit[1] - pos[1] * vit[0]
    gg = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)     # magnitude of angular momentum
    cis2 = np.sqrt(0.5 * (1 + gz / gg))
    q = -gy / (2 * gg * cis2)                   # sin(I/2)cos(OMEGA)
    p = gx / (2 * gg * cis2)                    # sin(I/2)sin(OMEGA)
    tp = 1 - 2 * p ** 2
    tq = 1 - 2 * q ** 2
    dg = 2 * p * q

    # Cartesian to plane orbital coordinates
    x1 = tp * pos[0] + dg * pos[1] - 2 * p * cis2 * pos[2]
    y1 = dg * pos[0] + tq * pos[1] + 2 * q * cis2 * pos[2]
    vx1 = tp * vit[0] + dg * vit[1] - 2 * p * cis2 * vit[2]
    vy1 = dg * vit[0] + tq * vit[1] + 2 * q * cis2 * vit[2]
    k = gg * vy1 / GM - x1 / rayon                          # ecos(omega + Omega)
    h = -gg * vx1 / GM - y1 / rayon                         #esin(omega + Omega)
    
    # auxilary variables for inclination and node calculations
    psi = 1 / (1 + np.sqrt(1 - k ** 2 - h ** 2))
    ach = 1 - psi * h ** 2
    ack = 1 - psi * k ** 2
    adg = psi * h * k
    det = ach * ack - adg ** 2

    # Compute true anomaly
    sm1 = x1 / a + k
    sm2 = y1 / a + h
    cf = (sm1 * ack - sm2 * adg) / det
    sf = (ach * sm2 - adg * sm1) / det
    fle = np.arctan2(sf, cf)
    lam = fle - k * sf + h * cf                     # Mean longitude

    e = np.sqrt(k ** 2 + h ** 2)                    # Ellipticity
    i = 2 * np.arcsin(np.sqrt(p ** 2 + q ** 2))     # Inclination
    omega = np.arctan2(p, q)                        # Longitude of ascending node
    varpi = np.arctan2(h, k)                        # argument of periapsis

    ell[0] = a
    ell[1] = e
    ell[2] = i
    ell[3] = omega
    ell[4] = varpi
    ell[5] = lam

    return ell                 # Return the orbital elements


def eletocar(ell, GM):
    # Initialize position and velocity vectors
    pos = np.zeros(3)
    vit = np.zeros(3)

    # orbital elements
    a = ell[0]     # Semi-major axis
    k = ell[1]     # e*cos(omega)
    h = ell[2]     # e*sin(omega)
    q = ell[3]     # sin(i/2)*cos(Omega)
    p = ell[4]     # sin(i/2)*sin(Omega)
    lam = ell[5]   # Mean longitude

    # Mean motion
    n = np.sqrt(GM / (a * a * a))
    
    # Compute the eccentric anomaly using iterative method
    fle = lam - k * np.sin(lam) + h * np.cos(lam)
    corf = (lam - fle + k * np.sin(fle) - h * np.cos(fle)) / (1.0 - k * np.cos(fle) - h * np.sin(fle))
    shift = 1.0e-8
    
    while abs(corf) > shift:
        fle += corf
        corf = (lam - fle + k * np.sin(fle) - h * np.cos(fle)) / (1.0 - k * np.cos(fle) - h * np.sin(fle))
        shift *= 1.1
    
    # true anomaly
    lf = -k * np.sin(fle) + h * np.cos(fle)
    sam1 = -k * np.cos(fle) - h * np.sin(fle)
    asr = 1.0 / (1.0 + sam1)
    phi = np.sqrt(1.0 - k * k - h * h)
    psi = 1.0 / (1.0 + phi)
    
    # position in the orbital plane
    x1 = a * (np.cos(fle) - k - psi * h * lf)
    y1 = a * (np.sin(fle) - h + psi * k * lf)
    
    # velocity in the orbital plane
    vx1 = n * asr * a * (-np.sin(fle) - psi * h * sam1)
    vy1 = n * asr * a * (np.cos(fle) + psi * k * sam1)
    
    # Convert to 3D Cartesian coordinates
    cis2 = 2.0 * np.sqrt(1.0 - p * p - q * q)
    tp = 1.0 - 2.0 * p * p
    tq = 1.0 - 2.0 * q * q
    dg = 2.0 * p * q
    
    pos[0] = x1 * tp + y1 * dg
    pos[1] = x1 * dg + y1 * tq
    pos[2] = (-x1 * p + y1 * q) * cis2
    
    vit[0] = vx1 * tp + vy1 * dg
    vit[1] = vx1 * dg + vy1 * tq
    vit[2] = (-vx1 * p + vy1 * q) * cis2

    return pos, vit  # Return the position and velocity vectors