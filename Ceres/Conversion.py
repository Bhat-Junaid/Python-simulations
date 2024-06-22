import numpy as np
from numba import njit
@njit
def cartoelt(pos, vit, GM):
    pos = np.array(pos)
    vit = np.array(vit)

    ell = np.zeros(6)
    rayon = np.sqrt(np.sum(pos ** 2))
    v2 = np.sum(vit ** 2)
    a = GM * rayon / (2 * GM - rayon * v2)       # semi major axis
    gx = pos[1] * vit[2] - pos[2] * vit[1]
    gy = pos[2] * vit[0] - pos[0] * vit[2]
    gz = pos[0] * vit[1] - pos[1] * vit[0]
    gg = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    cis2 = np.sqrt(0.5 * (1 + gz / gg))
    q = -gy / (2 * gg * cis2)                   # sin(I/2)cos(OMEGA)
    p = gx / (2 * gg * cis2)                    # sin(I/2)sin(OMEGA)
    tp = 1 - 2 * p ** 2
    tq = 1 - 2 * q ** 2
    dg = 2 * p * q
    x1 = tp * pos[0] + dg * pos[1] - 2 * p * cis2 * pos[2]
    y1 = dg * pos[0] + tq * pos[1] + 2 * q * cis2 * pos[2]
    vx1 = tp * vit[0] + dg * vit[1] - 2 * p * cis2 * vit[2]
    vy1 = dg * vit[0] + tq * vit[1] + 2 * q * cis2 * vit[2]
    k = gg * vy1 / GM - x1 / rayon                          # ecos(omega + Omega)
    h = -gg * vx1 / GM - y1 / rayon                         #esin(omega + Omega)
    psi = 1 / (1 + np.sqrt(1 - k ** 2 - h ** 2))
    ach = 1 - psi * h ** 2
    ack = 1 - psi * k ** 2
    adg = psi * h * k
    det = ach * ack - adg ** 2
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

    return ell
