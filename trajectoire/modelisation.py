#!/usr/bin/env
# -*- coding: utf-8 -*-

###########
# À faire #
###########

# Approximations
# - Diminution de la masse linéaire
# - Poussé du moteur constante
# - Marge statique constante
# - La matrice d'inertie de la fusée est un cylindre plein avec une masse uniformément répartie
# - Densité de l'air constante

# FUTUR
# - Réalisation de table
#   - Poussé
#   - Masse
#   - Marge statique
# - Intégration des calculs des constantes
#   - CzAlpha, CyBeta
#   - Centre de gravité
#   - Marge statique

import matplotlib.pyplot as plt
from numpy import array, sin, cos, pi, dot, arctan, sqrt, linspace
from scipy.integrate import odeint

def euler(F, a, b, y0, h):
    """Solution de y’=F(y,t) sur [a,b], y(a) = y0, pas h"""
    y = y0
    t = a
    les_y = [y0]
    les_t = [a]
    while t+h <= b:
        y = y + h * F(y, t)
        # print(y)
        # input()
        les_y.append(y)
        t += h
        les_t.append(t)
    return les_t, array(les_y)

g = 9.81
rho0 = 101325
Cx = 0.3
Cyb = 18.15  # Cy_beta
Cza = 18.15  # Cz_alpha
LongueurTube = 0.75
D = 0.08  # Diamètre de la fusée
Dogive = D  # Diamètre de l'ogive
L = 0.06634  # Envergure des ailerons
EPaileron = 0.002  # Épaisseur des ailerons
Strainee = pi*D**2/4+4*L*EPaileron
Sreference = pi*Dogive**2/4
LongueurRampe = 2

MS = 3.08  # Marge de stabilité
Yf, Zf, Lf, Mf, Nf = [0]*5
Vvent, epsilon = [0]*2

def F(Vect, t):
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect
    m = 1.5-t*0.0869 if t < 0.97 else 1.5-0.0843
    A = m*D*D**2/2
    B = m*LongueurTube**2/12
    C = m*LongueurTube**2/12
    lx = 0
    ly = MS*D
    lz = MS*D

    Vvent_x = Vvent*cos(epsilon)*cos(psi)*cos(theta)-Vvent*sin(epsilon)*sin(psi)
    Vvent_y = Vvent*cos(epsilon)*(sin(psi)*cos(phi)+sin(phi)*sin(theta))+Vvent*sin(epsilon)*cos(psi)*cos(phi)
    Vvent_z = Vvent*cos(epsilon)*(-sin(psi)*sin(phi)+sin(theta)*cos(phi))-Vvent*sin(epsilon)*cos(psi)*sin(phi)
    Vrx = -u+Vvent_x
    Vry = -v+Vvent_y
    Vrz = -w+Vvent_z

    rho = rho0*(20000-z)/(20000+z)
    # Xa = -rho*Strainee*Vrx**2*Cx/2
    Xa = 0
    Xf = 146.7 if t < 0.97 else 0

    du = 1/m*(Xa+Xf-m*g*sin(theta)+r*v-q*w)
    if (sqrt(x**2+y**2+z**2) < LongueurRampe):
        dv, dw, dp, dq, dr = [0]*5
    else:
        # alpha = -arctan(w/u)
        # beta = arctan(v/u*cos(arctan(w/u)))
        # Ya = rho*Sreference*Vry**2*Cyb*beta/2
        # Za = -rho*Sreference*Vrz**2*Cza*alpha/2
        # La = rho*Strainee*Vrx**2*Cx*lx/2
        # Ma = rho*Sreference*Vrz**2*Cza*alpha*lz/2
        # Na = rho*Sreference*Vry**2*Cyb*beta*ly/2
        Ya, Za, La, Ma, Na = [0]*5

        dv = 1/m*(Ya+Yf+m*g*cos(theta)*sin(phi)+p*w-r*u)
        dw = 1/m*(Za+Zf+m*g*cos(theta)*cos(phi)+q*u-p*v)
        dp = 1/A*(La+Lf+(B-C)*q*r)
        dq = 1/B*(Ma+Mf+(C-A)*p*r)
        dr = 1/C*(Na+Nf+(A-B)*p*q)

    return array([*RtoR0(Vect), du, dv, dw, p, q, r, dp, dq, dr])

def RtoR0(Vect):
    cphi, ctheta, cpsi = cos(Vect[6:9])
    sphi, stheta, spsi = sin(Vect[6:9])
    T = array([[cpsi*ctheta, -spsi*cphi+cpsi*stheta*sphi, spsi*sphi+cpsi*stheta*cphi],
               [spsi*ctheta, cpsi*cphi+spsi*stheta*sphi, -cpsi*sphi+spsi*stheta*cphi],
               [-stheta, ctheta*sphi, ctheta*cphi]])
    return dot(T, Vect[3:6])


if __name__ == "__main__":
    t = linspace(0, 10, 100)
    res = odeint(F, array([0, 0, 0, 0, 0, 0, 0, 1.396, 0, 0, 0, 0]), t)
    plt.plot(t, res[:, 3])
    plt.show()
