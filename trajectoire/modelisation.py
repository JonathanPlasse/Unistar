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
from numpy import array, sin, cos, pi, dot, arctan, linspace
from scipy.integrate import odeint

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


MS = 3.08  # Marge de stabilité
Yf, Zf, Lf, Mf, Nf = [0]*5
Vvent, epsilon = [0]*2

def Xf(t):
    return 146.7 if t < 0.97 else 0

def m(t):
    return 1.5-t*0.0869 if t < 0.97 else 1.5-0.0843

def F(Vect, t):
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect

    A = m(t)*D*D**2/2
    B = m(t)*LongueurTube**2/12
    C = m(t)*LongueurTube**2/12
    lx = 0
    ly = MS*D
    lz = MS*D

    Vvent_x = Vvent*cos(epsilon)*cos(psi)*cos(theta)-Vvent*sin(epsilon)*sin(psi)
    Vvent_y = Vvent*cos(epsilon)*(sin(psi)*cos(phi)+sin(phi)*sin(theta))+Vvent*sin(epsilon)*cos(psi)*cos(phi)
    Vvent_z = Vvent*cos(epsilon)*(-sin(psi)*sin(phi)+sin(theta)*cos(phi))-Vvent*sin(epsilon)*cos(psi)*sin(phi)
    Vrx = -u+Vvent_x
    Vry = -v+Vvent_y
    Vrz = -w+Vvent_z

    alpha = -arctan(w/u)
    beta = arctan(v/u*cos(arctan(w/u)))
    rho = rho0*(20000-z)/(20000+z)
    Xa = -rho*Strainee*Vrx**2*Cx/2
    Ya = rho*Sreference*Vry**2*Cyb*beta/2
    Za = -rho*Sreference*Vrz**2*Cza*alpha/2
    La = rho*Strainee*Vrx**2*Cx*lx/2
    Ma = rho*Sreference*Vrz**2*Cza*alpha*lz/2
    Na = rho*Sreference*Vry**2*Cyb*beta*ly/2

    du = 1/m*(Xa+Xf(t)-m(t)*g*sin(theta)+r*v-q*w)
    dv = 1/m*(Ya+Yf+m(t)*g*cos(theta)*sin(phi)+p*w-r*u)
    dw = 1/m*(Za+Zf+m(t)*g*cos(theta)*cos(phi)+q*u-p*v)
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
    t = linspace(0, 30, 1000)
    res = odeint(F, array([0, 0, 0, 0, 1.396, 0, 0, 0, 0, 0, 0, 0]), t)
