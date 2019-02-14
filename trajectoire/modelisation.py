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
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

tMax = 15
mTotal = 1.8
mCombustible = 0.0843
tCombustion = 0.97

def P(t):
    if t < 0:
        return 0
    elif t < .04:
        return 250-50/.04*t
    elif t < .68:
        return 200-50/(.68-.4)*(t-.4)
    elif t < 0.84:
        return 150-100/(.84-.68)*(t-.68)
    elif t < 1.04:
        return 50-50/(1.04-.84)*(t-.84)
    else:
        return 0

g = 9.81
rho0 = 1.013  # La densité de l'air
Cx = 0.3
Cyb = 18.15  # Cy_beta
Cza = 18.15  # Cz_alpha
LongueurTube = 0.75
D = 0.08  # Diamètre de la fusée
Dogive = D  # Diamètre de l'ogive
L = 0.06634  # Envergure des ailerons
EPaileron = 0.002  # Épaisseur des ailerons
Strainee = pi*D**2/4+4*L*EPaileron  # Surface de trainée
Sreference = pi*Dogive**2/4  # Surface de référence
LongueurRampe = 2

MS = 3.08  # Marge de stabilité
Yf, Zf, Lf, Mf, Nf = [0]*5  # Force et moment appliqué à la fusée
Vvent = 10  # Vitesse du vent
epsilon = 0  # Direction du vent

def F(Vect, t):
    """Fonction de dérivation du vecteur d'état Vect en fonction de t"""
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect
    m = mTotal-t*mCombustible/tCombustion if t < tCombustion else mTotal-mCombustible
    A = m*D*D**2/2
    B = m*LongueurTube**2/12
    C = m*LongueurTube**2/12
    ly = MS*D
    lz = MS*D

    Vvent_x = Vvent*cos(epsilon)*cos(psi)*cos(theta)-Vvent*sin(epsilon)*sin(psi)
    Vvent_y = Vvent*cos(epsilon)*(sin(psi)*cos(phi)+sin(phi)*sin(theta))+Vvent*sin(epsilon)*cos(psi)*cos(phi)
    Vvent_z = Vvent*cos(epsilon)*(-sin(psi)*sin(phi)+sin(theta)*cos(phi))-Vvent*sin(epsilon)*cos(psi)*sin(phi)
    Vrx = -u+Vvent_x
    Vry = -v+Vvent_y
    Vrz = -w+Vvent_z

    rho = rho0*(20000+z)/(20000-z)
    Xa = -rho*Strainee*Vrx**2*Cx/2
    # Xf = 146.7 if t < 0.97 else 0
    Xf = P(t)

    du = 1/m*(Xa+Xf-m*g*sin(theta)+r*v-q*w)
    if sqrt(x**2+y**2+z**2) < LongueurRampe:  # Test si la fusée est toujours dans la rampe de lancement
        dv, dw, dp, dq, dr = [0]*5
    else:
        alpha = -arctan(w/u)
        beta = arctan(v/u*cos(arctan(w/u)))
        Ya = -rho*Sreference*Vry**2*Cyb*beta/2
        Za = -rho*Sreference*Vrz**2*Cza*alpha/2
        Ma = rho*Sreference*Vrz**2*Cza*alpha*lz/2
        Na = rho*Sreference*Vry**2*Cyb*beta*ly/2

        dv = 1/m*(Ya+Yf+m*g*cos(theta)*sin(phi)+p*w-r*u)
        dw = 1/m*(Za+Zf+m*g*cos(theta)*cos(phi)+q*u-p*v)
        dp = 1/A*(Lf+(B-C)*q*r)
        dq = 1/B*(Ma+Mf+(C-A)*p*r)
        dr = 1/C*(Na+Nf+(A-B)*p*q)

    return array([*RtoR0(Vect), du, dv, dw, p, q, r, dp, dq, dr])

def RtoR0(Vect):
    """Renvoie les vitesses de la fusée dans le repère terrestre"""
    cphi, ctheta, cpsi = cos(Vect[6:9])
    sphi, stheta, spsi = sin(Vect[6:9])
    T = array([[cpsi*ctheta, -spsi*cphi+cpsi*stheta*sphi, spsi*sphi+cpsi*stheta*cphi],
               [spsi*ctheta, cpsi*cphi+spsi*stheta*sphi, -cpsi*sphi+spsi*stheta*cphi],
               [-stheta, ctheta*sphi, ctheta*cphi]])
    return dot(T, Vect[3:6])


if __name__ == "__main__":
    # Simulation du vol de la fusée avec différentes direction de vent

    ################
    # Affichage 2D #
    ################
    # t = linspace(0, tMax, 100)
    # for epsilon in linspace(0, pi, 10):
    #     res = odeint(F, array([0, 0, 0, 0, 0, 0, 0, 1.396, 0, 0, 0, 0]), t)
    #     plt.plot(res[:, 0], -res[:, 2])
    # plt.title("Altitude as a function of time")
    # plt.xlabel("t (in seconds)")
    # plt.ylabel("z (in meters)")
    # plt.show()

    ################
    # Affichage 3D #
    ################
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    t = linspace(0, tMax, 100)
    for epsilon in linspace(0, pi, 10):
        res = odeint(F, array([0, 0, 0, 0, 0, 0, 0, pi*80/180, 0, 0, 0, 0]), t)
        ax.plot(res[:, 0], res[:, 1], -res[:, 2])
    ax.set_xlabel("x (in meters)")
    ax.set_ylabel("y (in meters)")
    ax.set_zlabel(r"z (in meters)")
    plt.title("Rocket Trajectory")
    plt.show()
