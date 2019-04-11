#!/usr/bin/env python3
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
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import yaml


with open('caracteristique.yml', 'r') as f:
    c = yaml.load(f.read())

Xcg_sans = c['Xcg_sans']
Xcg_avec = c['Xcg_avec']

# valeurs connues de la fusee
l = c['l']
dref = c['D']
dail = c['D']
m = c['m']
n = c['n']
e = c['e']
p = c['p']
Q = c['nb_ailerons']
Xail = c['Xail']
taille = 1000

# valeurs qui seront calculees plus tard
Cna = 0
Xcpa = 0

CnaAil = 0
CnaCoiffe = 2

XcpaAil = 0
XcpaCoiffe = 0

MSCna = 0

# calculs stabilite
finesse = taille/dref

f = np.sqrt(e**2+(p+(n-m)/2)**2)
CnaAil = (1+dail/(2*e+dail))*(4*Q*(e/dref)**2)/(1+np.sqrt(1+(2*f/(m+n))**2))
XcpaAil = Xail + p*(m+2*n)/(3*(m+n))+1/6*(m+n-m*n/(m+n))
XcpaCoiffe = 7/15*l

# calculs de X CPA et du Cn_alpha
Cna = 2 + CnaAil
Xcpa = (XcpaCoiffe*CnaCoiffe + XcpaAil*CnaAil)/(CnaCoiffe + CnaAil)


def MS(Xcg):
    return abs(Xcg_sans-Xcpa)/dref


MSCna = MS(Xcg_sans)*Cna

# tests stabilite
if finesse < 10 or finesse > 20:
    if taille < 10:
        print("Fusée trop petite")
    else:
        print("fusée trop grande")
elif Cna < 15 or Cna > 30:
    if Cna < 15:
        print("criteres de stabilite non respecte (Cna<15)")
    elif Cna > 30:
        print("criteres de stabilite non respecte (Cna>30)")
elif MS(Xcg_sans) < 1.5 or MS(Xcg_sans) > 6:
    print(f"criteres de stabilite non respecte : MS = {MS(Xcg_sans)}")
elif MSCna < 30 or MSCna > 100:
    print(f"criteres de stabilite non respecte : MSCna = {MSCna}")
else:
    print("fusee stable")
    print(f"finesse : {finesse} ; Cna : {Cna}")
    print(f"MS : {MS(Xcg_sans)} ; MS.Cna : {MSCna}")

tMax = 14
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
Cyb = Cna  # Cy_beta
Cza = Cna  # Cz_alpha
LongueurTube = 0.75
D = 0.08  # Diamètre de la fusée
Dogive = D  # Diamètre de l'ogive
L = 0.06634  # Envergure des ailerons
EPaileron = 0.002  # Épaisseur des ailerons
Strainee = np.pi*D**2/4+4*L*EPaileron  # Surface de trainée
Sreference = np.pi*Dogive**2/4  # Surface de référence
LongueurRampe = 2

Yf, Zf, Lf, Mf, Nf = [0]*5  # Force et moment appliqué à la fusée
Vvent = 10  # Vitesse du vent
epsilon = 0  # Direction du vent


def F(Vect, t):
    """Fonction de dérivation du vecteur d'état Vect en fonction de t"""
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect
    if t < tCombustion:
        m = mTotal - t*mCombustible/tCombustion
    else:
        m = mTotal - mCombustible
    A = m*D*D**2/2
    B = m*LongueurTube**2/12
    C = m*LongueurTube**2/12
    ly = MS(Xcg_sans)*D
    lz = MS(Xcg_sans)*D

    Vvent_x = Vvent*np.cos(epsilon)*np.cos(psi)*np.cos(theta)\
        - Vvent*np.sin(epsilon)*np.sin(psi)
    Vvent_y = Vvent*np.cos(epsilon)\
        * (np.sin(psi)*np.cos(phi) + np.sin(phi)*np.sin(theta))\
        + Vvent*np.sin(epsilon)*np.cos(psi)*np.cos(phi)
    Vvent_z = Vvent*np.cos(epsilon)\
        * (-np.sin(psi)*np.sin(phi)+np.sin(theta)*np.cos(phi))\
        - Vvent*np.sin(epsilon)*np.cos(psi)*np.sin(phi)
    Vrx = -u+Vvent_x
    Vry = -v+Vvent_y
    Vrz = -w+Vvent_z

    rho = rho0*(20000+z)/(20000-z)
    Xa = -rho*Strainee*Vrx**2*Cx/2
    # Xf = 146.7 if t < 0.97 else 0
    Xf = P(t)

    du = 1/m*(Xa+Xf-m*g*np.sin(theta)+r*v-q*w)
    if np.sqrt(x**2+y**2+z**2) < LongueurRampe:  # Test si la fusée est toujours dans la rampe de lancement
        dv, dw, dp, dq, dr = [0]*5
    else:
        alpha = -np.arctan(w/u)
        beta = np.arctan(v/u*np.cos(np.arctan(w/u)))
        Ya = -rho*Sreference*Vry**2*Cyb*beta/2
        Za = -rho*Sreference*Vrz**2*Cza*alpha/2
        Ma = rho*Sreference*Vrz**2*Cza*alpha*lz/2
        Na = rho*Sreference*Vry**2*Cyb*beta*ly/2

        dv = 1/m*(Ya+Yf+m*g*np.cos(theta)*np.sin(phi)+p*w-r*u)
        dw = 1/m*(Za+Zf+m*g*np.cos(theta)*np.cos(phi)+q*u-p*v)
        dp = 1/A*(Lf+(B-C)*q*r)
        dq = 1/B*(Ma+Mf+(C-A)*p*r)
        dr = 1/C*(Na+Nf+(A-B)*p*q)

    return np.array([*RtoR0(Vect), du, dv, dw, p, q, r, dp, dq, dr])

def RtoR0(Vect):
    """Renvoie les vitesses de la fusée dans le repère terrestre"""
    cphi, ctheta, cpsi = np.cos(Vect[6:9])
    sphi, stheta, spsi = np.sin(Vect[6:9])
    T = np.array([[cpsi*ctheta, -spsi*cphi+cpsi*stheta*sphi, spsi*sphi+cpsi*stheta*cphi],
               [spsi*ctheta, cpsi*cphi+spsi*stheta*sphi, -cpsi*sphi+spsi*stheta*cphi],
               [-stheta, ctheta*sphi, ctheta*cphi]])
    return np.dot(T, Vect[3:6])


if __name__ == "__main__":
    # Simulation du vol de la fusée avec différentes direction de vent

    ################
    # Affichage 2D #
    ################
    # t = np.linspace(0, tMax, 100)
    # for epsilon in np.linspace(0, np.pi, 10):
    #     res = odeint(F, np.array([0, 0, 0, 0, 0, 0, 0, 1.396, 0, 0, 0, 0]), t)
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
    t = np.linspace(0, tMax, 100)
    for epsilon in np.linspace(0, np.pi, 10):
        res = odeint(F, np.array([0, 0, 0, 0, 0, 0, 0, np.pi*80/180, 0, 0, 0, 0]), t)
        ax.plot(res[:, 0], res[:, 1], -res[:, 2])
    ax.set_xlabel("x (in meters)")
    ax.set_ylabel("y (in meters)")
    ax.set_zlabel(r"z (in meters)")
    plt.title("Rocket Trajectory")
    plt.show()
