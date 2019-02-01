#!/usr/bin/env
# -*- coding: utf-8 -*-

###########
# À faire #
###########
# Rentrer table de poussé de la fusée
# Calcul de Xf, Yf, Zf, Lf, Mf, Nf en fonction du temps
# Calcul de la masse en fonction du temps
# Calcul du centre de gravité en fonction du temps
#

# import matplotlib.pyplot as plt
from numpy import array, sin, cos, pi, dot, atan
import sys
# from scipy.integrate import odeint

u0, v0, w0, phi0, theta0, psi0, p0, q0, r0 = [0]*9

g = 9.81
rho0 = 1015
Cx = 0
Cyb = 0 # Cy_beta
Cza = 0 # Cz_alpha
LongueurTube = 0
D = 0 # Diamètre de la fusée
Dogive = D # Diamètre de l'ogive
L = 0 # Envergure des ailerons
EPaileron = 0 # Épaisseur des ailerons
lx = 0
ly = MS*D
lz = MS*D

Xf = lambda t:146.7 if t<0.97 else 0
Yf, Zf, Lf, Mf, Nf, Vvent, epsilon = [0]*7
m = lambda t:1.5-t*0.0869 if t<0.97 else 1.5-0.0843
MS = 0 # Marge de stabilité

A = m*D*D**2/2
B = m*LongueurTube**2/12
C = m*LongueurTube**2/12
Strainee = pi*D**2/4+4*L*EPaileron
Sreference = pi*Dogive**2/4

# Vect et U sont des vecteurs dans le repère R de la fusée
Vect0 = array([u0, v0, w0, phi0, theta0, psi0, p0, q0, r0])
U0 = array([Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0])


def F(Vect, t):
    u, v, w, phi, theta, psi, p, q, r = Vect

    Vvent_x = Vvent*cos(epsilon)*cos(psi)*cos(theta)-Vvent*sin(epsilon)*sin(psi)
    Vvent_y = Vvent*cos(epsilon)*(sin(psi)*cos(phi)+sin(phi)*sin(theta))+Vvent*sin(epsilon)*cos(psi)*cos(phi)
    Vvent_z = Vvent*cos(epsilon)*(-sin(psi)*sin(phi)+sin(theta)*cos(phi))-Vvent*sin(epsilon)*cos(psi)*sin(phi)
    Vrx = -u+Vvent_x
    Vry = -v+Vvent_y
    Vrz = -w+Vvent_z

    alpha = -atan(w/u)
    beta = atan(v/u*cos(atan(w/u)))
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

    return array([du, dv, dw, p, q, r, dp, dq, dr])

def RtoR0(Vect):
    cphi, ctheta, cpsi = cos(Vect[3:6])
    sphi, stheta, spsi = sin(Vect[3:6])
    T = array([[cpsi*ctheta, -spsi*cphi+cpsi*stheta*sphi, spsi*sphi+cpsi*stheta*cphi],
           [spsi*ctheta, cpsi*cphi+spsi*stheta*sphi, -cpsi*sphi+spsi*stheta*cphi],
           [-stheta, ctheta*sphi, ctheta*cphi]])
    return dot(T, Vect[3:6])

if __name__ == "__main__":
    print(RtoR0(array([*map(float, sys.argv[1:7]), p0, q0, r0])))
