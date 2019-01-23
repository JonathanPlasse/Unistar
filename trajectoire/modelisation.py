#!/usr/bin/env
# -*- coding: utf-8 -*-

# import matplotlib.pyplot as plt
from numpy import array, sin, cos, pi  # , linspace, sign
# from scipy.integrate import odeint

x0, y0, z0, u0, v0, w0, phi0, theta0, psi0, p0, q0, r0 = 0
Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0 = 0
m, g, Cx, Cy, Cz, lx, ly, lz, rho = 0
LongueurTube = 0
D = 0
Dogive = 0
L = 0
EPaileron = 0

A = m*D*D**2/2
B = m*LongueurTube**2/12
C = m*LongueurTube**2/12
Strainee = pi*D**2/4+4*L*EPaileron
Sreference = pi*Dogive**2/4

# Vect et U sont des vecteurs dans le repère R de la fusée
Vect0 = array([x0, y0, z0, u0, v0, w0, phi0, theta0, psi0, p0, q0, r0])
U0 = array([Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0])


def F(Vect, U):
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect
    Xf, Yf, Zf, Lf, Mf, Nf, Vvent, epsilon = U

    Vvent_x = Vvent*cos(epsilon)*cos(psi)*cos(theta)-Vvent*sin(epsilon)*sin(psi)
    Vvent_y = Vvent*cos(epsilon)*(sin(psi)*cos(phi)+sin(phi)*sin(theta))+Vvent*sin(epsilon)*cos(psi)*cos(phi)
    Vvent_z = Vvent*cos(epsilon)*(-sin(psi)*sin(phi)+sin(theta)*cos(phi))-Vvent*sin(epsilon)*cos(psi)*sin(phi)
    Vrx = -u+Vvent_x
    Vry = -v+Vvent_y
    Vrz = -w+Vvent_z

    Xa = -rho*Strainee*Vrx**2*Cx/2
    Ya = rho*Strainee*Vry**2*Cy/2
    Za = -rho*Strainee*Vrz**2*Cz/2
    La = rho*Strainee*Vrx**2*Cx*lx/2
    Ma = rho*Strainee*Vrz**2*Cz*lz/2
    Na = rho*Strainee*Vry**2*Cy*ly/2

    du = 1/m*(Xa+Xf-m*g*sin(theta)+r*v-q*w)
    dv = 1/m*(Ya+Yf+m*g*cos(theta)*sin(phi)+p*w-r*u)
    dw = 1/m*(Za+Zf+m*g*cos(theta)*cos(phi)+q*u-p*v)
    dp = 1/A*(La+Lf+(B-C)*q*r)
    dq = 1/B*(Ma+Mf+(C-A)*p*r)
    dr = 1/C*(Na+Nf+(A-B)*p*q)

    return array([u, v, w, du, dv, dw, p, q, r, dp, dq, dr])
