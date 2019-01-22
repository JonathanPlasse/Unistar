#!/usr/bin/env
# -*- coding: utf-8 -*-

# import matplotlib.pyplot as plt
# from math import pi, exp
from numpy import array, sin, cos  # , linspace, sign
# from scipy.integrate import odeint

x0, y0, z0, u0, v0, w0, phi0, theta0, psi0, p0, q0, r0 = 0
Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0 = 0
m, g, Xa, Ya, Za, La, Ma, Na, A, B, C = 0

# Vect et U sont des vecteurs dans le repère R de la fusée
Vect0 = array([x0, y0, z0, u0, v0, w0, phi0, theta0, psi0, p0, q0, r0])
U0 = array([Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0])


def F(Vect, U):
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect
    Xf, Yf, Zf, Lf, Mf, Nf, Vvent, epsilon = U
    du = 1/m*(Xa+Xf-m*g*sin(theta)+r*v-q*w)
    dv = 1/m*(Xa+Yf-m*g*cos(theta)*sin(phi)+p*w-r*u)
    dw = 1/m*(Xa+Zf-m*g*cos(theta)*cos(phi)+q*u-p*v)
    dp = 1/A*(La+Lf+(B-C)*q*r)
    dq = 1/B*(La+Mf+(C-A)*p*r)
    dr = 1/C*(La+Nf+(A-B)*p*q)
    return array([u, v, w, du, dv, dw, p, q, r, dp, dq, dr])
