#!/usr/bin/env
# -*- coding: utf-8 -*-

# import matplotlib.pyplot as plt
# from math import pi, exp
from numpy import array  # , linspace, sign
# from scipy.integrate import odeint

x0, y0, z0, u0, v0, w0, phi0, theta0, psi0, p0, q0, r0 = 0
Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0 = 0
m, g, Xa, Ya, Za,

# Vect et U sont des vecteurs dans le repère R de la fusée
Vect0 = array([x0, y0, z0, u0, v0, w0, phi0, theta0, psi0, p0, q0, r0])
U0 = array([Xf0, Yf0, Zf0, Lf0, Mf0, Nf0, Vvent0, epsilon0])

def F(Vect, U):
    du=1/m
