#!/usr/bin/env
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from math import pi, exp
from numpy import array, linspace, sign
from scipy.integrate import odeint

tmax = 20
n = 1000
M = 0.028966  # en kg.mol^-1
g = 9.81  # en m.s^-2
R = 8.315  # en J.mol^-1.K^-1
T = 30+273.15  # en K
Cx = 0.35  # approximation
hs = R*T/(M*g)
rho0 = 1.292*273.15/T
mFusee = 0.9
mCarbu = 0.0843
tempsComb = 0.97
D = mCarbu/tempsComb
pousse = 146.7
rayon = 0.08
S = pi*rayon**2
z0 = 0
dz0 = 0


def m(t):
	if t < tempsComb:
		return mFusee-D*t
	else:
		return mFusee-mCarbu
	return


def P(t):
	if t < 0.8:
		return 220-(220-150)*t/0.8
	elif t < tempsComb:
		return 150-150*(t-0.8)/(tempsComb-0.8)
	else:
		return 0


def F(Y, t):
	z, dz = Y
	if z < 0:
		ddz = -dz*(n-1)/tmax
	else:
		ddz = (P(t)+D*dz-sign(dz)*rho0*exp(-z/hs)*S*Cx*dz**2/2)/m(t)-g
	return array([dz, ddz])


t = linspace(0, tmax, n)
z = odeint(F, array([0, 0]), t)

plt.plot(t, z[:,0], label='h')
plt.xlabel('Temps (en s)')
plt.ylabel('Hauteur (en m)')
plt.legend()
plt.show()
