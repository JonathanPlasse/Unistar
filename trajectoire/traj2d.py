#!/usr/bin/env
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from math import pi, exp, sqrt, sin, cos, asin
from numpy import array, linspace, sign
from scipy.integrate import odeint

def euler(F, a, b, y0, h):
	"""Solution de y’=F(y,t) sur [a,b], y(a) = y0, pas h"""
	y = y0
	t = a
	les_y = [y0]
	les_t = [a]
	while t+h <= b:
		y = y + h * F(y, t)
		les_y.append(y)
		t += h
		les_t.append(t)
	return les_t, array(les_y)

tmax = 30
n = 1000
M = 0.028966  # en kg.mol^-1
g = 9.81  # en m.s^-2
R = 8.315  # en J.mol^-1.K^-1
T = 30+273.15  # en K
Cx = 0.35  # approximation
hs = R*T/(M*g)
rho0 = 1.292*273.15/T
mFusee = 1.6
mCarbu = 0.0843
tempsComb = 0.97
D = mCarbu/tempsComb
pousse = 146.7
rayon = 0.04
S = pi*rayon**2
z0 = 0
dz0 = 0
x0 = 0
dx0 = 0
theta = 80*pi/180

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


def rho(z, dz, dx):
	return rho0*exp(-z/hs)*S*Cx*v(dz,dx)**2/2


def v(dz, dx):
	return sqrt(dz**2+dx**2)


def F(Y, t):
	global theta
	z, dz, x, dx = Y
	if z < 0:
		ddz = -dz*(n-1)/tmax
		ddx = -dx*(n-1)/tmax
	else:
		ddz = (P(t)-sign(dz)*rho(z, dz, dx))/m(t)*sin(theta)+D*dz-g
		ddx = (P(t)-sign(dx)*rho(z, dz, dx))/m(t)*cos(theta)+D*dx
	if v(dx,dz) != 0:
		theta = asin(dz/v(dx,dz))
	return array([dz, ddz, dx, ddx])


t, z = euler(F, 0, tmax, array([0,0,0,0]), tmax/(n-1))
plt.plot(z[:,2], z[:,0], label='trajectoire')
plt.xlabel('x (en m)')
plt.ylabel('y (en m)')
plt.legend()
plt.show()
