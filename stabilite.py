#!/usr/bin/env
# -*- coding: utf-8 -*-

import yaml
from math import sqrt

Xcg_without = 505
Xcg_with = 539

# valeurs connues de la fusee
L = 150
dref = 80
dail = 80
m = 120
n = 100
e = 120
p = 120
Q = 4
Xail = 750
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

f = sqrt(e**2+(p+(n-m)/2)**2)
CnaAil = (1+dail/(2*e+dail))*(4*Q*(e/dref)**2)/(1+sqrt(1+(2*f/(m+n))**2))
XcpaAil = Xail + p*(m+2*n)/(3*(m+n))+1/6*(m+n-m*n/(m+n))
XcpaCoiffe = 7/15*L

# calculs de X CPA et du Cn_alpha
Cna = 2 + CnaAil
Xcpa = (XcpaCoiffe*CnaCoiffe + XcpaAil*CnaAil)/(CnaCoiffe + CnaAil)

MS = abs(Xcg_without-Xcpa)/dref

MSCna = MS*Cna

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
elif MS < 1.5 or MS > 6:
    print("criteres de stabilite non respecte : MS = ", MS)
elif MSCna < 30 or MSCna > 100:
    print("criteres de stabilite non respecte : MSCna = ", MSCna)
else:
    print("fusee stable")
    print("finesse : ", finesse, " ; Cna : ", Cna, " ; MS : ", MS, " ; MS.Cna : ", MSCna)
