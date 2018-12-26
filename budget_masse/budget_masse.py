#!/usr/bin/env
# -*- coding: utf-8 -*-

import pandas

# masse en g
# hauteur en m

# lecture du fichier masses.csv
# f est un dictionnaire avec comme clé :
# - masse (en g)
# - hauteur (en m)
f = pandas.read_csv('masses.csv')

cg = 0  # hauteur du centre de gravité (en m)
mTotal = 0  # masse total de la fusée (en g)

# calcul du centre de gravité
for i in range(len(f)):
    cg += f["masse"][i] * f["hauteur"][i]
    mTotal += f["masse"][i]
cg /= mTotal

print(cg)
# print(f)
