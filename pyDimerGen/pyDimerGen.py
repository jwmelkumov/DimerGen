#!/usr/bin/env python3

# MIT License
#
# Copyright (c) 2024 John W. Melkumov
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ pyDimerGen #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Author: John W. Melkumov (jmelkumov@gmail.com).                          #
# Date: 06.01.2024.                                                        #
# Description:                                                             #
# This program takes as input two XYZ files (See Ref. * below.)            #
# with each file containing atom labels and coordinates of one of          #
# two monomers, A and B, as well as the desired center-of-mass to          #
# center-of-mass separation between them in angstroms, and 5 Euler         #
# angles in the ZYZ convention.                                            #
# Note: Although 6 Euler angles are used, (alphaA, betaA, gammaA) and      #
# (alphaB, betaB, gammaB), alphaA is hardcoded and set to 0.               #
# Ref.:                                                                    #
# * www.ccl.net/chemistry/resources/messages/1996/10/21.005-dir/index.html #
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

import os
import sys
import numpy as np

def compute_inertia_tensor(xyz, masses):
    inertia_tensor = np.zeros((3, 3))
    for i in range(len(masses)):
        r = xyz[i]
        mass = masses[i]
        inertia_tensor += mass * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    return inertia_tensor

def euler_rotation_matrix(alpha, beta, gamma):
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)

    # Define rotation matrices
    Rz_alpha = np.array([[np.cos(alpha_rad), -np.sin(alpha_rad), 0],
                         [np.sin(alpha_rad), np.cos(alpha_rad), 0],
                         [0, 0, 1]])

    Ry_beta = np.array([[np.cos(beta_rad), 0, np.sin(beta_rad)],
                        [0, 1, 0],
                        [-np.sin(beta_rad), 0, np.cos(beta_rad)]])

    Rz_gamma = np.array([[np.cos(gamma_rad), -np.sin(gamma_rad), 0],
                         [np.sin(gamma_rad), np.cos(gamma_rad), 0],
                         [0, 0, 1]])

    # Combine the rotation matrices
    return np.dot(Rz_gamma, np.dot(Ry_beta, Rz_alpha))

if len(sys.argv) >= 6:
    xyzfileA = sys.argv[1]  
    xyzfileB = sys.argv[2]
    sep = float(sys.argv[3])  
    alphaA = 0.0
    betaA = float(sys.argv[4])
    gammaA = float(sys.argv[5]) 
    alphaB = float(sys.argv[6]) 
    betaB = float(sys.argv[7])
    gammaB = float(sys.argv[8])
else:
    print("Usage: ./pyDimerGen.py xyzfileA xyzfileB sep betaA gammaA alphaB betaB gammaB")

monAfile = xyzfileA
monBfile = xyzfileB
labelsA = []
xyzA = []
labelsB = []
xyzB = []
with open(monAfile, 'r') as f:
    next(f)
    next(f)
    for line in f:
        if len(line) >= 4:
            col = line.split()
            labelsA.append(col[0]) 
            floats = [float(col[i]) for i in range(1,4)]
            xyzA.append(floats)
with open(monBfile, 'r') as f:
    next(f)
    next(f)
    for line in f:
        if len(line) >= 4:
            col = line.split()
            labelsB.append(col[0]) 
            floats = [float(col[i]) for i in range(1,4)]
            xyzB.append(floats)

massdict = {
'O': 16.000, 
'H': 1.000, 
'C': 12.000, 
'N': 14.000, 
'M': 0.000, # TIP4P Off-Atomic Site
'L': 0.000, # Force Field "Lone Pair"
'He': 4.003,
'Li': 6.940,
'Be': 9.012,
'B': 10.810,
'F': 18.998,
'Ne': 20.180,
'Na': 22.990,
'Mg': 24.305,
'Al': 26.982,
'Si': 28.085,
'P': 30.974,
'S': 32.060,
'Cl': 35.450,
'Ar': 39.948,
'K': 39.098,
'Ca': 40.078,
'Sc': 44.956,
'Ti': 47.867,
'V': 50.942,
'Cr': 51.996,
'Mn': 54.938,
'Fe': 55.845,
'Co': 58.933,
'Ni': 58.693,
'Cu': 63.546,
'Zn': 65.380,
'Ga': 69.723,
'Ge': 72.630,
'As': 74.922,
'Se': 78.971,
'Br': 79.904,
'Kr': 83.798,
'Rb': 85.468,
'Sr': 87.620,
'Y': 88.906,
'Zr': 91.224,
'Nb': 92.906,
'Mo': 95.950,
'Tc': 98.000,
'Ru': 101.070,
'Rh': 102.906,
'Pd': 106.420,
'Ag': 107.868,
'Cd': 112.414,
'In': 114.818,
'Sn': 118.710,
'Sb': 121.760,
'Te': 127.600,
'I': 126.904,
'Xe': 131.293,
'Cs': 132.905,
'Ba': 137.327,
'La': 138.905,
'Ce': 140.116,
'Pr': 140.907,
'Nd': 144.242,
'Pm': 145.000,
'Sm': 150.360,
'Eu': 151.964,
'Gd': 157.250,
'Tb': 158.925,
'Dy': 162.500,
'Ho': 164.930,
'Er': 167.259,
'Tm': 168.934,
'Yb': 173.045,
'Lu': 174.966,
'Hf': 178.490,
'Ta': 180.948,
'W': 183.840,
'Re': 186.207,
'Os': 190.230,
'Ir': 192.217,
'Pt': 195.084,
'Au': 196.967,
'Hg': 200.592,
'Tl': 204.383,
'Pb': 207.200,
'Bi': 208.980,
'Po': 209.000,
'At': 210.000,
'Rn': 222.000,
'Fr': 223.000,
'Ra': 226.000,
'Ac': 227.000,
'Th': 232.038,
'Pa': 231.036,
'U': 238.029,
'Np': 237.000,
'Pu': 244.000,
'Am': 243.000,
'Cm': 247.000,
'Bk': 247.000,
'Cf': 251.000,
'Es': 252.000,
'Fm': 257.000,
'Md': 258.000,
'No': 259.000,
'Lr': 266.000,
'Rf': 267.000,
'Db': 270.000,
'Sg': 271.000,
'Bh': 270.000,
'Hs': 269.000,
'Mt': 278.000,
'Ds': 281.000,
'Rg': 282.000,
'Cn': 285.000,
'Nh': 286.000,
'Fl': 289.000,
'Mc': 290.000,
'Lv': 293.000,
'Ts': 294.000,
'Og': 294.000,
}

massA = []
for i in labelsA:
    if (i[0] == 'O'):
        massA.append(massdict['O'])
    elif (i[0] == 'H'):
        massA.append(massdict['H'])
    elif (i[0] == 'C'):
        massA.append(massdict['C'])
    elif (i[0] == 'N'):
        massA.append(massdict['N'])
    elif (i[0] == 'M'):
        massA.append(massdict['M'])
    elif (i[0] == 'L') and i[1].isupper():
        massA.append(massdict['L'])
    elif (i[0] == 'B'):
        massA.append(massdict['B'])
    elif (i[0] == 'F'):
        massA.append(massdict['F'])
    elif (i[0] == 'P'):
        massA.append(massdict['P'])
    elif (i[0] == 'S'):
        massA.append(massdict['S'])
    elif (i[0] == 'K'):
        massA.append(massdict['K'])
    elif (i[0] == 'V'):
        massA.append(massdict['V'])
    elif (i[0] == 'Y'):
        massA.append(massdict['Y'])
    elif (i[0] == 'I'):
        massA.append(massdict['I'])
    elif (i[0] == 'W'):
        massA.append(massdict['W'])
    elif (i[0] == 'U'):
        massA.append(massdict['U'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'He':
            massA.append(massdict['He'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Li':
            massA.append(massdict['Li'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Be':
            massA.append(massdict['Be'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ne':
            massA.append(massdict['Ne'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Na':
            massA.append(massdict['Na'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Mg':
            massA.append(massdict['Mg'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Al':
            massA.append(massdict['Al'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Si':
            massA.append(massdict['Si'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cl':
            massA.append(massdict['Cl'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ar':
            massA.append(massdict['Ar'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ca':
            massA.append(massdict['Ca'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Sc':
            massA.append(massdict['Sc'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ti':
            massA.append(massdict['Ti'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cr':
            massA.append(massdict['Cr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Mn':
            massA.append(massdict['Mn'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Fe':
            massA.append(massdict['Fe'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Co':
            massA.append(massdict['Co'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ni':
            massA.append(massdict['Ni'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cu':
            massA.append(massdict['Cu'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Zn':
            massA.append(massdict['Zn'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ga':
            massA.append(massdict['Ga'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ge':
            massA.append(massdict['Ge'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'As':
            massA.append(massdict['As'])            
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Se':
            massA.append(massdict['Se'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Br':
            massA.append(massdict['Br'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Kr':
            massA.append(massdict['Kr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Rb':
            massA.append(massdict['Rb'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Sr':
            massA.append(massdict['Sr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Zr':
            massA.append(massdict['Zr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Nb':
            massA.append(massdict['Nb'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Mo':
            massA.append(massdict['Mo'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Tc':
            massA.append(massdict['Tc'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ru':
            massA.append(massdict['Ru'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Rh':
            massA.append(massdict['Rh'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pd':
            massA.append(massdict['Pd'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ag':
            massA.append(massdict['Ag'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cd':
            massA.append(massdict['Cd'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'In':
            massA.append(massdict['In'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Sn':
            massA.append(massdict['Sn'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Sb':
            massA.append(massdict['Sb'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Te':
            massA.append(massdict['Te'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Xe':
            massA.append(massdict['Xe'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cs':
            massA.append(massdict['Cs'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ba':
            massA.append(massdict['Ba'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'La':
            massA.append(massdict['La'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ce':
            massA.append(massdict['Ce'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pr':
            massA.append(massdict['Pr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Nd':
            massA.append(massdict['Nd'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pm':
            massA.append(massdict['Pm'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Sm':
            massA.append(massdict['Sm'])           
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Eu':
            massA.append(massdict['Eu'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Gd':
            massA.append(massdict['Gd'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Tb':
            massA.append(massdict['Tb'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Dy':
            massA.append(massdict['Dy'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ho':
            massA.append(massdict['Ho'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Er':
            massA.append(massdict['Er'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Tm':
            massA.append(massdict['Tm'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Yb':
            massA.append(massdict['Yb'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Lu':
            massA.append(massdict['Lu'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Hf':
            massA.append(massdict['Hf'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ta':
            massA.append(massdict['Ta'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Re':
            massA.append(massdict['Re'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Os':
            massA.append(massdict['Os'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ir':
            massA.append(massdict['Ir'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pt':
            massA.append(massdict['Pt'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Au':
            massA.append(massdict['Au'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Hg':
            massA.append(massdict['Hg'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Tl':
            massA.append(massdict['Tl'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pb':
            massA.append(massdict['Pb'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Bi':
            massA.append(massdict['Bi'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Po':
            massA.append(massdict['Po'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'At':
            massA.append(massdict['At'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Rn':
            massA.append(massdict['Rn'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Fr':
            massA.append(massdict['Fr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ra':
            massA.append(massdict['Ra'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ac':
            massA.append(massdict['Ac'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Th':
            massA.append(massdict['Th'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pa':
            massA.append(massdict['Pa'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Np':
            massA.append(massdict['Np'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Pu':
            massA.append(massdict['Pu'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Am':
            massA.append(massdict['Am'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cm':
            massA.append(massdict['Cm'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Bk':
            massA.append(massdict['Bk'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cf':
            massA.append(massdict['Cf'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Es':
            massA.append(massdict['Es'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Fm':
            massA.append(massdict['Fm'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Md':
            massA.append(massdict['Md'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'No':
            massA.append(massdict['No'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Lr':
            massA.append(massdict['Lr'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Rf':
            massA.append(massdict['Rf'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Db':
            massA.append(massdict['Db'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Sg':
            massA.append(massdict['Sg'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Bh':
            massA.append(massdict['Bh'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Hs':
            massA.append(massdict['Hs'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Mt':
            massA.append(massdict['Mt'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ds':
            massA.append(massdict['Ds'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Rg':
            massA.append(massdict['Rg'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Cn':
            massA.append(massdict['Cn'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Nh':
            massA.append(massdict['Nh'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Fl':
            massA.append(massdict['Fl'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Mc':
            massA.append(massdict['Mc'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Lv':
            massA.append(massdict['Lv'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Ts':
            massA.append(massdict['Ts'])
    elif len(labelsA) > 1 and i[1].islower():
        if i[:2] == 'Og':
            massA.append(massdict['Og'])
massB = []
for i in labelsB:
    if (i[0] == 'O'):
        massB.append(massdict['O'])
    elif (i[0] == 'H'):
        massB.append(massdict['H'])
    elif (i[0] == 'C'):
        massB.append(massdict['C'])
    elif (i[0] == 'N'):
        massB.append(massdict['N'])
    elif (i[0] == 'M'):
        massB.append(massdict['M'])
    elif (i[0] == 'L') and i[1].isupper():
        massB.append(massdict['L'])
    elif (i[0] == 'B'):
        massB.append(massdict['B'])
    elif (i[0] == 'F'):
        massB.append(massdict['F'])
    elif (i[0] == 'P'):
        massB.append(massdict['P'])
    elif (i[0] == 'S'):
        massB.append(massdict['S'])
    elif (i[0] == 'K'):
        massB.append(massdict['K'])
    elif (i[0] == 'V'):
        massB.append(massdict['V'])
    elif (i[0] == 'Y'):
        massB.append(massdict['Y'])
    elif (i[0] == 'I'):
        massB.append(massdict['I'])
    elif (i[0] == 'W'):
        massB.append(massdict['W'])
    elif (i[0] == 'U'):
        massB.append(massdict['U'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'He':
            massB.append(massdict['He'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Li':
            massB.append(massdict['Li'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Be':
            massB.append(massdict['Be'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ne':
            massB.append(massdict['Ne'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Na':
            massB.append(massdict['Na'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Mg':
            massB.append(massdict['Mg'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Al':
            massB.append(massdict['Al'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Si':
            massB.append(massdict['Si'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cl':
            massB.append(massdict['Cl'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ar':
            massB.append(massdict['Ar'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ca':
            massB.append(massdict['Ca'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Sc':
            massB.append(massdict['Sc'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ti':
            massB.append(massdict['Ti'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cr':
            massB.append(massdict['Cr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Mn':
            massB.append(massdict['Mn'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Fe':
            massB.append(massdict['Fe'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Co':
            massB.append(massdict['Co'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ni':
            massB.append(massdict['Ni'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cu':
            massB.append(massdict['Cu'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Zn':
            massB.append(massdict['Zn'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ga':
            massB.append(massdict['Ga'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ge':
            massB.append(massdict['Ge'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'As':
            massB.append(massdict['As'])            
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Se':
            massB.append(massdict['Se'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Br':
            massB.append(massdict['Br'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Kr':
            massB.append(massdict['Kr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Rb':
            massB.append(massdict['Rb'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Sr':
            massB.append(massdict['Sr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Zr':
            massB.append(massdict['Zr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Nb':
            massB.append(massdict['Nb'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Mo':
            massB.append(massdict['Mo'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Tc':
            massB.append(massdict['Tc'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ru':
            massB.append(massdict['Ru'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Rh':
            massB.append(massdict['Rh'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pd':
            massB.append(massdict['Pd'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ag':
            massB.append(massdict['Ag'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cd':
            massB.append(massdict['Cd'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'In':
            massB.append(massdict['In'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Sn':
            massB.append(massdict['Sn'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Sb':
            massB.append(massdict['Sb'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Te':
            massB.append(massdict['Te'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Xe':
            massB.append(massdict['Xe'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cs':
            massB.append(massdict['Cs'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ba':
            massB.append(massdict['Ba'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'La':
            massB.append(massdict['La'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ce':
            massB.append(massdict['Ce'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pr':
            massB.append(massdict['Pr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Nd':
            massB.append(massdict['Nd'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pm':
            massB.append(massdict['Pm'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Sm':
            massB.append(massdict['Sm'])           
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Eu':
            massB.append(massdict['Eu'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Gd':
            massB.append(massdict['Gd'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Tb':
            massB.append(massdict['Tb'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Dy':
            massB.append(massdict['Dy'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ho':
            massB.append(massdict['Ho'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Er':
            massB.append(massdict['Er'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Tm':
            massB.append(massdict['Tm'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Yb':
            massB.append(massdict['Yb'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Lu':
            massB.append(massdict['Lu'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Hf':
            massB.append(massdict['Hf'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ta':
            massB.append(massdict['Ta'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Re':
            massB.append(massdict['Re'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Os':
            massB.append(massdict['Os'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ir':
            massB.append(massdict['Ir'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pt':
            massB.append(massdict['Pt'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Au':
            massB.append(massdict['Au'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Hg':
            massB.append(massdict['Hg'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Tl':
            massB.append(massdict['Tl'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pb':
            massB.append(massdict['Pb'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Bi':
            massB.append(massdict['Bi'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Po':
            massB.append(massdict['Po'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'At':
            massB.append(massdict['At'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Rn':
            massB.append(massdict['Rn'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Fr':
            massB.append(massdict['Fr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ra':
            massB.append(massdict['Ra'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ac':
            massB.append(massdict['Ac'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Th':
            massB.append(massdict['Th'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pa':
            massB.append(massdict['Pa'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Np':
            massB.append(massdict['Np'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Pu':
            massB.append(massdict['Pu'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Am':
            massB.append(massdict['Am'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cm':
            massB.append(massdict['Cm'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Bk':
            massB.append(massdict['Bk'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cf':
            massB.append(massdict['Cf'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Es':
            massB.append(massdict['Es'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Fm':
            massB.append(massdict['Fm'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Md':
            massB.append(massdict['Md'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'No':
            massB.append(massdict['No'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Lr':
            massB.append(massdict['Lr'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Rf':
            massB.append(massdict['Rf'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Db':
            massB.append(massdict['Db'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Sg':
            massB.append(massdict['Sg'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Bh':
            massB.append(massdict['Bh'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Hs':
            massB.append(massdict['Hs'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Mt':
            massB.append(massdict['Mt'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ds':
            massB.append(massdict['Ds'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Rg':
            massB.append(massdict['Rg'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Cn':
            massB.append(massdict['Cn'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Nh':
            massB.append(massdict['Nh'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Fl':
            massB.append(massdict['Fl'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Mc':
            massB.append(massdict['Mc'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Lv':
            massB.append(massdict['Lv'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Ts':
            massB.append(massdict['Ts'])
    elif len(labelsB) > 1 and i[1].islower():
        if i[:2] == 'Og':
            massB.append(massdict['Og'])

xyzmatA = np.array(xyzA)
xyzmatB = np.array(xyzB)

# x
xA = 0.0
tmassA = 0.0
for i in range(0, len(xyzmatA)):
    xA += xyzmatA[i,0] * massA[i]
    tmassA += massA[i]
comxA = xA / tmassA
# y
yA = 0.0
tmassA = 0.0
for i in range(0, len(xyzmatA)):
    yA += xyzmatA[i,1] * massA[i]
    tmassA += massA[i]
comyA = yA / tmassA

# z
zA = 0.0
tmassA = 0.0
for i in range(0, len(xyzmatA)):
    zA += xyzmatA[i,2] * massA[i]
    tmassA += massA[i]
comzA = zA / tmassA

comA = np.array([comxA, comyA, comzA])

# Get COM coordinates for monomer B:
# x
xB = 0.0
tmassB = 0.0
for i in range(0, len(xyzmatB)):
    xB += xyzmatB[i,0] * massB[i]
    tmassB += massB[i]
comxB = xB / tmassB
# y
yB = 0.0
tmassB = 0.0
for i in range(0, len(xyzmatB)):
    yB += xyzmatB[i,1] * massB[i]
    tmassB += massB[i]
comyB = yB / tmassB

# z
zB = 0.0
tmassB = 0.0
for i in range(0, len(xyzmatB)):
    zB += xyzmatB[i,2] * massB[i]
    tmassB += massB[i]
comzB = zB / tmassB

comB = np.array([comxB, comyB, comzB])

# translate monomer A to COM:
centered_monomer_A = xyzmatA - comA

# translate monomer B to COM:
centered_monomer_B = xyzmatB - comB

# compute inertia tensor
inertia_tensor_A = compute_inertia_tensor(centered_monomer_A, massA)
inertia_tensor_B = compute_inertia_tensor(centered_monomer_B, massB)

# compute eigenvalues and eigenvectors
# (eigenvalues are sorted in ascending order)
eigenvalues_A, eigenvectors_A = np.linalg.eigh(inertia_tensor_A)
eigenvalues_B, eigenvectors_B = np.linalg.eigh(inertia_tensor_B)

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
#$$$$$$$$$$$$$$$$$$$$$$$$$$ ROTATE MONOMER A #$$$$$$$$$$$$$$$$$$$$$$$$$$
principal_axes_A = eigenvectors_A

# transform monomer A coordinates to principal axes frame
centered_monomer_A_inPA = np.dot(centered_monomer_A, principal_axes_A.T)

# construct rotation matrix for monomer A
rotation_matA = euler_rotation_matrix(alphaA, betaA, gammaA)

# apply the rotation matrix in the principal axes frame
centered_monomer_A_rotated_inPA = np.dot(centered_monomer_A_inPA, rotation_matA)

# transform back to original frame
centered_monA_rotated = np.dot(centered_monomer_A_rotated_inPA, principal_axes_A)
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
#$$$$$$$$$$$$$$$$$$$$$$$$$$ ROTATE MONOMER B #$$$$$$$$$$$$$$$$$$$$$$$$$$
principal_axes_B = eigenvectors_B

# transform monomer B coordinates to principal axes frame
centered_monomer_B_inPA = np.dot(centered_monomer_B, principal_axes_B.T)

rotation_matB = euler_rotation_matrix(alphaB, betaB, gammaB)

# apply the rotation matrix in the principal axes frame
centered_monomer_B_rotated_inPA = np.dot(centered_monomer_B_inPA, rotation_matB)

# transform back to original frame
centered_monB_rotated = np.dot(centered_monomer_B_rotated_inPA, principal_axes_B)
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

# get vector from COM of monomer A to COM of monomer B
RAB = comB - comA
# normalize RAB
magRAB = np.linalg.norm(RAB)
RABuvec = RAB / magRAB

# translate monomer B away from COM of monomer A
translation = sep * RABuvec
trans_rot_monomer_B = centered_monB_rotated + translation

namesA = np.array(labelsA)
namesB = np.array(labelsB)
fmonA = np.array([[f"{num:.6f}" for num in row] for row in centered_monA_rotated.astype(float)])
fmonB = np.array([[f"{num:.6f}" for num in row] for row in trans_rot_monomer_B.astype(float)])
monomerA = np.insert(fmonA.astype(str), 0, namesA, axis=1)
monomerB = np.insert(fmonB.astype(str), 0, namesB, axis=1)

print(len(monomerA)+len(monomerB))
print(f"Dimer_{sep}_Å_COM-COM_separation_{alphaA}deg_{betaA}deg_{gammaA}deg_{alphaB}deg_{betaB}deg_{gammaB}deg")

for i in range(0, len(monomerA)):
    print(monomerA[i,0], 
          monomerA[i,1], 
          monomerA[i,2], 
          monomerA[i,3])
for i in range(0, len(monomerB)):
    print(monomerB[i,0], 
          monomerB[i,1], 
          monomerB[i,2], 
          monomerB[i,3])

