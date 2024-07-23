#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: andre
"""

# These are satellite calibration constants for GOES-8, -10, -12, -13, and -14
# Usage in a separate code:
#
# from goesconstants import *
#

# VIS channel 1 definitions
RadVISUnits='Radiance  [W/(m2 sr um)]'
RefVISUnits='Reflectance'
X0=29.0

mVIS13=0.6118208                                 
kVIS13=0.001160                                  
C13=1.255  # Valid for GOES-13 for FEB2013       

# IR channel definitions
RadIRUnits='Radiance  [mW/(m2 sr cm-1)]'
TUnits='Brightness Temperature  [K]'
C1=1.191066e-5   #  mW/(m2 sr cm-4)      # C1 e C2 são a constantes para converter radiância (que antes era GVAR 10-bit) em temperatura efetiva ou "de brilho".
C2=1.438833      #  K/(cm-1)             # Essa conversão será feita invertendo a função de Planck. C1 e C2 também estão no mesmo link oferecido abaixo

#GOES 13 IR (Channels 2 and 4)
m213=227.3889        ############################################################################# 
m313=38.8383         # Constantes dos canais 2, 3, 4 e 6 para converter GVAR 10-bit em           #
m413=5.2285          # radiâncias no infravermelho. (Passo 1)                                    #
m613=5.5297          #                                                                           #
bb213=68.2167        #                  Estão disponíveis na tabela 1.2 de                       #
bb313=29.1287        # http://www.ospo.noaa.gov/Operations/GOES/calibration/gvar-conversion.html #
bb413=15.6854        #                                                                           #
bb613=16.5892        #############################################################################
n213=2561.7421       #############################################################################
n313=1522.5182       #
n413=937.23449       # Constantes dos canais 2, 3, 4 e 6 para converter radiância para 
n613=749.82589       # temperatura de brilho (ou "Teff", Passo 2) e depois de Teff para a temp. 
a213=-1.4755462      # real em K (Passo 3 ou 4, que é o 3 melhorado). Os valores de n, a, b e g 
a313=-4.1556932      # estão na Tabela 3-6 (Side 1), detector "a", do mesmo endereço web acima ->   
a413=-0.52227011     # 
a613=-0.16089410     # http://www.ospo.noaa.gov/Operations/GOES/calibration/tables/table3_6.htm
b213=1.0028656       #
b313=1.0142082       #
b413=1.0023802       #
b613=1.0006896       #
g213=-5.8203946e-7   #
g313=-8.0255086e-6   #
g413=-2.0798856e-6   #
g613=-3.9853774e-7   #############################################################################
factor13=232.4/354.1531901968   # RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]
F0_13= 9.74384716655927  # Solar constant averaged over GOES-13 Ch2 (W/m2) (Platnick and Fontenla, 2008)
Bfunction13=[3.22252E-004, 2.29538E-049, 1.99043E+001, -4.15215E-006, 1.32596E-008, 8.85164E-001, -3.57508E-002, 5.90154E-004, -5.12428E-006, 2.47901E-008, -6.35867E-011, 6.77820E-014, 3.25625E-007, 5.56472E-052, 4.61646E-004, 2.56534E-009, 1.26165E-011, 4.35725E-006, 3.09653E-008, 1.29105E-010, 5.14528E-013, 1.99084E-015, 7.02852E-018, 2.04987E-020]
satlat13=-0.125907    #Satellite subpoint latitude
satlon13=-74.574043   #Satellite subpoint longitude

