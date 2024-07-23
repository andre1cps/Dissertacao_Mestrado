#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: andre
"""

#grupo=['CONTROLE']
grupo=['DESMATAMENTO']

# Para o grupo CONTROLE:
if grupo[0] == 'CONTROLE':
    dryseasonbegin = 244   # Sept  1st
    dryseasonend  = 334    # Nov 30st
    cloudVISthreshold  = 0.00  # VIS reflectance above (acima) this limit indicates cloud pixel - 14%, pois não há rios e a refletancia da floresta é abaixo de 14%
    cloudTempthreshold = 403.15  # Temp in Kelvin below (abaixo) this limit indicates cloud pixel - 18°C

# Para o GRUPO DESMATAMENTO:
if grupo[0] == 'DESMATAMENTO':
    dryseasonbegin = 152   # Jun  1st
    dryseasonend  = 273    # Sept 30st
    cloudVISthreshold  = 0.00  # VIS reflectance above (acima) this limit indicates cloud pixel - 20,5%, pois a cidade tem refletancia até 20%
    cloudTempthreshold = 403.15  # Temp in Kelvin below (abaixo) this limit indicates cloud pixel - 18°C

# Define VIS/IR thresholds to detect water, ice or mixed phase pixels
#cloudVISthreshold  = 0.125  # The limit of 0.125 (for the Amazon) works well when comparing to ISCCP data. Previously it was 0.15 
#cloudTempthreshold = 300.0  # The limit of 300 K (for the Amazon) works well when comparing to ISCCP data. Previously it was 283K
icedry   =  7.6176
icewet   =  9.6324
waterdry =  3.9559
waterwet =  4.2353


# Nominal Earth radius in km (usado na função viewzen, dentro da reflectance)
ER = 6371.0

# Nominal Geosynchronous orbit (sattelite radius) in km (usado na função viewzen, dentro da reflectance)
SR = 42164.17478

# Planck's function:
h  = 6.62606896e-34 #   Planck's constant in  J.s
ls = 299792458     #   light speed in  m/s
k  = 1.3806504e-23  #   Boltzmann constant in J/K

t0 = 0.75 # Bidirectional transmission function (at 3.9um) between GOES - cloud tops - GOES, t0=0.75 by Kaufman and Nakajima 1993 (é o t⁰_3,7, 0,7cm) 
F0 = 9.68465 # Solar descending irradiance at TOA for 3.9um, by Platnick and Fontenla 2008, in W/[m2 um]. This is only indicative, the F0 to be used is convoluted with sensor response function (see goesconstants.py)

# Define lower limit for the temperature of liquid supercooled water or hottest ice in K
# -38.1C = 235.05K
tlim = 235.05

# Define percentage of pixels to be considered as 'hot ice', e.g. 0.10 means 10% hottest ice pixels
hoticepercent = 0.050

# Define percentage of pixels to be considered as 'cold water', e.g. 0.10 means 10% coldest water pixels
coldwaterpercent = 0.050

# Define percentage of pixels to be considered as 'hot water', e.g. 0.10 means 10% hottest water pixels
hotwaterpercent = 0.050

# Define local noon in UTC time. Example: in Manaus, time is UTC-4, so LocNoon=160000
LocNoon=160000

# Define 0 deg C in Kelvin
zeroc=273.15

# Definition of NaN (not-a-number) or missing data
mynan=-9999
strnan='-9999'
