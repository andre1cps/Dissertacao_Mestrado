#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 17:03:14 2017

@author: andre
"""

import os
import sys
import pandas as pd
import glob

datadir=sys.argv[1]
curdir=os.getcwd()
outfile=sys.argv[2]

slash=''
if datadir[len(datadir)-1]!='/':
    slash='/'

flist = sorted(glob.glob(datadir+slash+'*reffready.csv'))
frame = pd.DataFrame()
list = []
for file in flist:
    str=''
    df = pd.read_csv(file,index_col=None, header=0)
    list.append(df)
frame = pd.concat(list)
csv_out = open(outfile, 'wb')
frame.to_csv(csv_out,index=False)
csv_out.close()

exit()
