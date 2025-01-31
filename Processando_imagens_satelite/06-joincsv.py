#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 16:38:33 2015

@author: andre
"""
# Joins CSV files with a single header line
# Usage: python ./joincsv.py datadir outfile
# Example: python ./joincsv.py /home/andre output.csv

import os
import sys
import glob
import pandas as pd

datadir=sys.argv[1]
curdir=os.getcwd()
outfile=sys.argv[2]

allFiles = sorted(glob.glob(datadir+'/*.csv'))
frame = pd.DataFrame()
list = []
for file in allFiles:
    str=''
    df = pd.read_csv(file,index_col=None, header=0)
    list.append(df)
frame = pd.concat(list)
csv_out = open(outfile, 'wb')
frame.to_csv(csv_out,index=False)
csv_out.close()

exit()

