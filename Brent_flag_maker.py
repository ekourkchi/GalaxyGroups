#!/home/ehsan/Ureka/Ureka/variants/common/bin/python


import sys
import os
import random
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from math import *
from time import time
import wl_to_rgb as col
from maxheap import *
import random
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py
from scipy import interpolate
from scipy import polyval, polyfit
from scipy import odr
import copy
import astropy.stats.funcs as st



  
def pgc_gal(file):
  #file = 'north.ML.64.v12.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  pgc = Gtable['pgc']
  flag = Gtable['flag']
  
  N = len(flag)
  filter = np.zeros((N,), dtype=np.int)
  
  f1 = np.zeros((N,), dtype=np.int)
  f1[np.where(flag==1)] = 1  # it's a galaxy in group
  
  f2 = np.zeros((N,), dtype=np.int)
  f2[np.where(flag==0)] = 1 # it's a single galaxy
  
  filter = f1 + f2 
  
  indices = np.where(filter==1)
  pgc = pgc[indices]
  
  indices = np.argsort(pgc)
  pgc = pgc[indices]
  
  return pgc  # it returns the sorted catalog of PGCs
  

def difference(A, B):  # A &  Bmust be both sorted
  
  na = len(A)
  nb = len(B)
  
  M = max(max(A),max(B))
  differ = []
  a = 0 ; b = 0
  for i in range(M+10):
    
    EHSAN = False
    BRENT = False
    if a<na and i == A[a]:
      a+=1
      EHSAN = True
    if b<nb and i == B[b]:
      b+=1
      BRENT = True
    
    if EHSAN==True and BRENT==False:
      differ.append(i)
  
  np.asarray(differ)
  return differ
  
  
if __name__ == '__main__':  
  
  
  my_pgc = pgc_gal('north.ML.64.v12.group')
  brent_pgc = pgc_gal('north.ML.64.v12.group.bb')
  pgc =  difference(my_pgc, brent_pgc)
  myTable = Table()
  myTable.add_column(Column(data=pgc,name='pgc', dtype=np.dtype(int)))
  myTable.write("brent_north_bad_stars.csv", format='ascii.fixed_width',delimiter='|', bookend=False)

  
  
  my_pgc = pgc_gal('south.ML.64.v12.group')
  brent_pgc = pgc_gal('south.ML.64.v12.group.b')
  pgc =  difference(my_pgc, brent_pgc)
  myTable = Table()
  myTable.add_column(Column(data=pgc,name='pgc', dtype=np.dtype(int)))
  myTable.write("brent_south_bad_stars.csv", format='ascii.fixed_width',delimiter='|', bookend=False)
  
  

  
  
  