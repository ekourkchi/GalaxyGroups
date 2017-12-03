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
import random
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py
from scipy import stats
from scipy import odr
from scipy.optimize import curve_fit

import astropy.stats.funcs as st



# **************************************
#    Global Variables
# **************************************
H0 = 75.           # hubble constant
sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude
G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.08567758E19
#t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]

# This is the val;use we are gonna adopt for the paper
t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr
# **************************************

################################################################# 

if __name__ == '__main__':
  

  file = 'north.ML.64.v12.group'
  #file = 'south.ML.64.v12.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  flag = Gtable['flag']
  Vls  = Gtable['Vls']
  noGal = Gtable['No_Galaxies']
  Dist = Gtable['mDist']
  dcf2 = Gtable['dcf2'] 
    
  N = len(flag)
  f = np.zeros((N,), dtype=np.int)
  f1 = np.zeros((N,), dtype=np.int)
  f2 = np.zeros((N,), dtype=np.int)
  f3 = np.zeros((N,), dtype=np.int)
  
  
  f1[np.where(flag==1)] = 1 # galaxies in groups
  f2[np.where(flag==0)] = 1 # individual galaxies
  f3[np.where(flag==2)] = 1 # group header
  
  d = np.zeros((N,), dtype=np.int)
  d[np.where(dcf2>0)] = 1 # gal with dcf2

  v_l = np.zeros((N,), dtype=np.int)
  v_l[np.where(Vls > 0)] = 1 
  v_u = np.zeros((N,), dtype=np.int)
  v_u[np.where(Vls < 4000)] = 1 
  
  
  f = f1 + v_l + v_u
  fd = f + d
  indices = np.where(f==3)
  indices_d = np.where(fd==4)
  print "All Gal in the Groups: ", len(flag[indices])  , len(flag[indices_d])  
  
  
  
  f = f2 + v_l + v_u
  indices = np.where(f==3)
  indices_d = np.where(fd==4)
  print "All Gal NOT in the Groups: ", len(flag[indices])   , len(flag[indices_d])    
  
  
  v_l = np.zeros((N,), dtype=np.int)
  v_u = np.zeros((N,), dtype=np.int)
  v_l[np.where(Vls > 100)] = 1 
  v_u[np.where(Vls < 3500)] = 1
  
  haveDist = np.zeros((N,), dtype=np.int)
  haveDist[np.where(Dist != 0)] = 1
  
  f = f3 + v_l + v_u
  f_dist = f + haveDist
  indices = np.where(f==3)
  indices_d = np.where(f_dist==4)
  print "All Groups: ", len(flag[indices])
  print "All Groups with associated distance: ", len(flag[indices_d]), np.sum(noGal[indices_d])
  noGal = noGal[indices]
  print "# of gal in these groups: ", np.sum(noGal)
  print "No. of groups with N>4: ", len(noGal[np.where(noGal>4)]), np.sum(noGal[np.where(noGal>4)])
  
  

  



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

