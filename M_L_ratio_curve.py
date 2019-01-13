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

from astropy import coordinates as coord
from astropy import units as unit
import subprocess
from scipy  import interpolate

#################################################################

# Returns mass in solar unit, given the absolute K-band luminosity
def Mass(L_k):

  L = L_k / 1.E10
  
  if L < 0.0927:
    MtoL = 32.0*(L**-0.5)
  elif L >= 0.0927 and L <= 4.423:
    MtoL = 58.0*(L**-0.25)
  elif L > 4.423:
    MtoL = 32*(L**0.15)
  
  
  Mass_out = L_k * MtoL
  
  return Mass_out

#################################################################
def Mass_lst(L_lst):
  
  MASS = []
  for L in L_lst:
    MASS.append(Mass(L))
  
  return MASS

#################################################################


if __name__ == '__main__':
  
  fig = plt.figure(figsize=(6.428, 6), dpi=100)
  ax = fig.add_axes([0.13, 0.1, 0.83,  0.85])
  
  
  L_node = [1.E5, 9.27E8] 
  M_node = Mass_lst(L_node)
  ax.plot(L_node, M_node, '-', color='blue') 
  
  
  L_node = [4.423E10, 1.E13] 
  M_node = Mass_lst(L_node)
  ax.plot(L_node, M_node, '-', color='blue')   


  L_node = [9.27E8, 4.423E10] 
  M_node = Mass_lst(L_node)
  ax.plot(L_node, M_node, ':', color='black')  
  
  ax.plot([1.E5,1.E13], [1.E7,1E15], ':', color='black')
  
  
  
  
  table = np.genfromtxt('M_L_curve.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  L   = 10**table['log10_L']
  M   = 10**table['log10_M']
  ax.plot(L, M, '-', color='red')
  
  f2 = interpolate.interp1d(L, M)
  x = [1.E9, 2.E9, 3.E9, 5.E9, 7.E9, 8.E9, 8.E9, 1.E10, 2.E10, 3.E10]
  ax.plot(x, f2(x), 'o', color='green')
  
  
 


  ax.set_ylabel('Mass ['+r'$M_o$'+']', fontsize=14)
  ax.set_xlabel('Lumonosity ['+r'$L_o$'+']', fontsize=14)
    
    
  plt.xscale('log')
  plt.yscale('log')
  
  
  
  plt.show()
  
  
  
  
  
  
  