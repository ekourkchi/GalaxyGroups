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
def Mass(L_k, high=False, low=False):

  L = L_k / 1.E10
  
  if L < 0.0927:
    MtoL = 32.0*(L**-0.5)
  elif L >= 0.0927 and L <= 4.423:
    MtoL = 58.0*(L**-0.25)
  elif L > 4.423:
    MtoL = 32*(L**0.15)
  
  
  if high:
    MtoL = 32*(L**0.15)
  elif low:
    MtoL = 32.0*(L**-0.5)
  
  
  Mass_out = L_k * MtoL
  
  return Mass_out

#################################################################
#################################################################

# Returns mass in solar unit, given the absolute K-band luminosity
def Mass2(L_k):

  L = L_k / 1.E10
  

  if L < 4.423:
    MtoL = 32.0*(L**-0.7)
  elif L >= 4.423:
    MtoL = 32*(L**0.15)
  

  
  
  Mass_out = L_k * MtoL
  
  return Mass_out

#################################################################

def Mass2_lst(L_lst, high=False, low=False):
  
  MASS = []
  for L in L_lst:
    MASS.append(Mass2(L))
  
  return MASS


def Mass_lst(L_lst, high=False, low=False):
  
  MASS = []
  for L in L_lst:
    MASS.append(Mass(L, high=high, low=low))
  
  return MASS

#################################################################


if __name__ == '__main__':
  
  fig = plt.figure(figsize=(45/7., 6), dpi=100)
  #ax = fig.add_axes([0.13, 0.1, 0.83,  0.85])  # m-L relation
  ax = fig.add_axes([0.15, 0.13, 0.80,  0.80])  
  
  L = [1.E5, 1.E7, 5.E7,1.E8, 9.27E8, 5.E9, 1.E10, 4.423E10, 1.E12, 1.E15] 
  M = Mass_lst(L)
  M2 = Mass2_lst(L)
  
  
  L  = np.asarray(L)
  M  = np.asarray(M)
  M2 = np.asarray(M2)
  ax.plot(L, M/L, 'r-o')
  ax.plot(L, M2/L, 'b-.')
  
  
  
  
  
  
  ax.set_ylabel('M'+r'$_v$'+r'$^{exp}$'+'/L'+r'$_{K_s}$'+' ['+r'$M_\odot/L_\odot$'+']', fontsize=16)
  ax.set_xlabel('K'+r'$_s$'+'-band Luminosity ['+r'$L_\odot$'+']', fontsize=16)
 
  #plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     

  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)

  ax.annotate(r'$L^{-0.5}$', (3.E7, 141), rotation=0, color='black', size=18)
  ax.annotate(r'$L^{0.15}$', (4.5E11, 35), rotation=0, color='black', size=18)
  
  ax.annotate(r'$L^{-0.7}$', (2.E8, 800), rotation=0, color='blue', size=18)
  
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(1.E7,1.E13)
  plt.ylim(1.E1,3.E3)

  
  
  
  plt.show()
  
  
  
  
  
  
  
