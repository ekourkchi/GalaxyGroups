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
  
  L_node = [1.E5, 9.27E8, 4.423E10, 1.E15] 
  M_node = Mass_lst(L_node)
   

  
  L_node1 = np.asarray(L_node[0:2])
  M_node1 = np.asarray(M_node[0:2])
  ax.plot(L_node1, M_node1, '-', color='red', linewidth=2)
  L_node2 = np.asarray(L_node[2:4])
  M_node2 = np.asarray(M_node[2:4])  
  ax.plot(L_node2, M_node2, '-', color='red', linewidth=2)  
  
  
   
  
  L = [1.E10, 1.E7]
  L = np.asarray(L)
  L10 = L / 1.E10
  MtoL = 32*(L10**0.15)
  M = L * MtoL
  print L, L10
  print M
  ax.plot(L, M, ':', color='black', linewidth=1) 


  M = []
  L = 9.27E8
  betta1 = Mass(L)*L**2
  for L in L_node:
    M.append(betta1/L/L)
    
  
  M = []
  L = 1.E10
  betta0 = 32E10*L**(4./3)
  for L in L_node:
    M.append(betta0*L**(-4./3))
  


  M = []
  L = 4.423E10
  betta2 = Mass(L)*L**(20./23)
  for L in L_node:
    M.append(betta2*L**(-20./23))

  # black-red
  L_br = (betta1/betta0)**(3./2)
  M_br = betta1/L_br/L_br


  # black-orange
  L_bo = (betta2/betta0)**(-69./32)
  M_bo = betta2*L_bo**(-20./23)
  
  # red
  L = 9.27E8
  logL = log10(L)
  logM = log10(Mass(L))
  logL_br = log10(L_br)
  logM_br = log10(M_br)
  r_br = sqrt((logL-logL_br)**2+(logM-logM_br)**2)
  
  x = np.arange(log10(9.27E8),10,0.001)
  x_lst =[]
  y_lst=[]
  for xi in x:
    yi = -1.*sqrt(r_br**2-(xi-logL_br)**2)+logM_br
    if 10**yi<betta0*(10**xi)**(-4./3):
      x_lst.append(xi)
      y_lst.append(yi)
      
  x_lst_r = np.asarray(x_lst)    
  y_lst_r = np.asarray(y_lst)
  
  ax.plot(10**x_lst_r, (10**y_lst_r), '-', color='red', linewidth=2)   # 'orange'
  
  # orange
  L = 4.423E10
  logL = log10(L)
  logM = log10(Mass(L))
  logL_bo = log10(L_bo)
  logM_bo = log10(M_bo)
  r_bo = sqrt((logL-logL_bo)**2+(logM-logM_bo)**2)
  
  x = np.arange(9.5,log10(4.423E10), 0.001)
  x_lst =[]
  y_lst=[]
  for xi in x:
    yi = -1.*sqrt(r_bo**2-(xi-logL_bo)**2)+logM_bo
    if 10**yi>=betta0*(10**xi)**(-4./3):
      x_lst.append(xi)
      y_lst.append(yi)
      
  x_lst_o = np.asarray(x_lst)    
  y_lst_o = np.asarray(y_lst)
  
  ax.plot(10**x_lst_o, (10**y_lst_o), '-', color='red', linewidth=2)   # 'orange'
  



  
  
  ax.set_ylabel('MassL ['+r'$M_\odot$'+']', fontsize=16)
  ax.set_xlabel('Luminosity ['+r'$L_\odot$'+']', fontsize=16)
 
  #plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     

  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)

  ax.annotate(r'$L^{-0.5}$', (1.5E8, 3.E2), rotation=0, color='black', size=18)
  ax.annotate(r'$L^{0.15}$', (4.E11, 7.E1), rotation=0, color='black', size=18)
  

  
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(1.E7,1.E13)
  plt.ylim(1.E9,3.E15)

  
  
  
  plt.show()
  
  
  
  
  
  
  
