

## Importing Important Python Libraries
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



def median_clip(B_lst, K_lst, Ty_min, Ty_max):

  delta = B_lst-K_lst
  N_galaxies = len(delta)
  q1 = np.zeros((N_galaxies,), dtype=np.int)
  q2 = np.zeros((N_galaxies,), dtype=np.int)
  
  q1[np.where(Ty_lst<=Ty_max)] = 1
  q2[np.where(Ty_lst>=Ty_min)] = 1
  qq = q1+q2
  indices = np.where(qq==2)
  x = delta[indices]
  y = Ty_lst[indices]
  
  sig = np.std(x)
  ave = np.median(x)
  N_galaxies = len(x)
  q1 = np.zeros((N_galaxies,), dtype=np.int)
  q2 = np.zeros((N_galaxies,), dtype=np.int)
  
  q1[np.where(x<=ave+2.*sig)] = 1
  q2[np.where(x>=ave-2.*sig)] = 1
  qq = q1+q2
  indices = np.where(qq==2)
  x = x[indices]
  y = y[indices]
  
  return np.median(x), np.median(y), np.std(x), np.std(y)
  
################################################################# 

if __name__ == '__main__':
  

  
  inFile = 'AllSky.v23.cf3.csv'
  
  table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
  
  id   = table['pgc']
  sgl  = table['sgl']
  sgb  = table['sgb']
  Vls  = table['Vls']
  Vhelio = table['Vhelio']
  Ks   = table['Ks'] 
  dcf2 = table['dcf2'] 
  ed   = table['ed'] 
  dm = table['DM'] 
  edm   = table['eDM'] 
  
  Ty   = table['Ty'] 
  
  gl  = table['gl']
  gb  = table['gb']  
  
  
  ra  =  table['ra']
  dec = table['dec']
  coordinate_src = table['coordinate_src']
  Ty_src  = table['Ty_src']
  B_mag   = table['B_mag']
  Ks_src  = table['Ks_src']
  Vls_src = table['Vls_src']
  objname = table['objname']
  
  
  fig = plt.figure(figsize=(45/7., 6), dpi=100)
  ax = fig.add_axes([0.13, 0.1, 0.83,  0.85])	
  
#########################################################################
  K_lst = []
  B_lst = []
  K_model =[]
  Ty_lst = []
  
  for i in range(len(id)):
    
    if Vls_src[i]=='     2MRS ':
      if Ks[i] > 0 and B_mag[i]>0 and Ty[i]>-10 and Ty[i]<2:
         K_lst.append(Ks[i])
         B_lst.append(B_mag[i])
         Ty_lst.append(Ty[i])
         if Ty[i] < 2:
           K_guess = B_mag[i]-4.10
         elif  Ty[i] >= 2 and Ty[i] <=9:
           K_guess = B_mag[i]-4.60+0.25*Ty[i]
         else:
           K_guess = B_mag[i]-2.35  
         K_model.append(K_guess)
         
  K_lst =  np.asarray(K_lst)
  B_lst =  np.asarray(B_lst)
  K_model = np.asarray(K_model)
  Ty_lst = np.asarray(Ty_lst)
  
  
  #ax.plot(B_lst-K_lst, Ty_lst, '.', color='blue')
  
  for ty in np.arange(-5, 2, 0.15):
    x, y, xerr, yerr = median_clip(B_lst, K_lst, ty, ty+0.15)
    #ax.plot([x], [y], '.', color='blue')
    #ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', color='blue')
  
  
  #ax.plot([4.1,4.1], [-6,2], '-', color='black')
  ax.plot(K_lst, K_model, '.', color='blue')
#########################################################################
#########################################################################
  K_lst = []
  B_lst = []
  K_model =[]
  Ty_lst = []
  
  for i in range(len(id)):
    
    if Vls_src[i]=='     2MRS ':
      if Ks[i] > 0 and B_mag[i]>0 and Ty[i]>=2 and Ty[i]<=9:
         K_lst.append(Ks[i])
         B_lst.append(B_mag[i])
         Ty_lst.append(Ty[i])
         if Ty[i] < 2:
           K_guess = B_mag[i]-4.10
         elif  Ty[i] >= 2 and Ty[i] <=9:
           K_guess = B_mag[i]-4.60+0.25*Ty[i]
         else:
           K_guess = B_mag[i]-2.35  
         K_model.append(K_guess)
         
  K_lst =  np.asarray(K_lst)
  B_lst =  np.asarray(B_lst)
  K_model = np.asarray(K_model)
  Ty_lst = np.asarray(Ty_lst)
  
  
  #ax.plot(B_lst-K_lst, Ty_lst, '.', color='red')
  
  
  for ty in np.arange(2, 9, 0.15):
    x, y, xerr, yerr = median_clip(B_lst, K_lst, ty, ty+0.15)
    #ax.plot([x], [y], '.', color='red')  
    #ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', color='red')
  
  #ax.plot([4.1,2.35], [2,9], '-', color='black')
  ax.plot(K_lst, K_model, '.', color='red')
######################################################################### 
#########################################################################
  K_lst = []
  B_lst = []
  K_model =[]
  Ty_lst = []
  
  for i in range(len(id)):
    
    if Vls_src[i]=='     2MRS ':
      if Ks[i] > 0 and B_mag[i]>0 and Ty[i]>9:
         K_lst.append(Ks[i])
         B_lst.append(B_mag[i])
         Ty_lst.append(Ty[i])
         if Ty[i] < 2:
           K_guess = B_mag[i]-4.10
         elif  Ty[i] >= 2 and Ty[i] <=9:
           K_guess = B_mag[i]-4.60+0.25*Ty[i]
         else:
           K_guess = B_mag[i]-2.35  
         K_model.append(K_guess)
         
  K_lst =  np.asarray(K_lst)
  B_lst =  np.asarray(B_lst)
  K_model = np.asarray(K_model)
  Ty_lst = np.asarray(Ty_lst)
  
  
  #ax.plot(B_lst-K_lst, Ty_lst, '.', color='green')
  
  for ty in np.arange(9, 10, 0.15):
    x, y, xerr, yerr = median_clip(B_lst, K_lst, ty, ty+0.15)
    #ax.plot([x], [y], '.', color='green') 
    #ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', color='green')
  
  
  #ax.plot([2.35,2.35], [9,11], '-', color='black')
  ax.plot(K_lst, K_model, '.', color='green')
#########################################################################  
  ax.plot([0,13], [0,13], ':', color='green')
  ax.set_xlabel('Ks ', fontsize=14)
  ax.set_ylabel('Ks (estimated from B-mag)', fontsize=14)
  plt.xlim(0,13)
  plt.ylim(0,13)  
  
  
  
  #ax.set_xlabel('B-Ks ', fontsize=14)
  #ax.set_ylabel('Type', fontsize=14)
  #plt.xlim(-4,8)
  #plt.ylim(-6,11)
  
  
  
  
  plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  