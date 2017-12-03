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

################################################################# 
sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude
# **************************************
#    Global Variables
# **************************************
H0 = 75.           # hubble constant

G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.09E19
t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]
#M_to_L = 43 
#Rho_cr = (Mpc_km*1000)**3/(4*pi*G*t0_age**2)/M_sun/M_to_L

#################################################################
 
def biweight_Rg(ra, dec, v, d):
   
   N = len(ra)
   dist_inverse = [] 
   for i in range(0,N-1):
    for j in range(i+1,N):
      distij = angle(ra[i], dec[i], ra[j], dec[j])
      if distij !=0 : 
         dist_inverse.append(  distij)
   
   dist_inverse = np.asarray(dist_inverse)
   n = len(dist_inverse)
   if n!=0:
     Med=np.median(dist_inverse)
     Rh_1 = Med
     for p in range(0,10):
        Rh_1 = st.biweight_location(dist_inverse, M=Med)
        Med = Rh_1
     
     #Rg_radian =  N*N/(n*st.biweight_location(dist_inverse))
     #Rg_radian =  N*N/dist_inverse.sum()
     Rg_radian = Rh_1#2*N*Rh_1/(N-1)
     #print dist_inverse
     #print N, '  ' ,st.biweight_location(dist_inverse), '  ', dist_inverse.sum()/n
   else: 
     Rg_radian = float('nan')
   
   if d == 0 : d = v / H0
   return  d * Rg_radian
  
  

#################################################################
# **************************************
def angle(l1, b1, l2, b2):
  
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)
   
   x1 = cl1 * cb1
   y1 = sl1 * cb1
   z1 = sb1
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = cl2 * cb2
   y2 = sl2 * cb2
   z2 = sb2   
   
   XdotY = x1*x2 + y1*y2 + z1*z2
   #X2 = sqrt(x1**2 + y1**2 + z1**2)
   #Y2 = sqrt(x2**2 + y2**2 + z2**2)
   
   if XdotY > 1 :
     theta12 = 0.
   elif XdotY < -1 :
     theta12 = -1.*pi
   else:
     theta12 = acos(XdotY)  
   return theta12 # radian
# **************************************


################################################################# 

if __name__ == '__main__':
  
  
  fig = plt.figure(figsize=(7.5, 7), dpi=100)
  ax = fig.add_axes([0.13, 0.1, 0.83,  0.85])
  
# **************************************
  #################################################################
  file = 'north.ML.64.v12.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  flag = Gtable['flag']
  No = Gtable['No_Galaxies']
  R2t_dyn = Gtable['R2t_dyn']
  R2t_lum = Gtable['R2t_lum']
  sigmaP_dyn = Gtable['sigmaP_dyn']
  sigmaP_lum = Gtable['sigmaP_lum']
  mDist = Gtable['mDist']
  Vls = Gtable['Vls']
  gl = Gtable['gl']
  gb = Gtable['gb']
  #################################################################
  file = 'south.ML.64.v12.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  flag0 = Gtable['flag']
  No0 = Gtable['No_Galaxies']
  R2t_dyn0 = Gtable['R2t_dyn']
  R2t_lum0 = Gtable['R2t_lum']
  sigmaP_dyn0 = Gtable['sigmaP_dyn']
  sigmaP_lum0 = Gtable['sigmaP_lum']
  mDist0 = Gtable['mDist']
  Vls0 = Gtable['Vls']
  gl0 = Gtable['gl']
  gb0 = Gtable['gb'] 
  
  flag = np.concatenate((flag,flag0))
  No = np.concatenate((No,No0))
  R2t_dyn = np.concatenate((R2t_dyn,R2t_dyn0))
  R2t_lum = np.concatenate((R2t_lum,R2t_lum0))
  sigmaP_dyn = np.concatenate((sigmaP_dyn,sigmaP_dyn0))
  sigmaP_lum = np.concatenate((sigmaP_lum,sigmaP_lum0))
  Vls = np.concatenate((Vls,Vls0))
  mDist = np.concatenate((mDist,mDist0))
  gl = np.concatenate((gl,gl0))
  gb = np.concatenate((gb,gb0))
#################################################################
  Bi_Weight = True
  if Bi_Weight:
    i = 0
    while i < len(flag):
      if flag[i] == 2:
	j = i 
	v_galaxies = []
	gl_galaxies = []
	gb_galaxies = []
	i+=1
	while flag[i] == 1:
	  v_galaxies.append(Vls[i])
	  gl_galaxies.append(gl[i])
	  gb_galaxies.append(gb[i])
	  i+=1
	mean = st.biweight_location(v_galaxies)  # bi-weight mean radial velocity
	sigma = st.biweight_midvariance(v_galaxies) # bi-weight radial velocity Dispersion
	if isnan(mean):  mean = Vls[j]  
	if not isnan(sigma): sigmaP_dyn[j]  = sigma
	Rg_bi = biweight_Rg(gl_galaxies, gb_galaxies, mean, mDist[j])
	#if not isnan(Rg_bi): R2t_dyn[j] = Rg_bi * (pi*0.5/1.05/sqrt(1.5))
	i-=1
      i+=1
#################################################################  
  N = len(flag)
  filter = np.zeros((N,), dtype=np.int)
  
  f1 = np.zeros((N,), dtype=np.int)
  f1[np.where(Vls<100)] = 1
  
  f2 = np.zeros((N,), dtype=np.int)
  f2[np.where(Vls>3500)] = 1 
  
  f3 = np.zeros((N,), dtype=np.int)
  f3[np.where(flag!=2)] = 1
  
  filter = f1 + f2 + f3
  indices = np.where(filter==0)
  
  
  flag = flag[indices]
  No = No[indices]
  R2t_dyn = R2t_dyn[indices]
  R2t_lum = R2t_lum[indices]
  sigmaP_dyn = sigmaP_dyn[indices]
  sigmaP_lum =  sigmaP_lum[indices]
  Vls = Vls[indices]
  #################################################################

  betta = 100/1.8
  
  M1 = []
  Counter = 0 # how many crossign tiems are larger than the age of universe
  for i in range(0, len(flag)):
    Tx = (sqrt(1.5)*R2t_dyn[i]*Mpc_km)/sigmaP_dyn[i]/(sqrt(2.5))/Yr/1.0E9
    
    if Tx > 0 and No[i] >= 5 :
      M1.append(log10(Tx))
      if log10(Tx) > log10(13.8): Counter+=1
      
      

  nbins = np.arange(-0.4,1.8,0.2)
  M1 = np.asarray(M1)
  #n, bins, patches = py.hist( M1, bins = nbins ,  color='grey', alpha=0.5, linewidth=1, normed=0, histtype='stepfilled')
  n, bins, patches = py.hist( 10**M1, bins = 10**nbins ,  color='green', alpha=0.5, label=r'$N \geq 5$', linewidth=3, normed=0, histtype='step')
  print "How many groups (N>=5), with Tx large than the Univ. age:", Counter, len(M1)

  M1 = []
  for i in range(0, len(flag)):
    Tx = (sqrt(1.5)*R2t_dyn[i]*Mpc_km)/sigmaP_dyn[i]/(sqrt(2.5))/Yr/1.0E9
    if Tx > 0 and No[i] >= 20 :
      M1.append(log10(Tx))

  nbins = np.arange(-0.4,1.8,0.2)
  M1 = np.asarray(M1)
  n, bins, patches = py.hist( 10**M1, bins = 10**nbins ,  color='black', alpha=1.0, label=r'$N \geq 20$', linewidth=3, normed=0, histtype='step', linestyle ='dashed')



  
  Tx0 = (sqrt(1.5)*Mpc_km)/368/(sqrt(2.5))/Yr/1.0E9
  plt.plot([Tx0,Tx0],[0,1.8*betta], ':', linewidth=3, color='red')
  ax.annotate(r'$\sigma_p/R_{2t}=368 \/\/ km/s/Mpc$', (10**0.33, 70), rotation=90, color='black', size=14)
  

  print bins
  plt.plot([13.8,13.8],[0,1.8*betta], 'b--', linewidth=2, label=r'$t_0=13.8 \/ Gyr$')
  
  #plt.xlim(-0.4,1.6)
  #plt.ylim(0,100)
  #ax.set_yscale('log')
  ax.set_xscale('log')
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)
  
  ax.annotate('Today', (10**1.025, 1.7*betta), rotation=90, color='blue', size=14)
  ax.annotate('Future', (10**1.35, 1.4*betta), rotation=0, color='black', size=14)
  #ax.annotate('Future', (10**1.3, 1.4), rotation=0, color='black')
  ax.arrow( 10**1.26,1.38*betta, 30, 0 , fc="k", ec="k", head_width=1, head_length=5 )
  # x1, y1, x2-x1, y2-y1

  plt.xlabel('Crossing Time: t'+r'$_x$'+' [Gyr]', fontsize=18)
  plt.ylabel('Number of groups', fontsize=18)
  lg = plt.legend( loc=2, numpoints = 1, prop={'size':15}, labelspacing=0.2, markerscale=3  )
  lg.draw_frame(False)
  
  labels = [0.1,1,10,100]
  ax.set_xticks(labels)
  labels_str = ['0.1','1','10','100']
  ax.set_xticklabels(labels_str)
  
  
  
  print "red line peaks at: ", Tx0, " Gyr"
  print "how many groups with N>=20:  ", len(M1)
  
  plt.show()

  
  