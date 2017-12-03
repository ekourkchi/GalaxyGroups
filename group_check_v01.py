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
from scipy import odr, stats
import copy
import astropy.stats.funcs as st
from sklearn import linear_model


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

#################################################################

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x#+ B[1]
#################################################################
#################################################################

def g(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x+ B[1]
#################################################################
################################################################# 

def logbin(x, y, deg=1):
  
  
  x_bin = []
  y_bin = []
  x_err = []
  y_err = []
  
  
  for j in range(-3,1):
    for i in np.arange(1, 10, 1):
        X1 = 1.0*i*(10**j)
        X2 = 1.0*(i+1)*(10**j)
        
        for m in [0]:#range(-1,3):
            for n in [0]:#range(1, 10):  
	      Y1 = 0.#1.0*n*(10**m)
	      Y2 = 10000.#1.0*(n+1)*(10**m)
	      
	      x0 = []
	      y0 = []
	      
	      for p in range(0, len(x)):
		if X1 <= x[p] and x[p] < X2 and y[p] >= Y1 and y[p] < Y2:
		  x0.append(x[p])
		  y0.append(y[p])
	      
	      if len(x0)>=deg:
                x_bin.append(np.average(x0))
                y_bin.append(np.average(y0))
                x_err.append(np.std(x0))
	        y_err.append(np.std(y0))
		

  
  return x_bin, y_bin, x_err, y_err



################################################################# 

def linearlogbin(x, y, x_min=0.001, x_max=2 , deg=1):
  
  
  x_bin = []
  y_bin = []
  x_err = []
  y_err = []
  
  X1 = x_min
  X2 = X1 * 10**0.2
  while(X2<=x_max):

        
        for m in [0]:#range(-1,3):
            for n in [0]:#range(1, 10):  
	      #Y1 = 0.#1.0*n*(10**m)
	      #Y2 = 10000.#1.0*(n+1)*(10**m)
	      
	      x0 = []
	      y0 = []
	      
	      for p in range(0, len(x)):
		if X1 <= x[p] and x[p] < X2: # and y[p] >= Y1 and y[p] < Y2:
		  x0.append(x[p])
		  y0.append(y[p])
	      
	      if len(x0)>=deg:
                x_bin.append(np.average(x0))
                y_bin.append(np.average(y0))
                x_err.append(np.std(x0))
	        y_err.append(np.std(y0))
	X1 = X2
	X2 = X1 * 10**0.2

  
  return x_bin, y_bin, x_err, y_err


#################################################################

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x
#################################################################
# **************************************
# returns radian
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
   return theta12   # radian
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
     #for p in range(0,10):
        #Rh_1 = st.biweight_location(dist_inverse, M=Med)
        #Med = Rh_1
     
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

if __name__ == '__main__':
  
  case =  int(sys.argv[1])
  
  fig = plt.figure(figsize=(7.5, 7), dpi=100)
  ax = fig.add_axes([0.14, 0.11, 0.85,  0.85])
  file = 'all.iter.2.v42.group'
# **************************************

  mytable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  #print mytable.dtype
  
  pgc  = mytable['pgc']
  gl = mytable['gl']
  gb = mytable['gb']
	  
  Vls =  mytable['Vls']
  mDist = mytable['mDist']
  mDistErr = mytable['mDistErr']
  sigmaP_dyn = mytable['sigmaP_dyn']
  dcf2 = mytable['dcf2'] 
  ed   = mytable['ed'] 
  
  N =  mytable['No_Galaxies']
  
  R_2t =  mytable['R2t_dyn'] 
  flag = mytable['flag']
#################################################################
################################################################ 
  X1 = R_2t
  X2 = sigmaP_dyn
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
	if isnan(sigma): sigma = sigmaP_dyn[j]  
	Rg_bi = biweight_Rg(gl_galaxies, gb_galaxies, mean, mDist[j])
	if isnan(Rg_bi): Rg_bi = R_2t[j] /(pi*0.5/1.05/sqrt(1.5))
	X1[j] = Rg_bi * (pi*0.5/1.05/sqrt(1.5))
	X2[j] = sigma
	i-=1
      i+=1
#################################################################
  R_2t = X1
  sigma  = X2
  Rg_bi = R_2t / (pi*0.5/1.05/sqrt(1.5))
  Mv_dyn = (1.E9 * (2.5*pi/G/2.) * (Mpc_km*Rg_bi) * sigma**2)/M_sun  # solar mass


  #################################################################

  p = len(flag)
  filter = np.zeros((p,), dtype=np.int)
  
  f1 = np.zeros((p,), dtype=np.int)
  f1[np.where(Vls<100)] = 1
  
  f2 = np.zeros((p,), dtype=np.int)
  f2[np.where(Vls>3500)] = 1 
  
  f3 = np.zeros((p,), dtype=np.int)
  f3[np.where(flag!=2)] = 1
  
  filter = f1 + f2 + f3
  indices = np.where(filter==0)
  
  
  flag = flag[indices]
  sigma = sigma[indices]
  R_2t = R_2t[indices]
  N = N[indices]
  Vls = Vls[indices]
  pgc = pgc[indices]
  Mv_dyn = Mv_dyn[indices]
  #################################################################

  
  if case == 2: 
    R_2t_cp = np.copy(R_2t)
    R_2t = np.copy(sigma)
    sigma = np.copy(Mv_dyn)
  ##########################################
  
  xx = np.arange(0,2,0.001)
  
  Sigma_lst = []
  R_2t_lst = []
  Mv_dyn_lst = []
  for i in range(0, len(pgc)):
    if N[i] > 0 and  N[i]< 5:
      Sigma_lst.append(sigma[i])
      R_2t_lst.append(R_2t[i])

  #plt.plot(R_2t_lst, Sigma_lst, '.', markersize = 1, color='gray')
  

  
  
  Sigma_lst = []
  R_2t_lst = []
  No = []
  XX = []
  kolors = ['gray', 'blue', 'green', 'orange', 'brown']
  sizes = [3,5,6,7,8]
  for i in range(0, len(pgc)):
    if N[i] >= 5 :
      Sigma_lst.append(sigma[i])
      R_2t_lst.append(R_2t[i])
      No.append(1./N[i]**2)
      XX.append([R_2t[i]])
    
    if N[i]>=5 and N[i]<10: j=0
    if N[i]>=10 and N[i]<20: j=1
    if N[i]>=20 and N[i]<30: j=2
    if N[i]>=30 and N[i]<60: j=3
    if N[i]>=60 : j=4
    if N[i]>=5:
       plt.plot(R_2t[i], sigma[i], 'o', markersize = sizes[j], color=kolors[j], markeredgecolor = kolors[j])
    
  #plt.plot(R_2t_lst, Sigma_lst, 'o', markersize = 2, color='blue', markeredgecolor = 'blue')  
  
###################################################  
  R_2t_lst = np.asarray(R_2t_lst)
  Sigma_lst = np.asarray(Sigma_lst)
  No = np.asarray(No)
  
  
  
  
  if case==2:
      #index = np.where(No>=30)
      x = np.log10(R_2t_lst)
      y = np.log10(Sigma_lst)
      weight = np.zeros(len(x))+1.

      # url: http://docs.scipy.org/doc/scipy/reference/odr.html
      # Orthogonal distance regression
      linear = odr.Model(g)
      mydata = odr.Data(x, y, we=weight)
      myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
      myoutput = myodr.run()
      betta = myoutput.beta
      
      slope = betta[0]
      intercept = betta[1]
      
      line, = plt.plot([1,1000],[10**intercept,(1000.**slope)*(10**intercept)], '-', markersize = 2, color='black') 
      line.set_dashes([6, 3])   
      print " SLOP: ", slope, intercept


  if case==1:
      #index = np.where(No>=30)
      x = R_2t_lst
      y = Sigma_lst
      weight = np.zeros(len(x))+1.
      
      # url: http://docs.scipy.org/doc/scipy/reference/odr.html
      # Orthogonal distance regression
      linear = odr.Model(f)
      mydata = odr.Data(x, y, we=weight)
      myodr = odr.ODR(mydata, linear, beta0=[1.])
      myoutput = myodr.run()
      betta = myoutput.beta
      
      slope = betta[0]
      intercept = 0
      
      line, = plt.plot([0.001,2],[slope*0.001+intercept,slope*2+intercept], '-', markersize = 2, color='black') 
      line.set_dashes([6, 3])   
      print " SLOP: ", slope, intercept

  plt.plot([-1], [-1], 'o', markersize = sizes[0], color=kolors[0], label=r'$5  \leqslant N  < 10$', markeredgecolor = kolors[0])
  plt.plot([-1], [-1], 'o', markersize = sizes[1], color=kolors[1], label=r'$10  \leqslant N  < 20$', markeredgecolor = kolors[1])
  plt.plot([-1], [-1], 'o', markersize = sizes[2], color=kolors[2], label=r'$20  \leqslant N  < 30$', markeredgecolor = kolors[2])
  plt.plot([-1], [-1], 'o', markersize = sizes[3], color=kolors[3], label=r'$30  \leqslant N  < 60$', markeredgecolor = kolors[3])
  plt.plot([-1], [-1], 'o', markersize = sizes[4], color=kolors[4], label=r'$60  \leqslant N$', markeredgecolor = kolors[4])
  lg = plt.legend(loc=4, numpoints = 1, prop={'size':15}, labelspacing=0.1, markerscale=1.5)
  lg.draw_frame(False)
  i = 0 
  for text in lg.get_texts():
    text.set_color(kolors[i])
    i+=1
 ##########################################
 
  if case==1: 
    xx = np.arange(0,2,0.001)
    plt.plot(xx, 368.25*xx, 'k-', linewidth=2, color='#CC4F1B')
  elif case==2: 
    xx = np.arange(0,1000,0.001)
    plt.plot(xx, (2.0E6)*xx**3, 'k-', linewidth=2, color='#CC4F1B')
 
 
  if case ==1:
    plt.minorticks_on()
    plt.tick_params(which='major', length=9, width=1.5)
    plt.tick_params(which='minor', length=6, color='#000033', width=1.0) 
    plt.xlim(0.01,2.0)
    plt.ylim(1.0,1000)
    plt.ylabel(r'$\sigma_p \/$'+'[km s'+r'$^{-1}$'+']', fontsize=18)
    plt.xlabel(r'$R_{2t}$'+' [Mpc]', fontsize=18)
    #lg = plt.legend( loc=2, numpoints = 1, prop={'size':17}, labelspacing=0.2, markerscale=1  )
    #lg.draw_frame(False)
  elif case ==2:
    plt.minorticks_on()
    plt.tick_params(which='major', length=9, width=1.5)
    plt.tick_params(which='minor', length=6, color='#000033', width=1.0) 
    plt.xlim(1,1000)
    plt.ylim(1.E8,2.E15)
    plt.xlabel(r'$\sigma_p \/$'+'[km s'+r'$^{-1}$'+']', fontsize=18)
    plt.ylabel('Virial Mass [M'+r'$_{\odot}$'+']', fontsize=18)
    #lg = plt.legend( loc=2, numpoints = 1, prop={'size':17}, labelspacing=0.2, markerscale=1  )
    #lg.draw_frame(False)    
    
  
  ax.set_yscale('log')
  ax.set_xscale('log') 
  
  if case ==1:
    labels = [0.01,0.1,1]
    ax.set_xticks(labels)
    labels_str = ['0.01','0.1','1']
    ax.set_xticklabels(labels_str)
    
    labels = [1,10,100]
    ax.set_yticks(labels)
    labels_str = ['1' ,'10','100']
    ax.set_yticklabels(labels_str)  
  elif case == 2:
    labels = [1,10,100]
    ax.set_xticks(labels)
    labels_str = ['1','10','100']
    ax.set_xticklabels(labels_str)

  
  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)  
  
  
  
  if case == 1:  
    ax.annotate(r'$\sigma_p/R_{2t}=368 \/\/ km/s/Mpc$', (0.021,1.85), color='black', fontsize=16)
    plt.plot([0.015,0.018], [2,2], 'k-', linewidth=2, color='#CC4F1B')
  elif case == 2:  
    ax.annotate(r'$M_V=2.0 \times 10^6 \sigma^3_p$', (2.15,4.1E12), color='black', fontsize=16)
    plt.plot([1.4,1.7], [5.0E12, 5.0E12], 'k-', linewidth=2, color='#CC4F1B') 
  
  plt.show()
  







