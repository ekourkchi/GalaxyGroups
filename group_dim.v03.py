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
  Bi_Weight = True
  #if len(sys.argv) == 2 : Bi_Weight = False
  #elif int(sys.argv[2]) == 0 : Bi_Weight = False
  
  #################################################################
  file = 'all.iter.2.v44.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  
  flag = Gtable['flag']
  No = Gtable['No_Galaxies']
  R_size = Gtable['R2t_lum']
  Vls = Gtable['Vls']
  R2t_dyn = Gtable['R2t_dyn']
  R2t_lum = Gtable['R2t_lum']
  sigmaP_dyn = Gtable['sigmaP_dyn']
  sigmaP_lum = Gtable['sigmaP_lum']
  mDist = Gtable['mDist']
  gl = Gtable['gl']
  gb = Gtable['gb']
  Mv_dyn = Gtable['Mv_dyn'] 
  Mv_lum = Gtable['Mv_lum'] 
  #################################################################

  #################################################################
  X1 = R2t_dyn /(pi*0.5/1.05/sqrt(1.5))
  Y1 = R2t_lum
   
  
  X2 = sigmaP_dyn
  Y2 = sigmaP_lum
  
  X3 = Mv_dyn
  Y3 = Mv_lum
#################################################################
  
  if Bi_Weight:
    X1 = R2t_dyn
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
	if isnan(Rg_bi): Rg_bi = R2t_dyn[j] /(pi*0.5/1.05/sqrt(1.5))
	X1[j] = Rg_bi
	X2[j] = sigma
	i-=1
      i+=1
    X3 = (1.E9 * (2.5*pi/G/2.) * (Mpc_km*X1) * X2**2)/M_sun  # solar mass
#################################################################
  
  
  if case==1: X = X1 ; Y = Y1   # R_2t   dyn vs. lum
  if case==2: X = X2 ; Y = Y2   # sigmaP dyn vs. lum
  if case==3: X = Y3 ; Y = X3   # Mv     dyn vs. lum
  
  
  
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
  X = X[indices]
  Y = Y[indices]
  No = No[indices]
  R_size = R_size[indices]
  Vls = Vls[indices]
  #################################################################

  
  
  
  minR = min(R_size)
  maxR = max(R_size)
  index = np.where(No>=5)
  x_size = X[index]
  y_size = Y[index]
  maxVel = max(max(x_size), max(y_size))


  
  fig = plt.figure(figsize=(7, 7), dpi=100)
  ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)  
  
  if case==2:
     plt.xlabel("Dynamical Velocity Dispersion ("+r'$\sigma_V$'+") [km/s]", fontsize=18)
     plt.ylabel("Luminosity Dispersion ("+r'$\sigma_L$'+") [km/s]", fontsize=18)
     if Bi_Weight: plt.xlabel("Bi-Weight Velocity Dispersion ("+r'$\sigma_V$'+") [km/s]", fontsize=18) 
     maxVel = 100.*ceil(maxVel/100.)+50
     plt.xlim(0,maxVel)
     plt.ylim(0,maxVel)      
     line, = plt.plot([0,maxVel],[0,maxVel], '-', markersize = 2, color='black')
     ax.annotate('slope=0.83'+r'$\pm$'+'0.04', xy=(533, 428), xytext=(582, 326), arrowprops=dict(facecolor='black', shrink=0.01))
     
  elif case==1:
     plt.xlabel("Gravitational Radius [Mpc]", fontsize=18)  
     plt.ylabel('2'+r'$^{nd}$'+' Turnaround Radius (R'+r'$_{2t}$'+') [Mpc]', fontsize=18)
     if Bi_Weight: plt.xlabel("Bi-Weight Gravitational Radius (R"+r'$_g$'+") [Mpc]", fontsize=18) 
     maxVel = 1.5*ceil(maxVel/1.5)
     plt.xlim(0,maxVel) 
     plt.ylim(0,maxVel)
     line, = plt.plot([0,maxVel],[0,maxVel], '-', markersize = 2, color='black')
     plt.xlim(0,1.8) 
     plt.ylim(0,1.8) 
     
     ax.annotate('slope=1.44'+r'$\pm$'+'0.06', xy=(0.75, 1.1), xytext=(0.21, 1.35), arrowprops=dict(facecolor='black', shrink=0.01))
     
  elif case==3:
     plt.xscale('log')
     plt.yscale('log')
     plt.xlabel('Group Luminosity Mass [M'+r'$_{\odot}$'+']', fontsize=18)  
     plt.ylabel('Group Virial Mass [M'+r'$_{\odot}$'+']', fontsize=18)  
     
     
     
     plt.xlim(1.e11,2.e15) 
     plt.ylim(1.e11,2.e15) 
     line, = plt.plot([1.e11,1.e16],[1.e11,1.e16], '-', markersize = 2, color='black') 
  
  line.set_dashes([2, 3]) 
  kolors = ['gray', 'blue', 'green', 'orange', 'brown']
  sizes = [3,5,6,7,8]
  for i in range(0,len(flag)):
     #random.seed(No[i])
     #Red, Blue, Green = random.random(), random.random(), random.random()
     if No[i]>=5 and No[i]<10: j=0
     if No[i]>=10 and No[i]<20: j=1
     if No[i]>=20 and No[i]<30: j=2
     if No[i]>=30 and No[i]<60: j=3
     if No[i]>=60 : j=4
     if No[i]>=5:
       #size = floor(5.*(R_size[i]-minR)/(maxR-minR)+3)
       plt.plot(X[i], Y[i], 'o', markersize = sizes[j], color=kolors[j], markeredgecolor = kolors[j])
  
  
  


     
  
  if case < 3:
      index = np.where(No>=30)
      x = X[index]
      y = Y[index] 
      weight = No[index] 
      #slopee, intercept, r_value, p_value, std_err = stats.linregress(x,y)
      
      # url: http://docs.scipy.org/doc/scipy/reference/odr.html
      # Orthogonal distance regression
      linear = odr.Model(f)
      mydata = odr.Data(x, y, we=weight)
      myodr = odr.ODR(mydata, linear, beta0=[1.])
      myoutput = myodr.run()
      slopee = myoutput.beta
      intercept = 0
      
      line, = plt.plot([0,maxVel],[intercept,slopee*maxVel+intercept], '-', markersize = 2, color='black') 
      line.set_dashes([6, 3])   
      print " SLOP: ", slopee, "STD: ", myoutput.sd_beta

  
  
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
    
  
  if case < 3 : 
    index = np.where(No>=10)
    x = X[index]
    y = Y[index]
    print "difference in radii: ", np.mean((y-x)/x) , '  Slope:', slopee
    myoutput.pprint()
  
  plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
