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
MM = []
LL = []
#################################################################
def schechter(x, rho, Ms):
    n = 0
    #x = c * x
    return (rho/sqrt(pi)/x/x)*(1.+n/3.)*((x/Ms)**(0.5+n/6.))*np.exp(-1.*(x/Ms)**(1.+n/3.))
  
  
def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    
    rho = B[0]
    Ms = B[1]
    n = 0 
    
    return (rho/sqrt(pi)/x/x)*(1.+n/3.)*((x/Ms)**(0.5+n/6.))*np.exp(-1.*(x/Ms)**(1.+n/3.))
################################################################# 

def numberlogbin(x, x_min=1.E11, x_max=1.E16, C=0.2):
  
  
  x_bin = []
  y_bin = []
  x_err = []
  y_err = []
  
  X1 = x_min
  X2 = X1 * 10**C
  bin = []
  count = []
  while(X2<=x_max):

        
        N = len(x)
        filter = np.zeros((N,), dtype=np.int)
        f1 = np.zeros((N,), dtype=np.int)
        f1[np.where(x<X1)] = 1
        f2 = np.zeros((N,), dtype=np.int)
        f2[np.where(x>=X2)] = 1
        filter = f1 + f2
        indices = np.where(filter==0)
        x_filt = x[indices]
        
        bin.append(sqrt(X1*X2))
        count.append(len(x_filt))

	X1 = X2
	X2 = X1 * 10**C

  count = np.asarray(count)
  bin = np.asarray(bin)
  return bin, count


#################################################################
# **************************************

def Mass(L_k):

  L = L_k / 1.E10
  
  if L < 4.:
    MtoL = 80.24*(L**-0.30)
  elif L > 1000.:
    MtoL = 121.19
  else:
    MtoL = 43*(L**0.15)
  
  h = 0.75  # hubble constant
  Mass_out = h * L_k * MtoL
  
  return Mass_out

# **************************************
def Lumin(M):

  M12 = M / 1.E12
  Lumin = 1.E10*((M/43.E10)**(1./1.15))*exp(-0.6/M12)
  return Lumin

def Mass2(L_k):
  
  global MM
  global LL
  if len(MM) == 0 or len(LL) == 0:
     mass = 1.E6
     while(mass<=1.E20):
	MM.append(mass)
        LL.append(Lumin(mass))
        mass = mass * 10**0.05
     MM = np.asarray(MM)
     LL = np.asarray(LL)
  
  i = 0
  while(LL[i]<L_k): i+=1
    
  Mass = MM[i-1] + (MM[i]-MM[i-1])*(L_k-LL[i-1])/(LL[i]-LL[i-1])
    
  return Mass
################################################################# 

if __name__ == '__main__':
  
  fig = plt.figure(figsize=(7, 7), dpi=100)
  ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0) 
  #################################################################
  file = 'north.ML.64.v09.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  flag = Gtable['flag']
  Mv_lum = Gtable['Mv_lum']
  Vls  = Gtable['Vls']
  logK = Gtable['logK']
  #################################################################
  file = 'south.ML.64.v09.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  flag0 = Gtable['flag']
  Mv_lum0 = Gtable['Mv_lum']
  Vls0  = Gtable['Vls']
  logK0 = Gtable['logK']
  
  flag = np.concatenate((flag,flag0))
  Mv_lum = np.concatenate((Mv_lum,Mv_lum0))
  Vls = np.concatenate((Vls,Vls0))
  logK = np.concatenate((logK,logK0))
  #################################################################
  for i in range(0,len(logK)):
    Mv_lum[i] = Mass(10**logK[i])
    
  
  
  
  
  N = len(flag)
  filter = np.zeros((N,), dtype=np.int)
  
  f3 = np.zeros((N,), dtype=np.int)
  f3[np.where(flag==1)] = 1 
  
  

  
  f2 = np.zeros((N,), dtype=np.int)
  f2[np.where(Vls>800)] = 1 
  f1 = np.zeros((N,), dtype=np.int)
  f1[np.where(Vls<100)] = 1
  filter = f1 + f2 + f3
  indices = np.where(filter==0)
  Mv0 = Mv_lum[indices]
  
  
  
  f2 = np.zeros((N,), dtype=np.int)
  f2[np.where(Vls>2000)] = 1 
  f1 = np.zeros((N,), dtype=np.int)
  f1[np.where(Vls<100)] = 1
  filter = f1 + f2 + f3
  indices = np.where(filter==0)
  Mv1 = Mv_lum[indices]
   
  
  f2 = np.zeros((N,), dtype=np.int)
  f2[np.where(Vls>3500)] = 1 
  f1 = np.zeros((N,), dtype=np.int)
  f1[np.where(Vls<100)] = 1
  filter = f1 + f2 + f3
  indices = np.where(filter==0)
  Mv2 = Mv_lum[indices]
    
  
  
  
 
  
  
  

  
  bins, number =  numberlogbin(Mv2, x_min=3.E8, x_max=max(Mv2), C=0.2)
  indices = np.where(bins>3.E10)
  bins = bins[indices]
  number = number[indices]
  n = number[3:9]
  n = np.median(n)
  bins0 = bins[3::]
  number0 = number[3::]
  factor = 3000/number[9]
  #plt.plot(bins[0:3], factor*number[0:3], 'o', markersize = 4, color='grey', markeredgecolor = 'none') 
  plt.errorbar(bins[0:3], factor*number[0:3], yerr=np.sqrt(np.asarray(factor*number[0:3])), capsize=2, ls='none', color='grey', elinewidth=1, fmt='o', markersize = 4, markeredgecolor = 'none')
  
 
  bins1, number1 =  numberlogbin(Mv1, x_min=3.E8, x_max=max(Mv2))
  indices = np.where(bins1>9.E9)
  bins1 = bins1[indices]
  number1 = number1[indices]
  #plt.plot(bins1, factor*number1, 'o', markersize = 4, color='green', label='800 <v<2000 km/s') 
  m = np.median(number1[6:12])  
  number1 = number1*n/m
  
  #plt.plot(bins1[0:4], factor*number1[0:4], '^', markersize = 4, color='grey', markeredgecolor = 'none')
  plt.errorbar(bins1[0:4], factor*number1[0:4], yerr=np.sqrt(np.asarray(factor*number1[0:4])), capsize=2, ls='none', color='grey', elinewidth=1, fmt='^', markersize = 4, markeredgecolor = 'none')
  #plt.plot(bins1[6::], factor*number1[6::], '^', markersize = 4, color='grey', markeredgecolor = 'none')
  plt.errorbar(bins1[6::], factor*number1[6::], yerr=np.sqrt(np.asarray(factor*number1[6::])), capsize=2, ls='none', color='grey', elinewidth=1, fmt='^', markersize = 4, markeredgecolor = 'none')


  
  number = np.concatenate((number1[4:6],number[3::]))
  bins = np.concatenate((bins1[4:6],bins[3::]))
  #plt.plot(bins, factor*number, '-', color='orange')
  
  
  bins2, number2 =  numberlogbin(Mv0, x_min=3.E8, x_max=max(Mv2))
  #plt.plot(bins2, factor*number2, 'o', markersize = 4, color='blue', label='100 <v<800  km/s')
  p = np.mean(number2[11:12])
  q = np.mean(number[0:1])
  number2 = number2*q/p
  
  #plt.plot(bins2[12::], factor*number2[12::], 's', markersize = 3, color='grey', markeredgecolor = 'none')
  plt.errorbar(bins2[12::], factor*number2[12::], yerr=np.sqrt(np.asarray(factor*number2[12::])), capsize=2, ls='none', color='grey', elinewidth=1, fmt='s', markersize = 4, markeredgecolor = 'none')
  
  
  number = np.concatenate((number2[0:11],number))
  bins = np.concatenate((bins2[0:11],bins))
  #plt.plot(bins, number*factor, '-', color='orange')

  ##############################################################
  
  indices = np.where(bins>4.E10)
  bins = bins[indices]
  number = number[indices]
  
  alf = 10**0.2
  delta = ((alf-1)/sqrt(alf))*bins

  
  popt, pcov = curve_fit(schechter, bins, number/delta, [1.1E16, 3.4E13], sigma = 1./bins/bins, maxfev = 1000000)
  print popt
  
  
  bi =[]
  n0 = []
  b = 1.E8
  while b<=1.E16:
    bi.append(b)
    delt = ((alf-1)/sqrt(alf))*b
    #n0.append(delt*schechter(b, 10.E15, 2.E13, 0.9))
    #n0.append(delt*schechter(b, popt[0], popt[1], popt[2]))
    n0.append(delt*schechter(b, popt[0], popt[1]))
    b*= alf
  bi = np.asarray(bi)
  n0 = np.asarray(n0)
  plt.plot(bi, factor*n0, '--', color='black', label='Press-Schechter')
  
  
  
  ## url: http://docs.scipy.org/doc/scipy/reference/odr.html
  ## Orthogonal distance regression
  #linear = odr.Model(f)
  #mydata = odr.Data(bins, number/delta, we=1./bins/bins/bins)
  #myodr = odr.ODR(mydata, linear, beta0=[1.1E15, 3.4E12])
  #myoutput = myodr.run()
  #slope = myoutput.beta
  #print slope
  #n0 = delt*f(slope, bi)
  #plt.plot(bi, factor*n0, ':', color='orange', label='Orthogonal')  
  
  ##############################################################

  #plt.plot(bins0, factor*number0, 'o', markersize = 5, color='red', label='100<v<3500 km/s', markeredgecolor = 'red') 

  plt.errorbar(bins0, factor*number0, yerr=np.sqrt(np.asarray(factor*number0)), capsize=2, ls='none', color='red', elinewidth=1, fmt='o', label='100<v<3500 km/s', markersize = 5, markeredgecolor = 'red') 
  
  plt.errorbar(bins1[4:6], factor*number1[4:6], yerr=np.sqrt(np.asarray(factor*number1[4:6])), fmt='^', markersize = 5, color='green', label='100 <v<2000 km/s', markeredgecolor = 'green') 
  
  plt.errorbar(bins2[0:11], factor*number2[0:11], yerr=np.sqrt(np.asarray(factor*number2[0:11])), fmt='s', markersize = 4, color='blue', label='100 <v<800  km/s', markeredgecolor = 'blue')
  
  plt.plot([bins1[4],bins1[4]], [factor*number1[4],factor*number1[4]/10**0.8], ':', color='green', linewidth=2)
  plt.plot([bins0[0],bins0[0]], [factor*number0[0],factor*number0[0]/10**0.8], ':', color='red', linewidth=2)
  
 
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('Group Mass [M'+r'$_{\odot}$'+']', fontsize=18)  
  plt.ylabel("Number", fontsize=18)  
  plt.xlim(1.e8,1.e16) 
  plt.ylim(0.1,50000) 
  
  labels = [0.1,1,10,100,1000,10000]
  ax.set_yticks(labels)
  labels_str = ['0.1','1','10','100','1000','10'+r'$^4$']
  ax.set_yticklabels(labels_str)
  
  
  
  
  
  lg = plt.legend( loc=3, numpoints = 1, prop={'size':14}, labelspacing=0.2, markerscale=2  )
  lg.draw_frame(False)
  plt.show()
  
  
  
  
  
  
  
  
  
  
  
  

