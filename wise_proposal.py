## Importing Important Python Libraries
import sys
import os
import random
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import  lines
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Polygon, Ellipse
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

import matplotlib.colors as colors
import matplotlib.cm as cmx

import matplotlib.patches as mpatches





from matplotlib import *


# **************************************
# Global Variables
# Physical Constants
# **************************************
H0 = 75.           # hubble constant
sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude
G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.08567758E19
#t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]
t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr
h = 0.75  # hubble constant
# **************************************
col_pallet = ['blue', 'green', 'orange', 'red']
vel_pallet = [-100000, 6000, 10000, 13000, 100000]

#################################################################
def color_table(Vls):
  
  for i in range(len(vel_pallet)-1):
    if vel_pallet[i] <= Vls and Vls < vel_pallet[i+1]:
      return col_pallet[i]
  return   colour  


def vel_limits(p, Vls):
  
  if vel_pallet[p] <= Vls and Vls < vel_pallet[p+1]:
    return True
  else: 
    return False


def vel_indx(Vls):
  
  for i in range(len(vel_pallet)-1):
    if vel_pallet[i] <= Vls and Vls < vel_pallet[i+1]:
      return i
  return 0
  
#################################################################
################################################################ 

def plot_border(l0):
  
  
  X = []
  Y = []

  
  l = 0.
  for b in np.arange(-90,90,0.01):
    x, y = xymap(l,b,l0)
    X.append(x)
    Y.append(y)
  l = 360
  for b in np.arange(90,-90,-0.01):
    x, y = xymap(l,b,l0)
    X.append(x)
    Y.append(y)
  
  X.append(X[0])
  Y.append(Y[0])

  plt.plot(X, Y, '-', markersize = 1, linewidth=2., color='black')   # '#0D1891'
  
################################################################ 
################################################################ 

def xymap(l,b,l0):
  
  t = theta(b)
  x = 180.-((2*90.)/180.)*(l0-l)*cos(t)
  y = 90.*sin(t)
  
  return x, y 


#################################################################
################################################################# 
# b is galactic latitude [deg]
# 2*theta + sin(2*theta) = Pi * sin(b)

def theta(b):
  
  if b == 90.: return pi/2.
  if b == -90.: return -pi/2.
  
  b = b*pi/180.
  theta0 = b
  theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b))/(2.+2.*cos(2.*theta0))
  
  while abs(theta1-theta0) > 0.01:
    theta0 = theta1
    theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b))/(2.+2.*cos(2.*theta0))
  
  return theta1
  
################################################################ 
def xymap_aitoff(x,y):
  
  while x > 360:
    x-=360
  while x < 0:
    x+=360
    
  x0 = (180.-x)*pi/180.
  y0 = y*pi/180.
  
  return x0, y0
################################################################ 
def xymap_aitoff_lst(X,Y):
  
  X_l = []
  X_r = []
  Y_l = []
  Y_r = []
  X_m = []
  Y_m = []
  for i in range(len(X)):
    x, y = xymap_aitoff(X[i],Y[i])
    if X[i]>360:
      X_r.append(x)
      Y_r.append(y)
    elif X[i]<0:
      X_l.append(x)
      Y_l.append(y)  
    else:
      X_m.append(x)
      Y_m.append(y)        
   
  OUT = []
  if len(X_l)>0:
    X_l = np.asarray(X_l)
    Y_l = np.asarray(Y_l)
    OUT.append([X_l, Y_l])
  if len(X_r)>0:
    X_r = np.asarray(X_r)
    Y_r = np.asarray(Y_r)
    OUT.append([X_r, Y_r])  
  if len(X_m)>0:
    X_m = np.asarray(X_m)
    Y_m = np.asarray(Y_m)
    OUT.append([X_m, Y_m])
    
  return OUT
################################################################ 


def esn_aitoff_patch(ax, x0, y0, d, color='blue'):

  
  
  vertices = []
  x, y = xymap_aitoff(x0,y0) 
  vertices.append([x,y])
  x, y = xymap_aitoff(x0,y0+d) 
  vertices.append([x,y])
  x, y = xymap_aitoff(x0+d,y0+d) 
  vertices.append([x,y])
  x, y = xymap_aitoff(x0+d,y0) 
  vertices.append([x,y])  
  
  
  
  ax.add_patch(Polygon(vertices, closed=True, fill=True, color=color))

##################################################################################################################################
def median_error(X):
  
  size = len(X)
  X    = np.sort(X) 
  mean = np.median(X)
  
  u_err = X[int(round(0.84*size))] - mean
  l_err = mean - X[int(round(0.16*size))]
  
  return mean, u_err, l_err

########################################################################## 

  
  
  
########################################################################## 
########################################################################## 
########################################################################## 
  
if __name__ == '__main__':
  
  l0 = 180

  deg = 0.8
  a = pyfits.open('EBV.0.8.deg.fits')    # source: dust_test.pro  --- for paper: (deg=0.8  EBV.0.8.deg.fits) 'EBV.10.deg.fits'
  d = a[1].data
  
  
  SGL  = d['SGL']
  SGB  = d['SGB']
  EBV = d['EBV']
  
  fig = plt.figure(figsize=(12, 7.5), dpi=100)
  ax = fig.add_subplot(111, projection="aitoff")
  #plt.title("Supergalactic Aitoff Projection", y=1.08)
  ax.grid(True)
  ax.set_xticklabels([])
  plt.subplots_adjust(top=0.98, bottom=0.02, right=0.95, left=0.05)
  

  
 
################################################################ 

  #ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
  #ax2.set_axis_off()
  #ax2.set_xticks([])
  #ax2.set_yticks([])
  #ax2.xaxis.set_ticks_position('none')
  #ax2.yaxis.set_ticks_position('none')
  #ax2.annotate(r"$\/\/ V_{ls} \/ < \/ 0 \/\/\/\/ km\/ s^{-1}  $", (0.1,0.2), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$ \/\/\/\/ 0 \/ \leq \/ V_{ls} \/ < \/ 250$', (0.1,0.2-0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$\/ 250 \/ \leq \/ V_{ls} \/ < \/ 500$', (0.1,0.2-2*0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$\/ 500 \/ \leq \/ V_{ls} \/ < \/ 750$', (0.1,0.2-3*0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$\/ 750 \/ \leq \/ V_{ls} \/ < \/ 1000$', (0.30,0.2), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$1000 \/ \leq \/ V_{ls} \/ < \/ 1250$', (0.30,0.2-0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$1250 \/ \leq \/ V_{ls} \/ < \/ 1500$', (0.30,0.2-2*0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$1500 \/ \leq \/ V_{ls} \/ < \/ 1750$', (0.30,0.2-3*0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$1750 \/ \leq \/ V_{ls} \/ < \/ 2000$', (0.5,0.2), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$2000 \/ \leq \/ V_{ls} \/ < \/ 2500$', (0.5,0.2-0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$2500 \/ \leq \/ V_{ls} \/ < \/ 3000$', (0.5,0.2-2*0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$3000 \/ \leq \/ V_{ls} \/ < \/ 3500$', (0.5,0.2-3*0.05), xycoords='figure fraction', size=12, color='black')
  #ax2.annotate(r'$3500 \/ \leq \/ V_{ls}$', (0.7,0.2), xycoords='figure fraction', size=12, color='black')
  
  #p = 0
  #for m in [0.1,0.3,0.5]:
      #for n in [0.2,0.2-0.05,0.2-2*0.05,0.2-3*0.05]:
          #ax2.add_patch(patches.Rectangle((m-0.05, n), 0.03, 0.03, color=col_pallet[p]))
          #p+=1
  #ax2.add_patch(patches.Rectangle((0.7-0.05, 0.2), 0.03, 0.03, color=col_pallet[p]))
################################################################ 
  
  median, u_err, l_err = median_error(EBV)
  
  X = []
  Y = []
  
  
  jet = cm = plt.get_cmap('Greys')
  cNorm  = colors.Normalize(vmin=0, vmax=5)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
  
  
  for i in range(len(SGB)):
    if EBV[i]>2.0: val = 3
    elif EBV[i]>1.0: val = 2.5
    elif EBV[i]>0.5: val= 2.0
    elif EBV[i]<=0.5 and EBV[i]>0.1: val= 1.0
    elif EBV[i]<0.1: val= 0
    else: val= 0

    colorVal = scalarMap.to_rgba(val)
    esn_aitoff_patch(ax, SGL[i], SGB[i], deg, color=colorVal)
    

  
################################################################ 
  #jet = cm = plt.get_cmap('jet')
  #cNorm  = colors.Normalize(vmin=0, vmax=3500)
  #scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
  #ax3 = fig.add_axes([0.7, 0.1, 0.2, 0.03])
  #cbar = colorbar.ColorbarBase(ax3, cmap=cm, norm=cNorm, orientation='horizontal',ticks=[0, 3500])
  ##cbar.set_label('km/s')
  #cbar.ax.set_xticklabels([r'$0$', r'$3500$'])  # horizontal colorbar
################################################################ 
################################################################ 
  jet = cm = plt.get_cmap('Greys')
  cNorm  = colors.Normalize(vmin=0, vmax=5)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
  ax3 = fig.add_axes([0.720, 0.11, 0.2, 0.03])
  cols = []
  bounds = [0,1,2,3]
  #cols.append(scalarMap.to_rgba(0))
  cols.append(scalarMap.to_rgba(1.0))
  cols.append(scalarMap.to_rgba(2.0))
  cols.append(scalarMap.to_rgba(2.5))
  #cols.append(scalarMap.to_rgba(3.0))
  
  
  
  cm = colors.ListedColormap(cols)
  
  cm.set_over(scalarMap.to_rgba(3.0))
  cm.set_under(scalarMap.to_rgba(0))
  
  cNorm  = colors.BoundaryNorm(bounds, cm.N)
  cbar = colorbar.ColorbarBase(ax3, cmap=cm, norm=cNorm, orientation='horizontal',ticks=bounds, boundaries=[-1]+bounds+[4],extend='both',extendfrac='auto')
  cbar.set_ticks(bounds)
  cbar.ax.tick_params(labelsize=10) 
  cbar.ax.set_xticklabels(['<0.1','0.5','1','2<'])  # horizontal colorbar
  
  cbar.set_label('E(B-V) '+r'$[mag]$')
################################################################   
  
  if True:
    



    mytable = np.genfromtxt('HIcand.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
    
    sgl   = mytable['sgl']
    sgb   = mytable['sgb'] 

    for i in range(len(sgl)):
      x, y = xymap_aitoff(sgl[i], sgb[i])
      l, = ax.plot([x], [y], 'o', markersize = 1.5, color='black', markeredgecolor = 'black')
      
      
      
    mytable = np.genfromtxt('cf3_nohead.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
    
    Vls   = mytable['Vls']
    sgl   = mytable['sgl']
    sgb   = mytable['sgb']
    

    for i in range(len(Vls)):
      x, y = xymap_aitoff(sgl[i], sgb[i])
      color = color_table(Vls[i])
      l, = ax.plot([x], [y], 'o', markersize = 1.5, color=color, markeredgecolor = color)      
################################################################ 
  my_color = 'black'
  
  ax.annotate(r'$0^o$', (pi-0.1,pi/3.), size=12, color=my_color)
  ax.annotate(r'$90^o$', (pi/2.-0.1,pi/3.), size=12, color=my_color, backgroundcolor='white')
  ax.annotate(r'$180^o$', (-0.2,pi/3.), size=12, color=my_color)
  ax.annotate(r'$270^o$', (-pi/2.-0.1,pi/3.), size=12, color=my_color, backgroundcolor='white')
  
  #plt.show()
  plt.savefig('wise_proposal.png', dpi=600)  # the best view and size happens in eps format









