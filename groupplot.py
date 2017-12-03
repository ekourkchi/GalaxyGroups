#!/home/ehsan/Ureka/Ureka/variants/common/bin/python




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
def plot_header():
  fig = plt.figure(figsize=(7, 7), dpi=100)
  ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
  ax.xaxis.set_major_locator(MultipleLocator(10))
  ax.yaxis.set_major_locator(MultipleLocator(10))
  ax.xaxis.set_minor_locator(MultipleLocator(1))
  ax.yaxis.set_minor_locator(MultipleLocator(1)) 
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)


def ClusterPlot(id, ra, dec, Vls, R_2t2, nest, flag, isVirgo=False, R_min=6.0, R_max=100000000.):
  
  N = len(id)
  NoGroups = len(id[np.where(flag==2)])
  print "Number of groups: ", NoGroups
  
   
  
  
  if not isVirgo:
    plt.ylim(-45,-10)
    plt.xlim(70,35)
    plt.xlabel("RA (deg)", fontsize=20)
    plt.ylabel("DEC (deg)", fontsize=20)
  if isVirgo:
    xmin = 5*floor((sglV+35)/5.)+3
    xmax = 5*ceil((sglV-35)/5.)-3
    #ymax = 5*floor((sgbV+35)/5.)+3
    ymin = 5*ceil((sgbV-35)/5.)-3
    ymax = ymin - (xmax-xmin)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
  
  
    plt.xlabel("SGL (deg)", fontsize=20)
    plt.ylabel("SGB (deg)", fontsize=20)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
  
    # plotting big slod circle
    theta = np.arange(0,1,0.001)
    Circlx = 30.*np.cos(theta*2*pi) + sglV
    Circly = 30.*np.sin(theta*2*pi) + sgbV
    plt.plot(Circlx, Circly, 'k-')
  
    # plotting small balck dotted circle
    theta = np.arange(0,1,0.010)
    Circlx = 6.8*np.cos(theta*2*pi) + sglV
    Circly = 6.8*np.sin(theta*2*pi) + sgbV
    line00, = plt.plot(Circlx, Circly, 'k.', markersize = 2)   
    line00.set_dashes([2, 2]) 
  
    # plotting small red dotted circle
    Circlx = R_min*np.cos(theta*2*pi) + sglV
    Circly = R_min*np.sin(theta*2*pi) + sgbV
    line0, = plt.plot(Circlx, Circly, 'r.', markersize = 2)   
    line0.set_dashes([2, 3]) 
  
    # plotting small red dotted circle
    Circlx = R_max*np.cos(theta*2*pi) + sglV
    Circly = R_max*np.sin(theta*2*pi) + sgbV
    line0, = plt.plot(Circlx, Circly, 'g.', markersize = 2)   
    line0.set_dashes([2, 3])  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  gr = flag[0]
  i = 0 
  while gr == 2:
    
    random.seed(nest[i])
    Red, Green, Blue = random.random(), random.random(), random.random()
    Dist = Vls[i]/H0
    if Dist<1: Dist = 1.
    r = 180.*atan(R_2t2[i]/Dist)/pi
    
    d_theta = 0.01
    u = np.arange(0,2*pi,d_theta)
    X = u*0.
    Y = u*0.
    
    RA0 = ra[i]
    DEC0 = dec[i]
    
    for q in range(0,len(u)):
       x = r*np.cos(u[q]) + RA0
       y = r*np.sin(u[q]) + DEC0
       X[q] = x
       Y[q] = y
       
    if r <= 10000:
	if nest[i] == 13418:   # NGC_1399
	    Blue = 0; Red = 0.2 ; Green = 0.8
    
        line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Green, Blue)) 
	line.set_dashes([8, 3]) 
	#plt.text(RA0, DEC0+r, int(nest[i]), fontsize=8, color=(Red, Green, Blue))
    
    
    if nest[i] == 13418:   # NGC_1399
      star_ra = ra[i]
      star_dec = dec[i]
    i+=1  
    X = []
    Y = []
    while i<N and flag[i]==1: 
      RA0 = ra[i]
      DEC0 = dec[i]
      X.append(RA0)
      Y.append(DEC0)
      plt.plot(ra[i], dec[i], 'o', markersize = 3, color=(Red, Green, Blue), markeredgecolor = 'none')
      i+=1
    #plt.plot(X, Y, 'o', markersize = 3, color=(Red, Green, Blue), markeredgecolor = 'none')
    
        
    
    
    if i<N and flag[i]==2: 
      gr = 2
    else:
       while i<N and flag[i]==0:
         #plt.plot(ra[i], dec[i], '.', markersize = 1, color='#696969')
	 i+=1
       break
  
  if not isVirgo:
     plt.plot(star_ra, star_dec, '*', markersize = 10, color= 'black')
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
def xymap(l,b,l0):
  
  t = theta(b)
  x = 180.-((2*90.)/180.)*(l0-l)*cos(t)
  y = 90.*sin(t)
  
  return x, y 
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

  plt.plot(X, Y, '-', markersize = 1, linewidth=2., color='#0D1891')
  
  
  #X = []
  #Y = []
  #b = 75.
  #for l in np.arange(0,360,0.01):
    #x, y = xymap(l,b,l0)
    #X.append(x)
    #Y.append(y)
  #plt.plot(X, Y, '-', markersize = 1, linewidth=2., color='#0D1891')
  
  
  
################################################################ 

def plot_galaxies(inFile, l0):
  
  table = np.genfromtxt( inFile , delimiter=',', filling_values=0, names=True, dtype='float64')
  
  #print table.dtype
  id   = table['pgc']
  gl  = table['sgl']
  gb  = table['sgb']
  Vls  = table['Vls']
  N_galaxies = len(id)
  
  X0 = []
  Y0 = []
  for i in range(0, N_galaxies):
     if Vls[i]<=3500:
        x, y = xymap(gl[i],gb[i],l0)
        X0.append(x)
        Y0.append(y)
  plt.plot(X0, Y0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  
################################################################ 
def groupPlot(north_catal, south_catal, galaxies_north, galaxies_south):
  
  l0 = 180
  
  fig = plt.figure(figsize=(12, 12*0.5), dpi=100)
  ax = fig.add_axes([0.13, 0.13, 0.80,  0.80]) 
  plt.ylim(-100,100)
  plt.xlim(380,-20)
  plt.xlabel("SGL (deg)", fontsize=20)
  plt.ylabel("SGB (deg)", fontsize=20)
  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)
  
  plot_border(l0)
  plot_galaxies(galaxies_north, l0)
  plot_galaxies(galaxies_south, l0)
  
  # Virgo Border (6.8 deg)
  u = np.arange(0,1,0.010)
  X = u*0.
  Y = u*0.
  for i in range(0,len(u)):
     x = 6.8*np.cos(u[i]*2*pi) + sglV
     y = 6.8*np.sin(u[i]*2*pi) + sgbV
     X[i], Y[i] = xymap(x,y, l0)
  
  line00, = plt.plot(X,Y, 'r.', markersize = 2)   
  line00.set_dashes([2, 2]) 
  
  groupPlotak(north_catal, l0)
  groupPlotak(south_catal, l0)   
  plot_border(l0)
##########

def groupPlotak(group_file, l0):
  
  mytable = np.genfromtxt(group_file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  id    = mytable['pgc']
  flag  = mytable['flag']
  sgl   = mytable['sgl']
  sgb   = mytable['sgb']
  Vls   = mytable['Vls']
  nest  = mytable['nest']
  
  V = Vls[0]
  for i in range(len(id)):  # for all gaalxies in groups
    if flag[i] == 2: V = Vls[i]
    if flag[i] == 1 and V<=3500:
      random.seed(nest[i])
      Red, Blue, Green = random.random(), random.random(), random.random()
      x, y = xymap(sgl[i], sgb[i], l0)
      plt.plot(x, y, 'o', markersize = 3, color=(Red, Blue, Green), markeredgecolor = 'none')
        

################################################################# 

################################################################ 

if __name__ == '__main__':
  
  
  galaxies_north = 'AllSky.north.v11.csv'
  galaxies_south = 'AllSky.south.v11.csv'
  
  
  #cluster = 'all'
  #cluster = 'virgo'
  cluster = 'fornax'

  version = 'v21'
  iteration = '9'
  
  north_catal = 'north.iter.'+iteration+'.'+version+'.group'
  south_catal = 'south.iter.'+iteration+'.'+version+'.group'

  virgo_catal  =  'virgo.iter.'+iteration+'.v21.group'
  fornax_catal =  'fornax.iter.'+iteration+'.v21.group'
  
  
  
  
  
  if cluster=='all':
    groupPlot(north_catal, south_catal, galaxies_north, galaxies_south)
  
  if cluster=='fornax' or cluster=='virgo':
    plot_header()   
  
  if cluster == 'virgo':
    file = virgo_catal
    inFile  = galaxies_north
    table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
    id   = table['pgc']
    sgl  = table['sgl']
    sgb  = table['sgb']
    Vls  = table['Vls']
    flag = np.where(Vls<=3500)
    sgl = sgl[flag]
    sgb = sgb[flag]
    Vls = Vls[flag]
    plt.plot(sgl, sgb, '.', markersize = 1, color='#696969')
  
  
  if cluster=='fornax':
    file = fornax_catal
    inFile  = galaxies_south
    table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
    id   = table['pgc']
    ra  = table['ra']
    dec  = table['dec']
    Vls  = table['Vls']
    flag = np.where(Vls<=3500)
    ra = ra[flag]
    dec = dec[flag]
    Vls = Vls[flag]
    plt.plot(ra, dec, '.', markersize = 1, color='#696969')  
  
  
  
  
  if cluster=='fornax' or cluster=='virgo':
    mytable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
    id    = mytable['pgc']
    flag  = mytable['flag']
    sgl   = mytable['sgl']
    sgb   = mytable['sgb']
    ra   = mytable['ra']
    dec  = mytable['dec']
    Vls   = mytable['Vls']
    R_2t2 = mytable['R2t_lum']
    nest  = mytable['nest']
  
  if cluster=='fornax':
    ClusterPlot(id, ra, dec, Vls, R_2t2, nest, flag)
  
  if cluster == 'virgo':
    ClusterPlot(id, sgl, sgb, Vls, R_2t2, nest, flag, isVirgo=True, R_max=50.)
  

  
  
  
  
  
  
  
  
  
  #plt.savefig('test.eps', dpi=600) 
  
  plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  





