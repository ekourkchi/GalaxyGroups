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
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
################################################################ 
l0 = 180.
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
  #if l>180: l=360.-l
  x = 180.-((2*90.)/180.)*(l0-l)*cos(t)
  y = 90.*sin(t)
  
  return x, y 
################################################################ 

def plot_border():
  
  
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

  plt.plot(X, Y, '-', markersize = 1)
  
################################################################ 

def plot_galaxies(inFile):
  
  table = np.genfromtxt( inFile , delimiter=',', filling_values=0, names=True, dtype='float64')
  
  #print table.dtype
  id   = table['pgc']
  gl  = table['sgl']
  gb  = table['sgb']
  N_galaxies = len(id)
  
  X0 = []
  Y0 = []
  for i in range(0, N_galaxies):
    x, y = xymap(gl[i],gb[i],l0)
    X0.append(x)
    Y0.append(y)
  plt.plot(X0, Y0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  
  
################################################################ 
  
def plot_3d():
  
  
  u = np.linspace(0, np.pi, 30)
  v = np.linspace(0, 2 * np.pi, 30)

  x = np.outer(np.sin(u), np.sin(v))
  y = np.outer(np.sin(u), np.cos(v))
  z = np.outer(np.cos(u), np.ones_like(v))


  fig = plt.figure()
  ax = plt.axes(projection='3d')

  #ax.plot_wireframe(x, y, z)
  surf = ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap=cm.coolwarm,  #facecolors=colors,
        linewidth=0, antialiased=False)

  
  inFile = 'AllSky.north.csv'
  table = np.genfromtxt( inFile , delimiter=',', filling_values=0, names=True, dtype='float64')
  
  #print table.dtype
  id   = table['pgc']
  gl  = table['gl']
  gb  = table['gb']
  N_galaxies = len(id)
  
  X0 = []
  Y0 = []
  Z0 = []
  for i in range(0, N_galaxies):
    #if gl[i]>180.:  gl[i] = 360.- gl[i]
    X0.append(cos(gb[i]*pi/180.)*cos(gl[i]*pi/180.))
    Y0.append(cos(gb[i]*pi/180.)*sin(gl[i]*pi/180.))
    Z0.append(sin(gb[i]*pi/180.))
    
  #plt.plot(X0, Y0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  ax.scatter(X0,Y0,Z0, color='#696969')
   

################################################################ 

if __name__ == '__main__':
  
  fig = plt.figure(figsize=(9, 9*0.5), dpi=100)
  ##ax = fig.add_axes([0.13, 0.1, 0.85,  0.85]) 
  plt.ylim(-90,90)
  plt.xlim(380,-20)
  
  
  plot_border()
  plot_galaxies('AllSky.north.csv')
  plot_galaxies('AllSky.south.csv')
  
  

  
  
  plt.show()
  
  
  
  
  
  
  
  
  
  