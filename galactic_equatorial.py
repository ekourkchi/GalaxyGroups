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

from astropy import coordinates as coord
from astropy import units as unit

################################################################# 
# l & b are in degree
def galeq_converter(l, b):
  
  l1 = (l-123) * pi/180.
  b0 = 27.4 * pi/180.
  b = b * pi/180.
  
  alpha = atan2(sin(l1),(cos(l1)*sin(b0)-tan(b)*cos(b0)))
  #print alpha
  alpha = (alpha)+ 12.25*pi/180.
  
  
  delta = asin(sin(b)*sin(b0)+cos(b)*cos(b0)*cos(l1))
  
  
  alpha = alpha + 0.640265*pi/180 + (0.278369*pi/180) * sin(alpha) * tan(delta)
  delta = delta + 0.278369*pi/180
  
  alpha = alpha * 180. / pi
  if alpha<0: alpha+=360.
  delta = delta * 180. / pi

  
  return alpha, delta
  
  
  
  
################################################################ 

if __name__ == '__main__':
  
  #point = coord.Galactic(96.8794600 , -59.5418000, unit=(u.degree, u.degree))
  #point.fk5.ra.degree
  #point.fk5.dec.degree

  table = np.genfromtxt( 'AllSky.south.csv' , delimiter=',', filling_values=0, names=True, dtype='float64')
  
  #print table.dtype
  id   = table['pgc']
  gl  = table['gl']
  gb  = table['gb']
  N_galaxies = len(id)
  
  RA0 = []
  DEC0 = []
  for i in range(0, N_galaxies):
    point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
    RA0.append(point.fk5.ra.degree)
    DEC0.append(point.fk5.dec.degree)

  plt.plot(RA0, DEC0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  plt.ylim(-45,-10)
  plt.xlim(70,35)
  
  plt.show()
  
  