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


from kapteyn import wcs



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
col_pallet = ['darkviolet', 'blue', 'deepskyblue', 'forestgreen', 'y', 'gold', 'darkorange', 'red', 'magenta', 'maroon', 'sienna', 'slategrey', 'black']
vel_pallet = [-100000, 0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 100000]

#################################################################

def add_plane(ax, color='black', plane=None, projection=None):
  
  if plane==None or projection==None:
    return
  
  alpha = np.arange(0.,360,2)
  delta = alpha*0.
  
  tran = wcs.Transformation(plane + " j2000 j2000", projection)
  alpha, delta = tran((alpha,delta))
  
  for i in range(len(alpha)):
    if alpha[i] >180:
      alpha[i] -= 360.
  
  ind = np.argsort(alpha)
  alpha = alpha[ind]
  delta = delta[ind]
  ax.plot(alpha*pi/180, delta*pi/180, '-', color=color)
  
  
########################################################################## 
########################################################################## 
########################################################################## 
  
if __name__ == '__main__':
  
  if len(sys.argv) < 2:
    print "\nEnter the cluster name as the input ..." 
    print "\nexample: \n > python " + str(sys.argv[0])+ " manual"
    print "\nPossible options: Equatorial, Galactic, Supergalactic, Ecliptic" 
    print "Use north/south for whole sky.\n"
    sys.exit(1)
  
  projection = str(sys.argv[1])
  
  
  
  
  fig = plt.figure(figsize=(12, 12*0.5), dpi=100)
  ax = fig.add_subplot(111, projection="aitoff")
  ax.grid(True)
  ax.set_xticklabels([])
  plt.title(projection+" - Aitoff Projection", y=1.08)
  plt.subplots_adjust(top=0.85, bottom=0.05, right=0.95, left=0.05)
  
  #ax.annotate(r'$0^o$', (pi-0.1,pi/3.), size=11, color='black')
  #ax.annotate(r'$90^o$', (pi/2.-0.2,pi/3.), size=11, color='black')
  #ax.annotate(r'$180^o$', (-0.2,pi/3.), size=11, color='black')
  #ax.annotate(r'$270^o$', (-pi/2.-0.2,pi/3.), size=11, color='black')
  
  
################################################################ 
  ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
  ax2.set_axis_off()
  ax2.set_xticks([])
  ax2.set_yticks([])


################################################################ 

  if projection == 'Equatorial':
    add_plane(ax, color='b', plane='galactic', projection='equatorial')
    add_plane(ax, color='r', plane='ecliptic', projection='equatorial')
    add_plane(ax, color='g', plane='supergalactic', projection='equatorial')
  
  if projection == 'Galactic':
    add_plane(ax, color='black', plane='equatorial', projection='galactic')
    add_plane(ax, color='r', plane='ecliptic', projection='galactic')
    add_plane(ax, color='g', plane='supergalactic', projection='galactic')
  
  if projection == 'Supergalactic':
    add_plane(ax, color='black', plane='equatorial', projection='supergalactic')
    add_plane(ax, color='r', plane='ecliptic', projection='supergalactic')
    add_plane(ax, color='b', plane='galactic', projection='supergalactic')  
    
  if projection == 'Ecliptic':
    add_plane(ax, color='black', plane='equatorial', projection='ecliptic')
    add_plane(ax, color='g', plane='supergalactic', projection='ecliptic')
    add_plane(ax, color='b', plane='galactic', projection='ecliptic')      
################################################################ 

  plt.show()









