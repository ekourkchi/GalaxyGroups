

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

import spherematch as spherematch
from NAM_dFinder import *
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

def search(pgc_table, pgc_id, p, q):
  
  if p > q:
    return None
  m = (p+q)/2
  
  if pgc_table[m] == pgc_id:
    return m
  elif pgc_table[m] < pgc_id:
    return search(pgc_table, pgc_id, m+1, q)
  else:
    return search(pgc_table, pgc_id, p, m-1)
  
#################################################################
def find_dist(NEST_table, pgc_id):
  
  dist = NEST_table['dist']
  pgc  = NEST_table['PGC1']
  
  
  m =  search(pgc, pgc_id, 0, len(pgc)-1)
  
  if m != None:
    return [m, dist[m]]
  else:
    return None
 
    
    



#################################################################

if __name__ == '__main__':
  
  nest_file = 'nests_VI38.1799extended_lob.csv'
  nam_file  = 'NAMtestV5_08.dat.csv'
  NEST_table = nest(nest_file, nam_file)
  
  print find_dist(NEST_table, 27777)
  
  catal = 'north.iter.2.v31.group'
  mytable = np.genfromtxt(catal , delimiter='|', filling_values="-100000", names=True, dtype=None )
  
  id   = mytable['pgc']
  flag = mytable['flag']
  
  N = len(id)
  NoGroups = len(id[np.where(flag==2)])
  print "Number of groups: ", NoGroups
  if NoGroups == 0:
    print "[Warning] No group found in data-base ..." 
    print "Check the input catalog and choose the right option ...\n" 
    sys.exit()
  
  i = 0 
  if NoGroups!=0:
    while flag[i] != 2:
      i+=1
      
  
  gr = flag[i]
  while gr == 2:
    
    grp_id =  id[i]
    dist = -10000
    i+=1
    j = 0
    while i<N and flag[i]==1: 
      
      
      if dist == -10000:
	result = find_dist(NEST_table, id[i])
	
	if result != None:
	  dist = result[1]
	  if j != 0 : print "     ", j
      i+=1
      j+=1
      
    
    print grp_id, dist
    if i<N and flag[i]==2: 
      gr = 2
    else:
       break
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



