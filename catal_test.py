





from math import *
from matplotlib.widgets import Cursor

import sys
import os
import subprocess
import numpy as np
from datetime import *
import time
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column 
import random
from kapteyn import wcs
from optparse import OptionParser



  
#################################################################

if __name__ == '__main__':
  
  file = 'all.iter.0.v33.group'
  
  
  
  mytable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )

  id         = mytable['pgc']
  flag       = mytable['flag']
  sgl        = mytable['sgl']
  sgb        = mytable['sgb']
  gl         = mytable['gl']
  gb         = mytable['gb']  
  ra         = mytable['ra']
  dec        = mytable['dec']
  
  N = len(id)
  NoGroups = len(id[np.where(flag==2)])
  i = 0 
  if NoGroups!=0:
     while flag[i] != 2:
        i+=1
        
  gr = flag[i]
  while gr == 2:
    
    gr_id = id[i]
    p = i
    i+=1
    while i<N and (flag[i]==1): 
      if gb[i]>0 and gb[p]<0:
	print gr_id
      if gb[i]<0 and gb[p]>0:
	print gr_id      
      i+=1
    if i<N and flag[i]==2: 
      gr = 2
    else:
      break      
  
  
  
  
  
  
  