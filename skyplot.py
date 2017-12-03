#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2016"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "v2.6"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"

## Importing Important Python Libraries

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


import ttk as ttk
import Tkinter as tk
from tkFileDialog   import askopenfilename
import tkMessageBox
import tkFileDialog

import matplotlib.backends.backend_tkagg as tkagg



#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#if sys.version_info[0] < 3:
    #import Tkinter as Tk
#else:
    #import tkinter as Tk
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

event_active = True

virgo_on = False

LARGE_FONT= ("Verdana", 12)

col_pallet = ['darkviolet', 'blue', 'deepskyblue', 'forestgreen', 'y', 'gold', 'darkorange', 'red', 'magenta', 'maroon', 'sienna', 'slategrey', 'black']
vel_pallet = [-100000, 0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 100000]



#################################################
def zoom(xy_lim, ratio=1):
     

     i1=xy_lim[0]; i2=xy_lim[1]
     j1=xy_lim[2]; j2=xy_lim[3]
     
 
     xc = 0.5*(i1+i2)
     yc = 0.5*(j1+j2)
     delta_x = abs(i2-i1); delta_y =  abs(j2-j1)
     dx = 0.5*ratio*delta_x; dy=0.5*ratio*delta_y
     
     if i2 < i1:
       i1 = xc+dx ; i2=xc-dx
     else:
       i1 = xc-dx ; i2=xc+dx
     
     if j2 < j1:
       j1 = yc+dy ; j2=yc-dy
     else:
       j1 = yc-dy ; j2=yc+dy 
     
     xy_lim = [i1,i2,j1,j2]
     
     return xy_lim
   
#################################################
def pan(xy_lim, xc, yc, ratio=1):
     
     i1=xy_lim[0]; i2=xy_lim[1]
     j1=xy_lim[2]; j2=xy_lim[3]

     delta_x = abs(i2-i1); delta_y =  abs(j2-j1)
     dx = 0.5*ratio*delta_x; dy=0.5*ratio*delta_y

     if i2 < i1:
       i1 = xc+dx ; i2=xc-dx
     else:
       i1 = xc-dx ; i2=xc+dx
     
     if j2 < j1:
       j1 = yc+dy ; j2=yc-dy
     else:
       j1 = yc-dy ; j2=yc+dy 
     
     xy_lim = [i1,i2,j1,j2]
     
     return  xy_lim
#################################################
class GroupBrowser(object):
    
    def __init__(self, A, B, table_gal, A_cluster, B_cluster, table_all, root_win, mode):
      

	
	id_gal       = table_gal['pgc']#
	flag_gal     = table_gal['flag']#
	Vls_gal      = table_gal['Vls']#
	R_2t2_gal    = table_gal['R2t_lum']#
	mDist_gal    = table_gal['mDist']#
	

	
	self.ra = A_cluster
	self.dec = B_cluster
	self.listplot = []	
	
	
	self.listplot_gal = []
	indices = np.where(flag_gal<2)
	self.ra_gal = A[indices]
	self.dec_gal = B[indices]
	self.R_2t2_gal = R_2t2_gal[indices]
	self.Vls_gal = Vls_gal[indices]
	self.id_gal = id_gal[indices]
	self.mDist_gal = mDist_gal[indices]
	
	self.root = root_win
	self.visible_list = None
	self.visible_gal = None
	
	id       = table_all['pgc']
	flag     = table_all['flag']
	sgl      = table_all['sgl']
	sgb      = table_all['sgb']
	gl       = table_all['gl']
	gb       = table_all['gb']  
	ra       = table_all['ra']
	dec      = table_all['dec']
	Ks       = table_all['Ks']
	Vls      = table_all['Vls']
	R_2t2    = table_all['R2t_lum']
	nest     = table_all['nest']
	dcf2     = table_all['dcf2']
	ed       = table_all['ed']
	objname  = table_all['objname']
	mDist    = table_all['mDist']
	mDistErr = table_all['mDistErr']
	sigmaP_lum = table_all['sigmaP_lum']

	

	
	d_theta = 0.01
	u = np.arange(0,2*pi,d_theta)
	

	for p in range(len(self.id_gal)):
	   Dist = self.Vls_gal[p]/H0
	   if Dist<1: Dist = 1.
	   if self.mDist_gal[p] != 0: Dist = self.mDist_gal[p]
	   
	   if mode!='cartesian':
	     r = 180.*atan(self.R_2t2_gal[p]/Dist)/pi
	   elif mode=='cartesian':
	     r = self.R_2t2_gal[p]
	   
	   
	   X = u*0.; Y = u*0.
	   RA0 = self.ra_gal[p]; DEC0 = self.dec_gal[p]
	   for q in range(0,len(u)):
		x = r*np.cos(u[q]) + RA0
		y = r*np.sin(u[q]) + DEC0
		X[q] = x
		Y[q] = y
	   line, = plt.plot(X,Y, '-.', markersize = 1, color='red', picker=False, visible=False)
	   line.set_dashes([8, 3]) 
	   self.listplot_gal.append([self.id_gal[p], line])
	

	
	N = len(id)
	NoGroups = len(id[np.where(flag==2)])
	i = 0 
        if NoGroups!=0:
         while flag[i] != 2:
           i+=1
        
        gr = flag[i]
        while gr == 2:
	      
	      nest_center = nest[i]
	      random.seed(nest[i])
	      Red, Green, Blue = random.random(), random.random(), random.random()
	      Dist = Vls[i]/H0
	      if Dist<1: Dist = 1.
	      
	      if mDist[i] != 0: Dist = mDist[i]
	      
	      if mode!='cartesian':
	        r = 180.*atan(R_2t2[i]/Dist)/pi
	      elif mode=='cartesian':
		r = R_2t2[i]
	      
	      X = u*0.; Y = u*0.
	      
	      RA0 = self.ra[i]
	      DEC0 = self.dec[i]
	      
	      for q in range(0,len(u)):
		x = r*np.cos(u[q]) + RA0
		y = r*np.sin(u[q]) + DEC0
		X[q] = x
		Y[q] = y
		
	      line, = plt.plot(X,Y, ':', markersize = 1, color='red', picker=False, visible=False)
	      Vls_list = []
	      Vls_list.append([Vls[i],sigmaP_lum[i], mDist[i], mDistErr[i]])

	      i+=1  
	      X = []
	      Y = []
	      while i<N and flag[i]==1: 
		RA0 = self.ra[i]
		DEC0 = self.dec[i]
		X.append(RA0)
		Y.append(DEC0)
		Vls_list.append(Vls[i])
		i+=1
	      group_plt = plt.plot(X, Y, 'o', markersize = 3, mfc='none', markeredgecolor = 'red', visible=False)
	      
	      
	      if line != None and group_plt!= None:
		#line.set_visible(True)
		#for galplot in group_plt:
		   #galplot.set_visible(True)
	        plot_list = [nest_center, line, group_plt, Vls_list]
	        self.listplot.append(plot_list)
	      
	      
	      if i<N and flag[i]==2: 
		gr = 2
	      else:
		break
    
    
    def set_visible_gal(self, pgc_no):  # This shows the R_2t of the galaxy
      if self.listplot_gal != None:
	for plot_list in self.listplot_gal:
	  if plot_list[0] == pgc_no:
	    self.visible_gal = plot_list
	    plot_list[1].set_visible(True)
	    break


    def set_invisible_gal(self):  # 
      if self.visible_gal != None:
	 self.visible_gal[1].set_visible(False)
	 self.visible_gal = None


    ######
    def set_visible(self, pgc_no):  # that's the nest pgc and all of its galaxies
      
      if self.listplot != None:
	for plot_list in self.listplot:
	  if plot_list[0] == pgc_no:
	    self.visible_list = plot_list
	    plot_list[1].set_visible(True)
	    for galplot in plot_list[2]:
	      galplot.set_visible(True)
	    self.root.onClick([pgc_no, plot_list[3]])
	    
	    break



    def set_invisible(self):  # that's the nest pgc
      
      if self.visible_list != None:
	  self.visible_list[1].set_visible(False)
	  for galplot in self.visible_list[2]:
	     galplot.set_visible(False)
	     self.visible_list = None

      
#################################################

class PointBrowser(object):
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """

    def __init__(self, A, B, table_gal, A_cluster, B_cluster, table_all, root_win, mode):
        self.lastind = 0

        self.selected, = plt.plot(A, B, '*', ms=6, alpha=0.8,
                                 color='red', markeredgecolor='red', visible=False)
	
	self.group_plot = GroupBrowser(A, B, table_gal, A_cluster, B_cluster, table_all, root_win, mode)
	
	
	self.id_gal       = table_gal['pgc']
	self.flag_gal     = table_gal['flag']
	self.sgl_gal      = table_gal['sgl']
	self.sgb_gal      = table_gal['sgb']
	self.gl_gal       = table_gal['gl']
	self.gb_gal       = table_gal['gb']  
	self.ra_gal       = table_gal['ra']
	self.dec_gal      = table_gal['dec']
	self.Ks_gal       = table_gal['Ks']
	self.Vls_gal      = table_gal['Vls']
	self.R_2t2_gal    = table_gal['R2t_lum']
	self.nest_gal     = table_gal['nest']
	self.dcf2_gal     = table_gal['dcf2']
	self.ed_gal       = table_gal['ed']
	self.objname_gal  = table_gal['objname']
	self.mDist_gal    = table_gal['mDist']
	self.mDistErr_gal = table_gal['mDistErr']
	
	self.A = A
	self.B = B
	
	self.line = root_win.line
	self.fig  = root_win.fig
	self.ax   = root_win.ax
	self.mode = mode
	
	self.Usgx = None
	self.Usgy = None
	self.Usgz = None
	
	if mode == 'cartesian':
	  self.Usgx = self.ax.annotate(" ", (0.016,0.53), xycoords='figure fraction', size=11)
	  self.Usgy = self.ax.annotate(" ", (0.016,0.48), xycoords='figure fraction', size=11)
	  self.Usgz = self.ax.annotate(" ", (0.016,0.43), xycoords='figure fraction', size=11)

    def onpress(self, event):
        if self.lastind is None:
            return
        if event.key not in ('n', 'p'):
            return
        if event.key == 'n':
            inc = 1
        else:
            inc = -1

        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(A) - 1)
        self.update()

    def onpick(self, event):
      
      global event_active
      if not event_active: return
      
      if event.mouseevent.button == 1 :   # It's just active at left mouse click
        if event.artist != self.line:
            return True

        N = len(event.ind)
        if not N:
            return True

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        distances = np.hypot(x - self.A[event.ind], y - self.B[event.ind])
        indmin = distances.argmin()
        dataind = event.ind[indmin]
        
        if N>1: 
	  ind = event.ind
	  print "\n\n*****************************************************"
	  print "***** All Neighbour Galaxies ************************"
	  print "*****************************************************"

	  myTable = Table()
	  myTable.add_column(Column(data=self.id_gal[ind], name='pgc'))
	  myTable.add_column(Column(data=self.flag_gal[ind], name='flag'))
	  myTable.add_column(Column(data=self.ra_gal[ind],  name='ra'))  
	  myTable.add_column(Column(data=self.dec_gal[ind],  name='dec')) 
	  #myTable.add_column(Column(data=self.gl_gal[ind],  name='gl')) 
	  #myTable.add_column(Column(data=self.gb_gal[ind],  name='gb'))    
	  myTable.add_column(Column(data=self.sgl_gal[ind],  name='sgl')) 
	  myTable.add_column(Column(data=self.sgb_gal[ind],  name='sgb')) 
	  myTable.add_column(Column(data=self.Ks_gal[ind],  name='Ks'))
	  myTable.add_column(Column(data=self.Vls_gal[ind],  name='Vls')) 
	  myTable.add_column(Column(data=self.dcf2_gal[ind],  name='dcf2')) 
	  myTable.add_column(Column(data=self.ed_gal[ind],  name='ed')) 
	  myTable.add_column(Column(data=self.mDist_gal[ind],  name='mDist')) 
	  myTable.add_column(Column(data=self.mDistErr_gal[ind],  name='mDistErr')) 
	  myTable.add_column(Column(data=self.nest_gal[ind],  name='nest')) 
	  myTable.add_column(Column(data=self.objname_gal[ind],  name='objname')) 
	  print myTable
	  
        self.single_info(dataind)
    
    #### Also it's activated when we find a pgc galaxies
    def single_info(self, dataind):
      
        ind = [dataind]
        print "\n\n------ Selected Galaxy -------------------------------"

        myTable = Table()
        myTable.add_column(Column(data=self.id_gal[ind], name='pgc'))
        myTable.add_column(Column(data=self.flag_gal[ind], name='flag'))
        myTable.add_column(Column(data=self.ra_gal[ind],  name='ra'))  
        myTable.add_column(Column(data=self.dec_gal[ind],  name='dec')) 
        #myTable.add_column(Column(data=self.gl_gal[ind],  name='gl')) 
        #myTable.add_column(Column(data=self.gb_gal[ind],  name='gb'))    
        myTable.add_column(Column(data=self.sgl_gal[ind],  name='sgl')) 
        myTable.add_column(Column(data=self.sgb_gal[ind],  name='sgb'))
        myTable.add_column(Column(data=self.Ks_gal[ind],  name='Ks'))
        myTable.add_column(Column(data=self.Vls_gal[ind],  name='Vls')) 
        myTable.add_column(Column(data=self.dcf2_gal[ind],  name='dcf2')) 
        myTable.add_column(Column(data=self.ed_gal[ind],  name='ed')) 
        myTable.add_column(Column(data=self.mDist_gal[ind],  name='mDist')) 
        myTable.add_column(Column(data=self.mDistErr_gal[ind],  name='mDistErr')) 
        myTable.add_column(Column(data=self.nest_gal[ind],  name='nest')) 
        myTable.add_column(Column(data=self.objname_gal[ind],  name='objname')) 
        print myTable        
        
        if self.mode=='cartesian':
	  sgx, sgy, sgz = xyz_list(self.sgl_gal[ind], self.sgb_gal[ind], self.Vls_gal[ind], self.dcf2_gal[ind], self.mDist_gal[ind], self.flag_gal[ind])
	  self.Usgx.set_text("[SGX]: "+'{:.4f}'.format(sgx[0])+ ' Mpc')
	  self.Usgy.set_text("[SGY]: "+'{:.4f}'.format(sgy[0])+ ' Mpc')
	  self.Usgz.set_text("[SGZ]: "+'{:.4f}'.format(sgz[0])+ ' Mpc')
        

        self.lastind = dataind
        self.update()
        
        return self.A[dataind], self.B[dataind]

    def update(self):
        if self.lastind is None:
            return

        dataind = self.lastind
        
        self.group_plot.set_invisible()
        
        
        self.group_plot.set_invisible_gal() #
        self.group_plot.set_visible_gal(self.id_gal[dataind]) #
        
        
        # set the desired point visible
        self.selected.set_visible(True)
        self.selected.set_data(self.A[dataind], self.B[dataind]) 

        
        self.group_plot.set_visible(self.nest_gal[dataind])        
        self.fig.canvas.draw()
        


    def unselect(self, event):
        

        if self.lastind is None:
            return
	  
        if event.dblclick:
          self.group_plot.set_invisible()
          self.group_plot.set_invisible_gal()
	  dataind = self.lastind
	  # set the desired point invisible
	  self.selected.set_visible(False)
	  self.selected.set_data(self.A[dataind], self.B[dataind]) 
	  
	  if self.mode=='cartesian':
	    self.Usgx.set_text(" ")
	    self.Usgy.set_text(" ")
	    self.Usgz.set_text(" ")
	  
	  self.fig.canvas.draw()
       
#################################################
def table_deliver(filee):
  
  

  try:
    mytable = np.genfromtxt(filee , delimiter='|', filling_values="-100000", names=True, dtype=None )
   
  except:
    print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)
  
  id         = mytable['pgc']
  flag       = mytable['flag']
  sgl        = mytable['sgl']
  sgb        = mytable['sgb']
  gl         = mytable['gl']
  gb         = mytable['gb']  
  ra         = mytable['ra']
  dec        = mytable['dec']
  Ks         = mytable['Ks']
  Vls        = mytable['Vls']
  R_2t2      = mytable['R2t_lum']
  nest       = mytable['nest']
  dcf2       = mytable['dcf2']
  ed         = mytable['ed']
  objname    = mytable['objname']   
  mDist      = mytable['mDist']
  mDistErr   = mytable['mDistErr']
  sigmaP_lum = mytable['sigmaP_lum']
  
  print "\n[Success] The catalog \""+filee+"\" was sucessfully loaded ..."
  print "and has "+str(len(id))+" entries ...\n"
  
  
  out_table = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'dcf2':dcf2, 'ed':ed, 'objname':objname, \
		'mDist':mDist, 'mDistErr':mDistErr, 'sigmaP_lum':sigmaP_lum}
  return out_table
#################################################
def file_deliver(filee):
  
  prefix = filee.split('.')[0] 
  
  if prefix == 'north' or prefix == 'south':
    return table_deliver(filee)
  
  Have_north = os.path.isfile('north.'+filee)
  Have_south = os.path.isfile('south.'+filee)
  
  if Have_north and not Have_south:
    print "\n[Warning] The catalog \""+'south.'+filee+"\" is not available, ..."
    return table_deliver('north.'+filee)
  
  elif Have_south and not Have_north:
    print "\n[Warning] The catalog \""+'north.'+filee+"\" is not available, ..."
    return table_deliver('south.'+filee)
  
  elif not Have_north and not Have_south:

    print "\n[Warning] The catalog \""+'south.'+filee+"\" is not available, ..."
    print "[Warning] The catalog \""+'north.'+filee+"\" is not available, ...\n"
    
    return table_deliver(filee)
   
  
  else:
    north_table = table_deliver('north.'+filee)
    south_table = table_deliver('south.'+filee)
    
    north_flag = north_table['flag']
    south_flag = south_table['flag']
    
    north_ind = np.where(north_flag>0)
    south_ind = np.where(south_flag>0)
    
    id = np.concatenate((north_table['pgc'][north_ind],south_table['pgc'][south_ind]))
    flag = np.concatenate((north_table['flag'][north_ind],south_table['flag'][south_ind]))
    sgl = np.concatenate((north_table['sgl'][north_ind],south_table['sgl'][south_ind]))    
    sgb = np.concatenate((north_table['sgb'][north_ind],south_table['sgb'][south_ind]))
    gl = np.concatenate((north_table['gl'][north_ind],south_table['gl'][south_ind]))
    gb = np.concatenate((north_table['gb'][north_ind],south_table['gb'][south_ind]))
    ra = np.concatenate((north_table['ra'][north_ind],south_table['ra'][south_ind]))
    dec = np.concatenate((north_table['dec'][north_ind],south_table['dec'][south_ind]))
    Ks = np.concatenate((north_table['Ks'][north_ind],south_table['Ks'][south_ind]))
    Vls = np.concatenate((north_table['Vls'][north_ind],south_table['Vls'][south_ind]))
    R_2t2 = np.concatenate((north_table['R2t_lum'][north_ind],south_table['R2t_lum'][south_ind]))
    nest = np.concatenate((north_table['nest'][north_ind],south_table['nest'][south_ind]))
    dcf2 = np.concatenate((north_table['dcf2'][north_ind],south_table['dcf2'][south_ind]))
    ed = np.concatenate((north_table['ed'][north_ind],south_table['ed'][south_ind]))
    objname = np.concatenate((north_table['objname'][north_ind],south_table['objname'][south_ind]))
    mDist = np.concatenate((north_table['mDist'][north_ind],south_table['mDist'][south_ind]))
    mDistErr = np.concatenate((north_table['mDistErr'][north_ind],south_table['mDistErr'][south_ind])) 
    sigmaP_lum = np.concatenate((north_table['sigmaP_lum'][north_ind],south_table['sigmaP_lum'][south_ind])) 
    
    north_ind = np.where(north_flag<=0)
    south_ind = np.where(south_flag<=0)
    
    
    id1 = np.concatenate((north_table['pgc'][north_ind],south_table['pgc'][south_ind]))
    flag1 = np.concatenate((north_table['flag'][north_ind],south_table['flag'][south_ind]))
    sgl1 = np.concatenate((north_table['sgl'][north_ind],south_table['sgl'][south_ind]))    
    sgb1 = np.concatenate((north_table['sgb'][north_ind],south_table['sgb'][south_ind]))
    gl1 = np.concatenate((north_table['gl'][north_ind],south_table['gl'][south_ind]))
    gb1 = np.concatenate((north_table['gb'][north_ind],south_table['gb'][south_ind]))
    ra1 = np.concatenate((north_table['ra'][north_ind],south_table['ra'][south_ind]))
    dec1 = np.concatenate((north_table['dec'][north_ind],south_table['dec'][south_ind]))
    Ks1 = np.concatenate((north_table['Ks'][north_ind],south_table['Ks'][south_ind]))
    Vls1 = np.concatenate((north_table['Vls'][north_ind],south_table['Vls'][south_ind]))
    R_2t21 = np.concatenate((north_table['R2t_lum'][north_ind],south_table['R2t_lum'][south_ind]))
    nest1 = np.concatenate((north_table['nest'][north_ind],south_table['nest'][south_ind]))
    dcf21 = np.concatenate((north_table['dcf2'][north_ind],south_table['dcf2'][south_ind]))
    ed1 = np.concatenate((north_table['ed'][north_ind],south_table['ed'][south_ind]))
    objname1 = np.concatenate((north_table['objname'][north_ind],south_table['objname'][south_ind]))
    mDist1 = np.concatenate((north_table['mDist'][north_ind],south_table['mDist'][south_ind]))
    mDistErr1 = np.concatenate((north_table['mDistErr'][north_ind],south_table['mDistErr'][south_ind]))
    sigmaP_lum1 = np.concatenate((north_table['sigmaP_lum'][north_ind],south_table['sigmaP_lum'][south_ind])) 
    
    id = np.concatenate((id, id1))
    flag = np.concatenate((flag, flag1))
    sgl = np.concatenate((sgl, sgl1))  
    sgb = np.concatenate((sgb, sgb1)) 
    gl = np.concatenate((gl, gl1)) 
    gb = np.concatenate((gb, gb1)) 
    ra = np.concatenate((ra, ra1)) 
    dec = np.concatenate((dec, dec1)) 
    Ks = np.concatenate((Ks, Ks1)) 
    Vls = np.concatenate((Vls, Vls1)) 
    R_2t2 = np.concatenate((R_2t2, R_2t21)) 
    nest = np.concatenate((nest, nest1)) 
    dcf2 = np.concatenate((dcf2, dcf21)) 
    ed = np.concatenate((ed, ed1)) 
    objname = np.concatenate((objname, objname1)) 
    mDist = np.concatenate((mDist, mDist1)) 
    mDistErr = np.concatenate((mDistErr, mDistErr1))
    sigmaP_lum = np.concatenate((sigmaP_lum, sigmaP_lum1))

    out_table = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'dcf2':dcf2, 'ed':ed, 'objname':objname, \
		 'mDist':mDist, 'mDistErr':mDistErr, 'sigmaP_lum':sigmaP_lum}    
    
    return out_table


#################################################

def arg_parser():
    parser = OptionParser(usage="""\
\n
 - A GUI for displaying spatial distribution of the groups ...
 - Single Click --> select an object
   -- if it's a single galaxy: shows its R_2t
   -- if it belongs to a group, highlights the entire group, and
       shows its physical R_2t using the actual measured distances,
       just to see how it is compared to that one just based on velocities
 - Double click --> unselect all objects
 

 - How to run: 
 
     1$ %prog -f [catalog_name] -c [coord_system] -a [alpha] -d [delta]
     
     or
     
     2$ %prog -f [catalog_name] -c [coord_system] -p [pgc_number]
     - It works, if the PGC object is available in the catalog


 - Use the -h option to see all options ...

 - "alpha" and "delta" are the right ascentation and declination of the 
   center of the projected coordinate system, in degree.
 
 - The specified center and its neighbouring points, is projected into 
   a new sphere, where the center lies on the equator with no declination
   
 - This way, a flat tangential projection is obtained, where the original 
   coordinate frame, (i.e. equatorial, galactic or supergalactic) is displayed
   by dotted lines.    
   
 - if coord_system is equatorial    -->  alpha=RA0  delta=DEC0
 - if coord_system is galactic      -->  alpha=l0   delta=b0 
 - if coord_system is supergalactic -->  alpha=SGL0 delta=SGB0


 - Example(s): 
    $ python %prog -f north.iter.2.v26.group -c supergalactic -a 50 -d 0
    
    If the [catalog_name] is not available, then 'north' and 'south' would be added as prefix
    and then it looks for any available catalog. 
    
    $ python %prog -f iter.2.v26.group -c supergalactic -a 360 -d 0 w 30
    - This would plot north and south (Galactic) catalogs together, where SGL=0 or SGL=360
    
    $ python %prog -f iter.2.v26.group -c supergalactic -p  34695
    - This would center the plot on the group with PGC ID = 34695
        
    $ python %prog -h 
      To see help and all available options.
 
 - Author: "Ehsan Kourkchi"
 - Version: v2.6
 - Copyright 2016

""")
    
    parser.add_option('-a', '--alpha',
                      type='float', action='store',
                      help="""right ascension of the coordinate system center [degree]""")
    parser.add_option('-d', '--delta',
                      type='float', action='store',
                      help="""declination of the coordinate system center [degree]""")    
    parser.add_option('-f', '--file',
                      type='string', action='store',
                      help="""The input file, e.g. \'north.iter.9.v23.group\' """)
    parser.add_option("-c", "--coordinate",
                      type='string', action='store',
                      help="""Coordinate system (i.e. equatorial, galactic, supergalactic, cartesian)""")
    parser.add_option("-p", "--pgc",
                      type='float', action='store',
                      help="""pgc number (optional, specify the plot center)""")    
    parser.add_option("-w", "--width",
                      type='float', action='store',
                      help="""The width of the plot in degrees""")    
    #parser.add_option('-n', "--north",
                  #action="store_true", dest="north", default=False,
                  #help="Selecting the north galactic sky")
      

    parser.add_option("-m", "--Vmin",
                      type='float', action='store',
                      help="""The lower limit on radial velocity""") 
    
    parser.add_option("-x", "--Vmax",
                      type='float', action='store',
                      help="""The upper limit on radial velocity""")         
    

    parser.add_option("-X", "--SGX",
                      type='float', action='store',
                      help="""SGX center""")

    parser.add_option("-Y", "--SGY",
                      type='float', action='store',
                      help="""SGY center""")

    parser.add_option("-Z", "--SGZ",
                      type='float', action='store',
                      help="""SGZ center""")
    
    parser.add_option("-j", "--projection",
                      type='int', action='store',
                      help="""projection type
                      
                      1  - SGX-SGY
                      2  - SGX-SGZ
                      3  - SGY-SGZ
                      10 - SGY-SGX
                      20 - SGZ-SGX
                      30 - SGZ-SGY
                      """)     

    parser.add_option("-u", "--XT",
                      type='float', action='store',
                      help="""X-axis thickness (Mpc)""")
    
    parser.add_option("-v", "--YT",
                      type='float', action='store',
                      help="""Y-axis thickness (Mpc)""")

    parser.add_option("-t", "--ZT",
                      type='float', action='store',
                      help="""Z-axis thickness (Mpc)""")    
    
    
    parser.add_option("-G", action="store_true", help="Plot Virgo Border", default=False)
    
    (opts, args) = parser.parse_args()
    return opts, args
########


# The interface is in degree
def tangent(alpha, delta, V):
  
  alpha *= pi/180.
  delta *= pi/180.
  
  C_alph = cos(alpha)
  S_alph = sin(alpha)
  C_delt = cos(delta)
  S_delt = sin(delta)
  
  x1 = V[0]*C_alph + V[1]*S_alph
  y1 = V[1]*C_alph - V[0]*S_alph
  z1 = V[2]
  
  V[0] = x1*C_delt + z1*S_delt
  V[1] = y1
  V[2] = z1*C_delt - x1*S_delt
  
  return V

# **************************************

# The interface is in degree
def tangent_inv(alpha, delta, V):
  
  alpha *= -pi/180.
  delta *= -pi/180.
  
  C_alph = cos(alpha)
  S_alph = sin(alpha)
  C_delt = cos(delta)
  S_delt = sin(delta)
  
  
  x1 = V[0]*C_delt + V[2]*S_delt
  y1 = V[1]
  z1 = V[2]*C_delt - V[0]*S_delt
  
  
  V[0] = x1*C_alph + y1*S_alph
  V[1] = y1*C_alph - x1*S_alph
  V[2] = z1  
  
  return V

# **************************************
def xyz_list(alpha_list, delta_list, Vls_list, dcf2_list, mDist_list, flag_list):
  
  N = len(alpha_list)
  sgx = np.zeros(N)
  sgy = np.zeros(N)
  sgz = np.zeros(N)
  for i in range(N):

    d=Vls_list[i]/H0
    if d<=0:
      d = 1
    if mDist_list[i] != 0:
      d = mDist_list[i]
    if dcf2_list[i] != 0:
      d = dcf2_list[i]
      
    
    V = d*xyz(alpha_list[i], delta_list[i])
    sgx[i] = V[0]
    sgy[i] = V[1]
    sgz[i] = V[2]
	
  return sgx, sgy, sgz

# **************************************
def xyz(alpha, delta):
  
  alpha *= pi/180.
  delta *= pi/180.
  
  C_alph = cos(alpha)
  S_alph = sin(alpha)
  C_delt = cos(delta)
  S_delt = sin(delta)
  
  V = [C_delt*C_alph, C_delt*S_alph, S_delt]
  return np.asarray(V)


def spherical(V):
  
  alpha = atan2(V[1], V[0])
  delta = asin(V[2])
  
  return alpha*180./pi, delta*180./pi
  

# **************************************
# returns angular separation of 
# **************************************
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

   
   if XdotY > 1 :
     theta12 = 0.
   elif XdotY < -1 :
     theta12 = -1.*pi
   else:
     theta12 = acos(XdotY)  
   return theta12*180./pi 


#################################################################
def transform_inv(alpha, delta, a, b):
  

  
  xvz = tangent_inv(alpha, delta, xyz(a, b))
  ra, dec = spherical(xvz)
  
  return ra, dec

#################################################################
def transform(alpha, delta, ra, dec):
  

  
  
  xvz = tangent(alpha, delta, xyz(ra, dec))
  a, b = spherical(xvz)
  
  return a, b

#################################################################
def transform_list(alpha, delta, ra, dec):
  
  
  if len(ra) != len(dec): 
    print "[Error] RA and DEC should be in the same size ..."
    return 0, 0 
  
  
  N = len(ra)
  A = np.zeros(N)
  B = np.zeros(N)
  
  for i in range(len(ra)):
    A[i], B[i] = transform(alpha, delta, ra[i], dec[i])
      
  return A, B

#################################################################

def main(alpha, delta, fig, ax, width=10, mode='equatorial', V_range=None, SGX=None, SGY=None, SGZ=None, Projection=None):
  global i1, i2, j1, j2
  
  R = 2*width
  step = 0.5

  ax.annotate("Center *", (0.016,0.9), xycoords='figure fraction', size=14, color='maroon')
  
  if mode=="equatorial":
     x0_text = 'RA0: '
     y0_text = 'DEC0: '
  elif mode=="galactic":
     x0_text = 'l0: '
     y0_text = 'b0: '     
  elif mode=="supergalactic":
     x0_text = 'SGL0: '
     y0_text = 'SGB0: '     
  elif mode=="cartesian":
    if Projection == 1:
      x0_text = 'SGX: '
      y0_text = 'SGY: '      
    elif Projection == 2:
      x0_text = 'SGX: '
      y0_text = 'SGZ: ' 
    elif Projection == 3:
      x0_text = 'SGY: '
      y0_text = 'SGZ: ' 
    if Projection == 10:
      x0_text = 'SGY: '
      y0_text = 'SGX: '      
    elif Projection == 20:
      x0_text = 'SGZ: '
      y0_text = 'SGX: ' 
    elif Projection == 30:
      x0_text = 'SGZ: '
      y0_text = 'SGY: ' 
      
  if mode!="cartesian":
    ax.annotate(x0_text+'{:.4f}'.format(alpha), (0.016,0.86), xycoords='figure fraction', size=11)
    ax.annotate(y0_text+'{:.4f}'.format(delta), (0.016,0.82), xycoords='figure fraction', size=11)
  elif mode=="cartesian":
    if SGX!=None and SGY!=None and Projection==1:
      ax.annotate(x0_text+'{:.4f}'.format(SGX), (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+'{:.4f}'.format(SGY), (0.016,0.82), xycoords='figure fraction', size=11)
    elif SGX!=None and SGZ!=None and Projection==2:
      ax.annotate(x0_text+'{:.4f}'.format(SGX), (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+'{:.4f}'.format(SGZ), (0.016,0.82), xycoords='figure fraction', size=11)  
    elif SGY!=None and SGZ!=None and Projection==3:
      ax.annotate(x0_text+'{:.4f}'.format(SGY), (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+'{:.4f}'.format(SGZ), (0.016,0.82), xycoords='figure fraction', size=11)        
    elif SGX!=None and SGY!=None and Projection==10:
      ax.annotate(x0_text+'{:.4f}'.format(SGY), (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+'{:.4f}'.format(SGX), (0.016,0.82), xycoords='figure fraction', size=11)
    elif SGX!=None and SGZ!=None and Projection==20:
      ax.annotate(x0_text+'{:.4f}'.format(SGZ), (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+'{:.4f}'.format(SGX), (0.016,0.82), xycoords='figure fraction', size=11)  
    elif SGY!=None and SGZ!=None and Projection==30:
      ax.annotate(x0_text+'{:.4f}'.format(SGZ), (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+'{:.4f}'.format(SGX), (0.016,0.82), xycoords='figure fraction', size=11)         
    else:
      ax.annotate(x0_text+' 0', (0.016,0.86), xycoords='figure fraction', size=11)
      ax.annotate(y0_text+' 0'.format(SGX), (0.016,0.82), xycoords='figure fraction', size=11)         
 
  ax.annotate("Coordinates ", (0.016,0.73), xycoords='figure fraction', size=14, color='maroon')
   
  Ux = ax.annotate(" ", (0.016,0.69), xycoords='figure fraction', size=12)
  Uy = ax.annotate(" ", (0.016,0.65), xycoords='figure fraction', size=12)
  
  Ura = ax.annotate(" ", (0.016,0.56), xycoords='figure fraction', size=11)
  Ualf = ax.annotate(" ", (0.016,0.54), xycoords='figure fraction', size=10)
  Udec = ax.annotate(" ", (0.016,0.50), xycoords='figure fraction', size=11)
  Udelt = ax.annotate(" ", (0.016,0.47), xycoords='figure fraction', size=10)  
      
  Ugl = ax.annotate(" ", (0.016,0.38), xycoords='figure fraction', size=11)
  Ugb = ax.annotate(" ", (0.016,0.34), xycoords='figure fraction', size=11)
      
  Usgl = ax.annotate(" ", (0.016,0.25), xycoords='figure fraction', size=11)
  Usgb = ax.annotate(" ", (0.016,0.21), xycoords='figure fraction', size=11)     
  
  if mode!='cartesian':
      ax.annotate("WCS ...", (0.016,0.60), xycoords='figure fraction', size=14, color='maroon')
      ax.annotate("Galactic", (0.016,0.42), xycoords='figure fraction', size=14, color='maroon')
      ax.annotate("Super-Galactic", (0.016,0.29), xycoords='figure fraction', size=14, color='maroon')
  elif mode=='cartesian':
      ax.annotate("Selected Gaalxy", (0.016,0.58), xycoords='figure fraction', size=14, color='maroon')


  N =  np.int(np.round(alpha))
  N = 5*(N/5)
  
  ra_min = ceil(np.round(alpha)-R/cos(delta*pi/180))
  ra_max = floor(np.round(alpha)+R/cos(delta*pi/180))
  ra_list = np.concatenate((np.arange(N-5, ra_min-5, -5), np.arange(N, ra_max+5, 5)))
  ra_list_max = max(ra_list)
  ra_list_min = min(ra_list)
  
  
  M =  np.int(np.round(delta))
  M = 5*(M/5)  
  
  dec_min = M-R
  dec_max = M+R

#--------------------------------------------------------
  if mode != "cartesian":
      # grid-Vertical
      for ra in ra_list:
	A=[]; B=[]
	for dec in np.arange(dec_min, dec_max+step, step):
	    a, b = transform(alpha, delta, ra, dec)
	    A.append(a)
	    B.append(b)
	plt.plot(A, B, ':', color='#696969')

      

      # grid-Horiontal
      for dec in np.arange(dec_min, dec_max+5, 5):
	A=[]; B=[]
	for ra in np.arange(ra_min, ra_list_max+step, step):
	    a, b = transform(alpha, delta, ra, dec)
	    A.append(a)
	    B.append(b)
	plt.plot(A, B, ':', color='#696969')  

  
      ax.set_xlim(width,-1*width)
      ax.set_ylim(-1*width,width)
      i1 = width
      i2 = -1*width
      j1 = -1*width
      j2 = width

      xvz = tangent(alpha, delta, xyz(alpha, delta))
      a, b = spherical(xvz)
      plt.plot([a], [b], '*', markersize = 10, color= 'black')
#--------------------------------------------------------
  elif mode == "cartesian":
      
      # grid
      for pp in np.arange(-100,100,5):
	plt.plot([pp,pp], [-100,100], ':', color='#696969')
	plt.plot([-100,100], [pp,pp], ':', color='#696969')
      
      if SGX!=None and SGY!=None and Projection==1:
	i1 = SGX-width
	i2 = SGX+width
	j1 = SGY-width
	j2 = SGY+width
      elif SGX!=None and SGZ!=None and Projection==2:
	i1 = SGX-width
	i2 = SGX+width
	j1 = SGZ-width
	j2 = SGZ+width
      elif SGY!=None and SGZ!=None and Projection==3:
	i1 = SGY-width
	i2 = SGY+width
	j1 = SGZ-width
	j2 = SGZ+width	
      elif SGX!=None and SGY!=None and Projection==10:
	j1 = SGX-width
	j2 = SGX+width
	i1 = SGY-width
	i2 = SGY+width
      elif SGX!=None and SGZ!=None and Projection==20:
	j1 = SGX-width
	j2 = SGX+width
	i1 = SGZ-width
	i2 = SGZ+width
      elif SGY!=None and SGZ!=None and Projection==30:
	j1 = SGY-width
	j2 = SGY+width
	i1 = SGZ-width
	i2 = SGZ+width		
      elif V_range[1]!=None:
	i1 = -1.*V_range[1]/H0
	i2 = V_range[1]/H0
	j1 = -1.*V_range[1]/H0
	j2 = V_range[1]/H0
      else:
	print "\n\n[Error] Given corrdinates and the projection are not consistent ... "
	print "not enough information ..."
	print "Please check your input parameters !\n"
	try:
	   print 'The minimum requirement is the upper limit for radial velocity.'
	   print 'Either enter a number here to continue, or anything else to abort ...'
	   input_var = input("Enter V_max (km/s): ")
	   try:
	     V_range[1] = float(input_var)
	     i1 = -1.*V_range[1]/H0
	     i2 = V_range[1]/H0
	     j1 = -1.*V_range[1]/H0
	     j2 = V_range[1]/H0	     
	   except:
	     exit(1)
	except:
	   print '\nSomething went wrong, \nTry again with command line ...\n'
	   exit(1)
	
      ax.set_xlim(i1,i2)
      ax.set_ylim(j1,j2)
      plt.plot([0.5*(i1+i2)], 0.5*(j1+j2), '*', markersize = 10, color= 'black')
  
  ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    #bottom='off',      # ticks along the bottom edge are off
    top='off')         # ticks along the top edge are off
    ##labelbottom='off') # labels along the bottom edge are off
  
  ax.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    #bottom='off',      # ticks along the bottom edge are off
    right='off')         # ticks along the top edge are off
    ##labelbottom='off') # labels along the bottom edge are off  
  

  def in_motion(event):
     #print 'you pressed', event.key, event.button, event.xdata, event.ydata, event.key, event.inaxes
     a = event.xdata
     b = event.ydata
     if event.inaxes == ax: 
        
        if mode != "cartesian":
	  Ux.set_text("X: "+'{:.2f}'.format(a))
	  Uy.set_text("Y: "+'{:.2f}'.format(b))
	else:
	  if Projection == 1: 
	    Ux.set_text("SGX: "+'{:.2f}'.format(a))
	    Uy.set_text("SGY: "+'{:.2f}'.format(b))
	  if Projection == 2: 
	    Ux.set_text("SGX: "+'{:.2f}'.format(a))
	    Uy.set_text("SGZ: "+'{:.2f}'.format(b))
	  if Projection == 3: 
	    Ux.set_text("SGY: "+'{:.2f}'.format(a))
	    Uy.set_text("SGZ: "+'{:.2f}'.format(b))	  
	  if Projection == 10: 
	    Ux.set_text("SGY: "+'{:.2f}'.format(a))
	    Uy.set_text("SGX: "+'{:.2f}'.format(b))
	  if Projection == 20: 
	    Ux.set_text("SGZ: "+'{:.2f}'.format(a))
	    Uy.set_text("SGX: "+'{:.2f}'.format(b))
	  if Projection == 30: 
	    Ux.set_text("SGZ: "+'{:.2f}'.format(a))
	    Uy.set_text("SGY: "+'{:.2f}'.format(b))
	    
        if mode=="equatorial":
	  RA, DEC = transform_inv(alpha, delta, a, b)
          tran = wcs.Transformation("equatorial j2000 j2000", "galactic")
          gl, gb = tran((RA,DEC))
          tran = wcs.Transformation("equatorial j2000 j2000", "supergalactic")
          sgl, sgb = tran((RA,DEC))

	elif mode=="galactic":
	  gl, gb  = transform_inv(alpha, delta, a, b)
	  tran = wcs.Transformation("galactic j2000 j2000", "equatorial")
	  RA, DEC = tran((gl,gb))
	  tran = wcs.Transformation("galactic j2000 j2000", "supergalactic")
	  sgl, sgb = tran((gl,gb))
	  
	elif mode=="supergalactic":
	  sgl, sgb  = transform_inv(alpha, delta, a, b)
	  tran = wcs.Transformation("supergalactic j2000 j2000", "equatorial")
	  RA, DEC = tran((sgl,sgb))
	  tran = wcs.Transformation("supergalactic j2000 j2000", "galactic")
	  gl, gb = tran((sgl,sgb))


       
        if mode != "cartesian":
	  c = SkyCoord(ra=RA, dec=DEC, unit=(u.degree, u.degree))
	  wcs_hmsdms =  c.to_string('hmsdms', precision=4, sep=':')
	  wcs_hmsdms = wcs_hmsdms.split(" ")        
	  Ura.set_text("RA: "+'{:.4f}'.format(RA))
	  Udec.set_text("DEC: "+'{:.4f}'.format(DEC))  
	  Ualf.set_text(r"$\alpha: $"+wcs_hmsdms[0])
	  Udelt.set_text(r"$\delta: $"+wcs_hmsdms[1]) 
	  Ugl.set_text("l: "+'{:.4f}'.format(gl))
	  Ugb.set_text("b: "+'{:.4f}'.format(gb)) 
	  if sgl < 0: sgl += 360.
	  Usgl.set_text("SGL: "+'{:.4f}'.format(sgl))
	  Usgb.set_text("SGB: "+'{:.4f}'.format(sgb))

        
        fig.canvas.draw()
     else:
        Ux.set_text(" ")
        Uy.set_text(" ")
        Ura.set_text(" ")
        Udec.set_text(" ")
        Ualf.set_text(" ")
	Udelt.set_text(" ")
	Ugl.set_text(" ")
	Ugb.set_text(" ")	
	Usgl.set_text(" ")
	Usgb.set_text(" ")
        fig.canvas.draw()
   
   
  fig.canvas.mpl_connect('motion_notify_event', in_motion)
  
  xy_lim = [i1, i2, j1, j2]
  
  return ra_list_min, ra_list_max, dec_min, dec_max, xy_lim 
  
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
   

def ClusterPlot(id, ra, dec, Vls, R_2t2, mDist, nest, flag, parent, mode='none'):
  
  N = len(id)
  NoGroups = len(id[np.where(flag==2)])
  print "Number of groups: ", NoGroups
  
  if NoGroups == 0:
    print "[Warning] No group found in data-base ..." 
    print "Check the input catalog and choose the right option ...\n" 
    print "You may be using the wrong catalog, not covering the desired coordiante system/range ...\n" 
    exit(1)
    return
  
  
######################################################  
  if virgo_on:
    ### M87 = PGC 41361
    for p in range(len(id)):
      if id[p] == 41361: break
    
    Dist = Vls[p]/H0
    Vls_grp = Vls[p]
    if Dist<1: Dist = 1.
    if mode!='cartesian':
       r = 6.8
    elif mode=='cartesian':
       r = Dist*tan(6.8*pi/180.)
    d_theta = 0.01
    u = np.arange(0,2*pi,d_theta)
    X = u*0.
    Y = u*0.
    
    RA0 = ra[p]
    DEC0 = dec[p]
    
    for q in range(0,len(u)):
       x = r*np.cos(u[q]) + RA0
       y = r*np.sin(u[q]) + DEC0
       X[q] = x
       Y[q] = y       
    line, = plt.plot(X,Y, '-', markersize = 2, color='black', picker=False) 
    line.set_dashes([8, 3]) 
    
    
    if mode == "supergalactic":
      xA, yB = transform_list(parent.alpha, parent.delta, [110,105,105,96,0], [-3,-3,-4,-4,0])
      line, = plt.plot(xA[0:2],yB[0:2], '-', markersize = 2, color='black', picker=False) 
      line, = plt.plot(xA[1:3],yB[1:3], '-', markersize = 2, color='black', picker=False)
      line, = plt.plot(xA[2:4],yB[2:4], '-', markersize = 2, color='black', picker=False)

 
      
      
######################################################  

  
  i = 0 
  if NoGroups!=0:
    while flag[i] != 2:
      i+=1
  

  gr = flag[i]
  while gr == 2:
    
    random.seed(nest[i])
    Red, Green, Blue = random.random(), random.random(), random.random()
    colour = color_table(Vls[i])
    
    
    Dist = Vls[i]/H0
    Vls_grp = Vls[i]
    if Dist<1: Dist = 1.
    
    if mode!='cartesian':
       r = 180.*atan(R_2t2[i]/Dist)/pi
    elif mode=='cartesian':
       r = R_2t2[i]
    
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
        #line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Green, Blue), picker=False) 
        line, = plt.plot(X,Y, '-', markersize = 2, color=colour, picker=False) 
	line.set_dashes([8, 3]) 
	parent.list_line.append([line, (Red, Green, Blue), colour, Vls_grp, True])

    i+=1  
    X = []
    Y = []
    while i<N and flag[i]==1: 
      RA0 = ra[i]
      DEC0 = dec[i]
      X.append(RA0)
      Y.append(DEC0)
      i+=1
    #plt.plot(X, Y, 'o', markersize = 3, color=(Red, Green, Blue), markeredgecolor = (Red, Green, Blue))
    point, = plt.plot(X, Y, 'o', markersize = 3, color=colour, markeredgecolor = colour)
    parent.list_point.append([point, (Red, Green, Blue), colour, Vls_grp, True])
    
    if i<N and flag[i]==2: 
      gr = 2
    else:
       break
  
  

#################################################
FILEOPENOPTIONS = dict(defaultextension='.bin',
               filetypes=[('Group catalogs','*.group'), ('All files','*.*')])



 #################################
class control_help:

    def __init__(self, parent, text):
        top = self.top = tk.Toplevel(parent)
        self.myLabel = tk.Label(top, text='Help')
        self.myLabel.pack()
        
        self.parent = parent
        
        
        scrollbar_x = tk.Scrollbar(top, orient=tk.HORIZONTAL)
        scrollbar_y = tk.Scrollbar(top, orient=tk.VERTICAL)
        
        self.mySubmitButton = tk.Button(top, text='Close', command=self.quit) # top means in the current frame
        self.mySubmitButton.pack(side = "bottom")
        
       
        top.resizable(width=False, height=False)
        top.title('How it works?')
        
        T = tk.Listbox(top, height=50, width=80)
        
        T.config(xscrollcommand=scrollbar_x.set)
        T.config(yscrollcommand=scrollbar_y.set)
        

        scrollbar_x.config(command=T.xview)
        scrollbar_y.config(command=T.yview)
        scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X, expand=1)
        scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y, expand=1)
        

        T.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        lines =  text.split("\n")
        for line in lines:
           T.insert(tk.END, line)
        
       
        self.top.protocol("WM_DELETE_WINDOW", self.quit)
    
    def quit(self):
        self.parent.inputDialog = None
        self.top.destroy()
        
 #################################
 #################################
 #################################

class MainWin(tk.Tk):
    
    
    def OpenFile(self):
        name = askopenfilename(parent=self,**FILEOPENOPTIONS)
        if name!= None and name != "":
	   self.filee = name
	   if self.frame.inputDialog != None:
	      self.frame.inputDialog.quit()
           self.show_frame()
  
    def about(self):

	#print "This is \"skyplop\" version v2.6"
	#print "Copyright 2016"
	#print "Author: \"Ehsan Kourkchi\""
	text = '''
This is \"skyplop\" version v2.6 
Copyright 2016
Author: \"Ehsan Kourkchi\"   
'''
	tkMessageBox.showinfo("About", text)
    
    def controls(self):
      
      text = '''
      * Mouse Actions:
        
        1) Left single-click: selecting an object
        2) Left double-click on a wihite background: unselect
        3) Middle click: re-center the field
        4) Right double-click: re-load the catalog around the current mouse pointer 
                               position and with the same initial field parameters
        5) Middle button: scroll-up   : zoom-out
        6) Middle button: scroll-down : zoom-in
        
        
      ** Keyboard Actions:
        
        1) arrow-keys (up,down,left,right): re-centering the field
        2) ctrl+z: zoom-in
        3) ctrl+x: zoom-out
        
        
      *** General Description:
        
        - When an object is selected, if it belongs to a  group, the entire group 
          would be highlighted, and its real physical secound turnaround radius is 
          displayed by a dotted circle
          
        - The corresponding 2nd turnaround radius of the selected individual galaxy 
          would also be displayed
        
        - If a group is selected, its internal velocity distribution and other information
          would be displayed in a pop-up window. If user closes the window, it is displayed 
          any time a new group is chosed. If it is left open, by choosing a new group, its 
          information would be displayed in the same window
        
        - All corresponding information would also be dsipalyed in the terminal
        
      **** Menu-bar and miscellaneous
        
        - File > Open Catalog
          It would allow the user to choose a new catalog that covers the same region
          All other parameters would remain the same as they are entered in the terminal
          when starting the program
      
     
        - Use '-h' option with this code to see all available options and flags
         
         $ python skyplot -h 
          

      '''
      
      
      control_help(self, text)
      
      
    
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        
        self.filee = None
        self.mode = None
        self.alpha = None
        self.width = None
        self.pgc = None
        
        self.Vmin = None
        self.Vmax = None
        
        self.SGX = None
        self.SGY = None
        self.SGZ = None
        
        self.Projection = None
        
        self.xt = None
        self.yt = None
        self.zt = None
        
        self.container = tk.Frame(self, width=768, height=576)
        self.container.pack(side="top", fill="both", expand = True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        self.wm_title("skyplot")
        
        menu = tk.Menu(self)
        self.config(menu=menu)
        filemenu = tk.Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open Catalog ...", underline=0, command=self.OpenFile)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.on_exit)
        
        helpmenu = tk.Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)
        helpmenu.add_command(label="Skyplot Help", command=self.controls)
        helpmenu.add_command(label="About...", command=self.about)
        
        #filee = None
 
        self.protocol("WM_DELETE_WINDOW", self.on_exit)
        
    
    def on_exit(self):
        """When you click to exit, this function is called"""
        if tkMessageBox.askyesno("Exit", "Do you really want to quit skyplot?"):
	    if self.frame.inputDialog != None:
	      self.frame.inputDialog.quit()
            sys.exit()
            
    def show_frame(self):

      self.frame = Plot_page(self.container, self, self.filee, self.mode, self.alpha, self.delta, self.width, self.pgc, self.Vmin, self.Vmax, self.SGX, self.SGY, self.SGZ, self.Projection, self.xt, self.yt, self.zt)
      self.frame.grid(row=0, column=0, sticky="nsew")
      self.frame.tkraise()
      
      
 #################################
class VelDist:

    def __init__(self, parent, velocities):
        top = self.top = tk.Toplevel(parent)
        self.myLabel = tk.Label(top, text='Radial Velocity Distribution')
        self.myLabel.pack()
        #self.myEntryBox = tk.Entry(top)
        #self.myEntryBox.pack()
        top.resizable(width=False, height=False) # to force it to be the same size
        self.fig = plt.figure(figsize=(7,3), dpi=100)
        self.ax = self.fig.add_axes([0.1, 0.2, 0.85,  0.7]) 
        self.canvas = FigureCanvasTkAgg(self.fig, self.top)

        
        self.show_histogram(velocities)
        self.parent = parent
        
        
        self.mySubmitButton = tk.Button(top, text='Close', command=self.quit)
        self.mySubmitButton.pack(side = "bottom")
        
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self.top)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self.top.protocol("WM_DELETE_WINDOW", self.quit)
        
    def show_histogram(self, velocities, replot=False):
        
        vel      = velocities[1][0][0]
        sigma    = velocities[1][0][1]
        dist     = velocities[1][0][2]
        dist_e   = velocities[1][0][3]
        pgc_code = velocities[0]
        vel_data = velocities[1][1:]
        
        #print "pgc:", pgc_code
        #print "vel:", vel
        #print "sigma:", sigma
        #print "data:", vel_data
        bins = np.arange(vel-3*sigma, vel+3*sigma, sigma/2)
        
        if replot:
	  self.ax.cla()
	  self.fig.clf()
	  self.ax = self.fig.add_axes([0.1, 0.2, 0.85,  0.7]) 
	  
	  n, bins, patches = self.ax.hist(vel_data , bins, histtype='step',
			    color=['green'],
			    label=['label'], fill=False)   # , stacked=True)
	  
	  
	else:
          n, bins, patches = self.ax.hist(vel_data , bins, histtype='step',
			    color=['green'],
			    label=['label'], fill=False)   # , stacked=True)
        
        

        self.ax.minorticks_on()
        self.ax.tick_params(which='major', length=7)
        self.ax.tick_params(which='minor', length=4) 
        
        
        self.ax.set_xlabel(r"$V_{ls}$"+" [km/s]")
        self.ax.set_ylabel("Number of Galaxies")
        
        
        x = self.ax.get_xlim()
        y = self.ax.get_ylim()

        
        
        if y[1]+2 > 5:
	  self.ax.yaxis.set_major_locator(MultipleLocator(5))
          self.ax.yaxis.set_minor_locator(MultipleLocator(1))
        else:
	  self.ax.yaxis.set_major_locator(MultipleLocator(1))
	  self.ax.yaxis.set_minor_locator(MultipleLocator(1))
	  
        self.ax.set_ylim([y[0], y[1]+2])
        
        self.ax.plot([vel, vel],[y[0], y[1]+2], '--', color='black')
        self.ax.plot([vel-sigma, vel-sigma],[y[0], y[1]+2], ':', color='red')
        self.ax.plot([vel+sigma, vel+sigma],[y[0], y[1]+2], ':', color='red')
        
        
        self.ax.annotate("PGC: "+'{:d}'.format(pgc_code), (0.15,0.8), xycoords='figure fraction')
        self.ax.annotate(r"$Vls$"+": "+'{:d}'.format(vel)+"  " +r"$(km/s)$", (0.15,0.7), xycoords='figure fraction')
        self.ax.annotate(r"$\sigma_v$"+": "+'{:.0f}'.format(sigma), (0.15,0.6), xycoords='figure fraction')
        
        if dist != 0:
          self.ax.annotate(r"$d$"+":  "+'{:.2f}'.format(dist)+r'$\pm$'+'{:.0f}'.format(dist_e*100)+'% Mpc', (0.15,0.5), xycoords='figure fraction')
        else:
	  self.ax.annotate(r"$d$"+":  ?", (0.15,0.5), xycoords='figure fraction')
        
        
        if replot:
	  self.canvas.show()
          self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
          self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        
        
        #print velocities


    def quit(self):
        self.parent.inputDialog = None
        self.top.destroy()

#################################

class Plot_page(tk.Frame):
  
  def onClick(self, velocities):
    
        if self.inputDialog == None:
           self.inputDialog = VelDist(self, velocities)
           self.wait_window(self.inputDialog.top)
        else:
	   self.inputDialog.show_histogram(velocities, replot=True)
	   
        
  
  
  def on_exit(self):
        """When you click to exit, this function is called"""
        if tkMessageBox.askyesno("Exit", "Do you really want to quit skyplot?"):
	    if self.inputDialog != None:
	      self.inputDialog.quit()
            sys.exit()
            
  def __init__(self, parent, controller, filee, mode, alpha, delta, width, pgc, Vmin, Vmax, SGX, SGY, SGZ, Projection, xt, yt, zt):
        root = tk.Frame.__init__(self, parent)


        self.win = tk.Frame(self)
        #self.win.grid(row=1, column=5, columnspan=10)  
        

        self.win.pack(side='top')


        button = tk.Button(self, text="Quit", command=self.on_exit)# 
        button.pack(side='bottom')
       

        self.controller = controller
	self.fig = plt.figure(figsize=(13,8 ), dpi=100)
	
	resetax = axes([0.01, 0.18, 0.16, 0.77])
	info_bts = Button(resetax, '', color='lightgray', hovercolor='lightgray')
	
	resetax = axes([0.01, 0.05, 0.16, 0.1])
	info_bts = Button(resetax, '', color='lightgray', hovercolor='lightgray')
	
	
	self.color_bts = []
	self.color_labels = []
	self.list_line = []
	self.list_point = []
	self.lline = []
	
	self.Vmin = Vmin
	self.Vmax = Vmax
	self.SGX  = SGX
	self.SGY  = SGY
        self.SGZ  = SGZ
        self.Projection  = Projection
        self.xt = xt
        self.yt = yt
        self.zt = zt


	
	self.ax = self.fig.add_axes([0.25, 0.07, 7.2/13,  0.9])  # main plot area
	subplots_adjust(left=0.25, bottom=0.25)
	self.fig.patch.set_facecolor('lightgray')


	
	self.color_labels.append(self.ax.annotate(r"$\/\/\/\/\/\/\/\/\/\/ V_{ls} \/ < \/ 0   $", (0.88,0.91), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$ \/\/\/\/ 0 \/ \leq \/ V_{ls} \/ < \/ 250$', (0.88,0.91-0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$\/ 250 \/ \leq \/ V_{ls} \/ < \/ 500$', (0.88,0.91-2*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$\/ 500 \/ \leq \/ V_{ls} \/ < \/ 750$', (0.88,0.91-3*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$\/ 750 \/ \leq \/ V_{ls} \/ < \/ 1000$', (0.88,0.91-4*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$1000 \/ \leq \/ V_{ls} \/ < \/ 1250$', (0.88,0.91-5*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$1250 \/ \leq \/ V_{ls} \/ < \/ 1500$', (0.88,0.91-6*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$1500 \/ \leq \/ V_{ls} \/ < \/ 1750$', (0.88,0.91-7*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$1750 \/ \leq \/ V_{ls} \/ < \/ 2000$', (0.88,0.91-8*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$2000 \/ \leq \/ V_{ls} \/ < \/ 2500$', (0.88,0.91-9*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$2500 \/ \leq \/ V_{ls} \/ < \/ 3000$', (0.88,0.91-10*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$3000 \/ \leq \/ V_{ls} \/ < \/ 3500$', (0.88,0.91-11*0.05), xycoords='figure fraction', size=12, color='black'))
	self.color_labels.append(self.ax.annotate(r'$3500 \/ \leq \/ V_{ls}$', (0.88,0.91-12*0.05), xycoords='figure fraction', size=12, color='black'))
	

	self.canvas = FigureCanvasTkAgg(self.fig, self)
     
        self.inputDialog = None
        self.browser = None
        
        
        self.filee = filee
        self.mode = mode
        self.alpha = alpha
        self.delta = delta
        self.width = width
        self.pgc = pgc
        
        self.show_plot()  
        
        

        v = tk.StringVar()
        e = tk.Entry(self.win, textvariable=v)
        v.set("")
        
        v_min = tk.StringVar()
        e_min = tk.Entry(self.win, textvariable=v_min)
        
        
        v_max = tk.StringVar()
        e_max = tk.Entry(self.win, textvariable=v_max)
        
        
        if Vmin == None:
            self.label_vmin = self.ax.annotate(r"$V_{min}: - $", (0.01,0.97), xycoords='figure fraction', size=12, color='black')
            v_min.set("V.min : -")
            Velmin = -100000.
        else:
	    self.label_vmin = self.ax.annotate(r"$V_{min}: $"+str(Vmin), (0.01,0.97), xycoords='figure fraction', size=12, color='black')
	    v_min.set("V.min : "+str(Vmin))
	    Velmin = Vmin
        
        if Vmax == None:
            self.label_vmax = self.ax.annotate(r"$V_{max}: - $", (0.13,0.97), xycoords='figure fraction', size=12, color='black')
            v_max.set("V.max : -")
            Velmax = 100000.
        else:
	    self.label_vmax = self.ax.annotate(r"$V_{max}: $"+str(Vmax), (0.13,0.97), xycoords='figure fraction', size=12, color='black')
	    v_max.set("V.max : "+str(Vmax))
	    Velmax = Vmax


        for obj in self.list_line:
	  if obj[3] < Velmin or obj[3] >= Velmax: 
	    obj[4] = False 
	    obj[0].set_visible(False)
	  else: 
	      obj[4] = True
	      p = vel_indx(obj[3])
	      if self.color_bts[p][1]:
	        obj[0].set_visible(True)
        for obj in self.list_point:
	  if obj[3] < Velmin or obj[3] >= Velmax: 
	    obj[4] = False 
	    obj[0].set_visible(False)
	  else: 
	      obj[4] = True    
	      p = vel_indx(obj[3])
	      if self.color_bts[p][1]:
	        obj[0].set_visible(True)
        for obj in self.lline:
	  if obj[2] < Velmin or obj[2] >= Velmax: 
	    obj[3] = False 
	    obj[0].set_visible(False)
	  else: 
	      obj[3] = True    
	      p = vel_indx(obj[2])
	      if self.color_bts[p][1]:
	        obj[0].set_visible(True)
        


        def callback_vel():
            s_min = e_min.get()
            s_max = e_max.get()
            vel_min = []
            vel_max = []
            mn = False
            mx = False
            
            if s_min[-1] == '-':
	      vel_min = -100000.
	      mn = True
	    else: 
	      for t in s_min.split():
		try:
		  vel_min.append(float(t))
		  if len(vel_min)>0:
		    vel_min = int(vel_min[0])
		    v_min.set("V.min : "+str(vel_min))
		    mn = True
		except ValueError:
		  v_min.set("V.min : err")

            if s_max[-1] == '-':
	      vel_max = 100000.
	      mx = True
	    else: 
	      for t in s_max.split():
		try:
		  vel_max.append(float(t))
		  if len(vel_max)>0:
		    vel_max = int(vel_max[0])
		    v_max.set("V.max : "+str(vel_max))
		    mx = True
		except ValueError:
		  v_max.set("V.max : err")
                
                
            if mn and mx and vel_max<vel_min:
	      print "\nWarning:V.min must be smaller than V.max ...."
	      print "Try again ....\n"
	      return
	    if mn and mx:
	      if vel_min == -100000.:
	         self.label_vmin.set_text(r"$V_{min}:  -$")
	      else: 
		 self.label_vmin.set_text(r"$V_{min}:  $"+str(vel_min))
		 
	      if vel_max == 100000.:
	         self.label_vmax.set_text(r"$V_{max}:  -$")
	      else:
		 self.label_vmax.set_text(r"$V_{max}:  $"+str(vel_max))
	      for obj in self.list_line:
		if obj[3] < vel_min or obj[3] >= vel_max: 
		  obj[4] = False 
		  obj[0].set_visible(False)
		else: 
		    obj[4] = True
		    p = vel_indx(obj[3])
		    if self.color_bts[p][1]:
		      obj[0].set_visible(True)
	      for obj in self.list_point:
		if obj[3] < vel_min or obj[3] >= vel_max: 
		  obj[4] = False 
		  obj[0].set_visible(False)
		else: 
		    obj[4] = True    
		    p = vel_indx(obj[3])
		    if self.color_bts[p][1]:
		      obj[0].set_visible(True)
	      for obj in self.lline:
		if obj[2] < vel_min or obj[2] >= vel_max: 
		  obj[3] = False 
		  obj[0].set_visible(False)
		else: 
		    obj[3] = True    
		    p = vel_indx(obj[2])
		    if self.color_bts[p][1]:
		      obj[0].set_visible(True)
            
        def callclear_vel():
            v_min.set("V.min : ")
            v_max.set("V.max : ")

        def callshowall_vel():
	      if Vmin == None: 
	         v_min.set("V.min : -")
	      else: 
		 v_min.set("V.min : "+str(Vmin))
	      if Vmax == None: 
	         v_max.set("V.min : -")
	      else: 
		 v_max.set("V.min : "+str(Vmax))
              self.label_vmin.set_text(r"$V_{min}: -       $")
	      self.label_vmax.set_text(r"$V_{max}: -       $")
	      print 'Please waint ...'
	      for obj in self.list_line:
		    obj[4] = True
		    p = vel_indx(obj[3])
		    if self.color_bts[p][1]:
		      obj[0].set_visible(True)
	      for obj in self.list_point:
		    obj[4] = True    
		    p = vel_indx(obj[3])
		    if self.color_bts[p][1]:
		      obj[0].set_visible(True)
	      for obj in self.lline:
		    obj[3] = True    
		    p = vel_indx(obj[2])
		    if self.color_bts[p][1]:
		      obj[0].set_visible(True)            
            
        
        def callback():
            string = e.get()
            try: 
	      s = float(string)
	      if self.browser!= None:
		  B = True
		  for i in range(len(self.browser.id_gal)):
		    if s == self.browser.id_gal[i]:
		      x0, y0 = self.browser.single_info(i)
		      dy = abs(self.xy_lim[3]-self.xy_lim[2])
		      dx = abs(self.xy_lim[0]-self.xy_lim[1])
		      self.xy_lim[0] = x0 + dx/2.
		      self.xy_lim[1] = x0 - dx/2.
		      self.xy_lim[2] = y0 - dy/2.
		      self.xy_lim[3] = y0 + dy/2.
		      self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	              self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	              self.fig.canvas.draw()
	              
	              time.sleep(0.7)
	              for pp in range(3):
			self.browser.group_plot.set_invisible_gal() 
	                self.fig.canvas.draw()
	                time.sleep(0.7)
	                self.browser.group_plot.set_visible_gal(self.browser.id_gal[i]) 
	                self.fig.canvas.draw()
	                time.sleep(0.7)
	              
		      B = False
		      break
		  if B: 
		    print "PGC "+string+" is not available ..."
		    v.set("Not Available")
            except:
	      v.set("Not a number")
	      print "Error: Please enter a valid PGC number"
	      
	     
            
        def callclear():
            v.set("")
        
        def call_reset():
	  self.xy_lim = self.xy_lim_0[:]
	  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	  self.fig.canvas.draw()
        
        def call_zoomout():
	  self.xy_lim = zoom(self.xy_lim, ratio = 10./9)
	  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	  self.fig.canvas.draw()
	  
	
	def call_zoomin():
	  self.xy_lim = zoom(self.xy_lim, ratio = 9/10.)
	  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	  self.fig.canvas.draw()
	  
        
        def call_up():
	  dy = abs(self.xy_lim[3]-self.xy_lim[2])
	  self.xy_lim[3] += 0.15*dy
	  self.xy_lim[2] += 0.15*dy
	  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	  self.fig.canvas.draw()
        
        def call_down():
	  dy = abs(self.xy_lim[3]-self.xy_lim[2])
	  self.xy_lim[3] -= 0.15*dy
	  self.xy_lim[2] -= 0.15*dy
	  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	  self.fig.canvas.draw()

        def call_left():
	  dx = abs(self.xy_lim[0]-self.xy_lim[1])
	  self.xy_lim[0] += 0.15*dx
	  self.xy_lim[1] += 0.15*dx
	  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	  self.fig.canvas.draw()

        def call_right():
	  dx = abs(self.xy_lim[0]-self.xy_lim[1])
	  self.xy_lim[0] -= 0.15*dx
	  self.xy_lim[1] -= 0.15*dx
	  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	  self.fig.canvas.draw()


	  

	  
        label = tk.Label(self.win, text="    ", font=LARGE_FONT)
        #label.pack(pady=10,padx=10)
       
        b = tk.Button(self.win, text="Find", command=callback)
        c = tk.Button(self.win, text="Clear", command=callclear)
        
        b_vel = tk.Button(self.win, text="Set", command=callback_vel)
        c_vel = tk.Button(self.win, text="Clear", command=callclear_vel)  
        d_vel = tk.Button(self.win, text="Show All", command=callshowall_vel)
        
        def callprompt():
	  global event_active
	  event_active = False
	  b.config(state=tk.DISABLED)
	  c.config(state=tk.DISABLED)
	  d.config(state=tk.DISABLED)
	  b_vel.config(state=tk.DISABLED)
	  c_vel.config(state=tk.DISABLED)
	  d_vel.config(state=tk.DISABLED)
	  e_vel.config(state=tk.DISABLED)
	  try:
	    print
	    input_var = input("Enter PGC #: ")
	    v.set(input_var)
	  except:
	    print
	    print "Error: ..."    
	    event_active = True
	  b.config(state=tk.NORMAL)
	  c.config(state=tk.NORMAL)
	  d.config(state=tk.NORMAL)
	  b_vel.config(state=tk.NORMAL)
	  c_vel.config(state=tk.NORMAL)
	  d_vel.config(state=tk.NORMAL)
	  e_vel.config(state=tk.NORMAL)	    
	    
	def callprompt_vel():
	  global event_active
	  event_active = False
	  b.config(state=tk.DISABLED)
	  c.config(state=tk.DISABLED)
	  d.config(state=tk.DISABLED)
	  b_vel.config(state=tk.DISABLED)
	  c_vel.config(state=tk.DISABLED)
	  d_vel.config(state=tk.DISABLED)
	  e_vel.config(state=tk.DISABLED)
	  try:
	    print
	    input_var = input("Enter V-min: ")
	    v_min.set("V.min (km/s): "+str(input_var))
	    input_var = input("Enter V-max: ")
	    v_max.set("V.max (km/s): "+str(input_var))
	  except:
	    print
	    print "Error: ..."
          event_active = True
	  b.config(state=tk.NORMAL)
	  c.config(state=tk.NORMAL)
	  d.config(state=tk.NORMAL)
	  b_vel.config(state=tk.NORMAL)
	  c_vel.config(state=tk.NORMAL)
	  d_vel.config(state=tk.NORMAL)
	  e_vel.config(state=tk.NORMAL)
	  
        d = tk.Button(self.win, text="...", command=callprompt)
        e_vel = tk.Button(self.win, text="...", command=callprompt_vel)
        
        
        b_up = tk.Button(self.win, text="^", command=call_up, font=LARGE_FONT)
        b_down = tk.Button(self.win, text="v", command=call_down, font=LARGE_FONT)
        b_left = tk.Button(self.win, text="<", command=call_left, font=LARGE_FONT)
        b_right = tk.Button(self.win, text=">", command=call_right, font=LARGE_FONT)
        b_zoomin = tk.Button(self.win, text="Zoom in", command=call_zoomin)
        b_zoomout = tk.Button(self.win, text="zoom out", command=call_zoomout)
        b_reset = tk.Button(self.win, text="reset view", command=call_reset)
        
        l1 = tk.Label(self.win, text="PGC #:")
        l2 = tk.Label(self.win, text="             ")
        l3 = tk.Label(self.win, text="        ")

        
        
        e_min.grid(row=0,column=0)
        e_max.grid(row=0,column=1)
        b_vel.grid(row=0,column=2)
        c_vel.grid(row=0,column=3)
        d_vel.grid(row=0,column=4)
        e_vel.grid(row=0,column=5)
        label.grid(row=0,column=6)
        
        l1.grid(row=0,column=7)
        e.grid(row=0,column=8)
        b.grid(row=0,column=9)
        c.grid(row=0,column=10)
        d.grid(row=0,column=11)
        
        l2.grid(row=0,column=12)  # blanck space
        
        b_left.grid(row=0,column=13)
        b_up.grid(row=0,column=14)
        b_down.grid(row=0,column=15)
        b_right.grid(row=0,column=16)
        
        l3.grid(row=0,column=17)  # blanck space
        
        b_zoomin.grid(row=0,column=18)
        b_zoomout.grid(row=0,column=19)
        b_reset.grid(row=0,column=20)








  def show_plot(self):
    
        
     
        filee  = self.filee
        mode  = self.mode
        alpha = self.alpha
        delta = self.delta
        width = self.width
        pgc   = self.pgc
        
    
	mytable = file_deliver(filee)
	id    = mytable['pgc']
	flag  = mytable['flag']
	sgl   = mytable['sgl']
	sgb   = mytable['sgb']
	gl   = mytable['gl']
	gb   = mytable['gb']  
	ra   = mytable['ra']
	dec  = mytable['dec']
	Ks   = mytable['Ks']
	Vls   = mytable['Vls']
	R_2t2 = mytable['R2t_lum']
	nest  = mytable['nest']
	dcf2  = mytable['dcf2']
	ed    = mytable['ed']
	objname = mytable['objname']
	mDist = mytable['mDist']
	mDistErr = mytable['mDistErr']
	sigmaP_lum = mytable['sigmaP_lum']
	
	N_galaxies = len(id)
	ra1    = np.zeros((N_galaxies,), dtype=np.int)
	ra2    = np.zeros((N_galaxies,), dtype=np.int)
	dec1   = np.zeros((N_galaxies,), dtype=np.int)
	dec2   = np.zeros((N_galaxies,), dtype=np.int)
	flag_p = np.zeros((N_galaxies,), dtype=np.int)
	Vls_p  = np.zeros((N_galaxies,), dtype=np.int)
	has_tick = False
	t1     = np.zeros((N_galaxies,), dtype=np.int)
	t2     = np.zeros((N_galaxies,), dtype=np.int)

	if mode == "equatorial":
	  x_cord = ra
	  y_cord = dec
	elif mode == "galactic":
	  x_cord = gl
	  y_cord = gb 
	elif mode == "supergalactic":
	  x_cord = sgl
	  y_cord = sgb
	elif mode == "cartesian":
	  x_cord = sgl
	  y_cord = sgb	

	if pgc != None:
	  print "Using the specified object coordinates: pgc", pgc 
	  index = np.where(nest == pgc)
	  alpha = x_cord[index[0][0]]
	  delta = y_cord[index[0][0]]
	  print alpha, delta
	
	x_min, x_max, y_min, y_max, self.xy_lim = main(alpha, delta, self.fig, self.ax, width=width, mode = mode, V_range=[self.Vmin, self.Vmax], SGX=self.SGX, SGY=self.SGY, SGZ=self.SGZ, Projection=self.Projection)
	self.xy_lim_0 = self.xy_lim[:]


	if mode!='cartesian':
	    if x_max>=360:
	      x_max-= 360; x_min-=360
	    
	    if x_min < 0:
	      ra2[np.where(x_cord>x_min+360)] = 2
	      ra1[np.where(x_cord<x_max)] = 2
	    else:
	      ra2[np.where(x_cord>x_min)] = 1
	      ra1[np.where(x_cord<x_max)] = 1
	    
	    dec1[np.where(y_cord<y_max)] = 1
	    dec2[np.where(y_cord>y_min)] = 1  
	    
        elif mode=='cartesian':
	    sgx, sgy, sgz = xyz_list(sgl, sgb, Vls, dcf2, mDist, flag)
	    
	    
	    
	    if self.Projection == 1:
		ra1[np.where(sgx<self.xy_lim_0[1])] = 1
		ra2[np.where(sgx>self.xy_lim_0[0])] = 1
		dec1[np.where(sgy<self.xy_lim_0[3])] = 1
		dec2[np.where(sgy>self.xy_lim_0[2])] = 1
		if self.zt != None and self.SGZ != None:
		  has_tick = True
		  t1[np.where(sgz<self.SGZ+0.5*abs(self.zt))] = 1
		  t2[np.where(sgz>self.SGZ-0.5*abs(self.zt))] = 1
	    if self.Projection == 2:
		ra1[np.where(sgx<self.xy_lim_0[1])] = 1
		ra2[np.where(sgx>self.xy_lim_0[0])] = 1
		dec1[np.where(sgz<self.xy_lim_0[3])] = 1
		dec2[np.where(sgz>self.xy_lim_0[2])] = 1	
		if self.yt != None and self.SGY != None:
		  has_tick = True
		  t1[np.where(sgy<self.SGY+0.5*abs(self.yt))] = 1
		  t2[np.where(sgy>self.SGY-0.5*abs(self.yt))] = 1
	    if self.Projection == 3:
		ra1[np.where(sgy<self.xy_lim_0[1])] = 1
		ra2[np.where(sgy>self.xy_lim_0[0])] = 1
		dec1[np.where(sgz<self.xy_lim_0[3])] = 1
		dec2[np.where(sgz>self.xy_lim_0[2])] = 1
		if self.xt != None and self.SGX != None:
		  has_tick = True
		  t1[np.where(sgx<self.SGX+0.5*abs(self.xt))] = 1
		  t2[np.where(sgx>self.SGX-0.5*abs(self.xt))] = 1
	    if self.Projection == 10:
		ra1[np.where(sgy<self.xy_lim_0[1])] = 1
		ra2[np.where(sgy>self.xy_lim_0[0])] = 1
		dec1[np.where(sgx<self.xy_lim_0[3])] = 1
		dec2[np.where(sgx>self.xy_lim_0[2])] = 1
		if self.zt != None and self.SGZ != None:
		  has_tick = True
		  t1[np.where(sgz<self.SGZ+0.5*abs(self.zt))] = 1
		  t2[np.where(sgz>self.SGZ-0.5*abs(self.zt))] = 1
	    if self.Projection == 20:
		ra1[np.where(sgz<self.xy_lim_0[1])] = 1
		ra2[np.where(sgz>self.xy_lim_0[0])] = 1
		dec1[np.where(sgx<self.xy_lim_0[3])] = 1
		dec2[np.where(sgx>self.xy_lim_0[2])] = 1
		if self.yt != None and self.SGY != None:
		  has_tick = True
		  t1[np.where(sgy<self.SGY+0.5*abs(self.yt))] = 1
		  t2[np.where(sgy>self.SGY-0.5*abs(self.yt))] = 1
	    if self.Projection == 30:
		ra1[np.where(sgz<self.xy_lim_0[1])] = 1
		ra2[np.where(sgz>self.xy_lim_0[0])] = 1
		dec1[np.where(sgy<self.xy_lim_0[3])] = 1
		dec2[np.where(sgy>self.xy_lim_0[2])] = 1	
		if self.xt != None and self.SGX != None:
		  has_tick = True
		  t1[np.where(sgx<self.SGX+0.5*abs(self.xt))] = 1
		  t2[np.where(sgx>self.SGX-0.5*abs(self.xt))] = 1
		
	#########
	flag_p[np.where(flag<2)] = 1
	#Vls_p[np.where(Vls!=0)] = 1   # 
        cond = ra1 + ra2 + dec1 + dec2 + flag_p +  t1 + t2
        if has_tick:
          indices = np.where(cond==7)
        else:
	  indices = np.where(cond==5)
        #########  	
        
	id_gal      = id[indices]
	flag_gal    = flag[indices]
	ra_gal      = ra[indices]
	dec_gal     = dec[indices]
	sgl_gal     = sgl[indices]
	sgb_gal     = sgb[indices]
	gl_gal      = gl[indices]
	gb_gal      = gb[indices]
	Ks_gal      = Ks[indices]
	Vls_gal     = Vls[indices]
	R_2t2_gal   = R_2t2[indices]
	nest_gal    = nest[indices]
	dcf2_gal    = dcf2[indices]
	ed_gal      = ed[indices]
	objname_gal = objname[indices]
	mDist_gal   = mDist[indices]
	mDistErr_gal   = mDistErr[indices]
	sigmaP_lum_gal = sigmaP_lum[indices]
	
	table_gal = {'pgc':id_gal, 'flag':flag_gal, 'sgl': sgl_gal, 'sgb':sgb_gal, 'gl':gl_gal, 'gb':gb_gal, 'ra':ra_gal, 'dec':dec_gal, 'Ks':Ks_gal, 'Vls':Vls_gal,  \
               'R2t_lum':R_2t2_gal, 'nest':nest_gal, 'dcf2':dcf2_gal, 'ed':ed_gal, 'objname':objname_gal, \
		 'mDist':mDist_gal, 'mDistErr':mDistErr_gal, 'sigmaP_lum':sigmaP_lum_gal}
	
	
	if mode == "equatorial":
	  A, B = transform_list(alpha, delta, ra_gal, dec_gal)
	elif mode == "galactic":
	  A, B = transform_list(alpha, delta, gl_gal, gb_gal)
	elif mode == "supergalactic":
	  A, B = transform_list(alpha, delta, sgl_gal, sgb_gal)
 	elif mode == "cartesian":
	  if self.Projection == 1:
	    A = sgx[indices]
	    B = sgy[indices]
	  if self.Projection == 2:
	    A = sgx[indices]
	    B = sgz[indices]
	  if self.Projection == 3:
	    A = sgy[indices]
	    B = sgz[indices]
	  if self.Projection == 10:
	    A = sgy[indices]
	    B = sgx[indices]
	  if self.Projection == 20:
	    A = sgz[indices]
	    B = sgx[indices]
	  if self.Projection == 30:
	    A = sgz[indices]
	    B = sgy[indices]
	    
	self.line, = plt.plot(A, B, '.', markersize = 1, color='white', markeredgecolor='white', picker=5)   # color='#696969'
	

	for i in range(len(id_gal)):
	  l, = plt.plot([A[i]], [B[i]], 'o', markersize = 2.5, color='white', markeredgecolor=color_table(Vls_gal[i]))
	  self.lline.append([l,color_table(Vls_gal[i]),Vls_gal[i], True, '#696969'])
	#########

	#########
	cond = ra1 + ra2 + dec1 + dec2 +  t1 + t2
        if has_tick:
          indices = np.where(cond==6)
        else:
	  indices = np.where(cond==4)


	id      = id[indices]
	flag    = flag[indices]
	ra      = ra[indices]
	dec     = dec[indices]
	sgl     = sgl[indices]
	sgb     = sgb[indices]
	gl      = gl[indices]
	gb      = gb[indices]
	Ks      = Ks[indices]
	Vls     = Vls[indices]
	R_2t2   = R_2t2[indices]
	nest    = nest[indices]
	dcf2    = dcf2[indices]
	ed      = ed[indices]
	objname = objname[indices]
	mDist   = mDist[indices]
	mDistErr   = mDistErr[indices]
	sigmaP_lum = sigmaP_lum[indices]
	
	table_all = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'dcf2':dcf2, 'ed':ed, 'objname':objname, \
		 'mDist':mDist, 'mDistErr':mDistErr, 'sigmaP_lum':sigmaP_lum}
	#########  
	

	
	if mode=="equatorial":
	  A_cluster, B_cluster = transform_list(alpha, delta, ra, dec)
	elif mode=="galactic":
	  A_cluster, B_cluster = transform_list(alpha, delta, gl, gb)
	elif mode=="supergalactic":
	  A_cluster, B_cluster = transform_list(alpha, delta, sgl, sgb)
	elif mode=="cartesian":
	  if self.Projection == 1:
	    A_cluster = sgx[indices]
	    B_cluster = sgy[indices]
	  if self.Projection == 2:
	    A_cluster = sgx[indices]
	    B_cluster = sgz[indices]
	  if self.Projection == 3:
	    A_cluster = sgy[indices]
	    B_cluster = sgz[indices]
	  if self.Projection == 10:
	    A_cluster = sgy[indices]
	    B_cluster = sgx[indices]
	  if self.Projection == 20:
	    A_cluster = sgz[indices]
	    B_cluster = sgx[indices]
	  if self.Projection == 30:
	    A_cluster = sgz[indices]
	    B_cluster = sgy[indices]
	  
	ClusterPlot(id, A_cluster, B_cluster, Vls, R_2t2, mDist, nest, flag, self, mode=mode)


	self.browser = PointBrowser(A, B, table_gal, A_cluster, B_cluster, table_all, self, mode)

	

	
	self.ax.minorticks_on()
	self.ax.tick_params(which='major', length=7, width=1.5)
	self.ax.tick_params(which='minor', length=4, color='#000033', width=1.0) 
	
	
	if mode!='cartesian':
	   self.ax.set_xlabel(r'$\Delta$' + 'X [deg]', fontsize=14)
	   self.ax.set_ylabel(r'$\Delta$' + 'Y [deg]', fontsize=14)
	elif mode=='cartesian':
	  if self.Projection ==1:
	    self.ax.set_xlabel('SGX [Mpc]', fontsize=14)
	    self.ax.set_ylabel('SGY [Mpc]', fontsize=14)
	  if self.Projection ==2:
	    self.ax.set_xlabel('SGX [Mpc]', fontsize=14)
	    self.ax.set_ylabel('SGZ [Mpc]', fontsize=14)
	  if self.Projection ==3:
	    self.ax.set_xlabel('SGY [Mpc]', fontsize=14)
	    self.ax.set_ylabel('SGZ [Mpc]', fontsize=14)
	  if self.Projection ==10:
	    self.ax.set_xlabel('SGY [Mpc]', fontsize=14)
	    self.ax.set_ylabel('SGX [Mpc]', fontsize=14)
	  if self.Projection ==20:
	    self.ax.set_xlabel('SGZ [Mpc]', fontsize=14)
	    self.ax.set_ylabel('SGX [Mpc]', fontsize=14)
	  if self.Projection ==30:
	    self.ax.set_xlabel('SGZ [Mpc]', fontsize=14)
	    self.ax.set_ylabel('SGY [Mpc]', fontsize=14)	    
	
	
	self.ax.annotate("Catalog Name:", (0.016,0.12), xycoords='figure fraction', size=12, color='maroon')
	my_label = filee
	spl = my_label.split('/')
	if spl[-1] == '':
	  label = spl[-2]
	else:
	  label = spl[-1]
	self.ax.annotate(label, (0.016,0.08), xycoords='figure fraction', size=10, color='black')
	


	  
	def scroll_event(event):
	      #print 'you pressed', event.key, event.button, event.xdata, event.ydata, event.key

	      if event.inaxes == self.ax: 
		if event.key is None and event.button == 'up':
		  self.xy_lim = zoom(self.xy_lim, ratio = 10./9)
		  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
		  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
		  self.fig.canvas.draw()

		elif event.key is None and event.button == 'down':
		  self.xy_lim = zoom(self.xy_lim, ratio = 9/10.)
		  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
		  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
		  self.fig.canvas.draw()

	      
	
	def on_click(event):

	  if event.dblclick and event.button==1:
	    self.browser.unselect(event)
	  if event.dblclick and event.button==3:
	    
	    # close any velocity window
	    if self.inputDialog != None:
	        self.inputDialog.quit()

	    alpha, delta = transform_inv(self.alpha, self.delta, event.xdata, event.ydata)
	    self.controller.alpha = alpha
	    self.controller.delta = delta
	    self.controller.show_frame()
	  elif event.inaxes == self.ax and event.button == 2:
	    self.xy_lim = pan(self.xy_lim, event.xdata, event.ydata)
	    self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	    self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	    self.fig.canvas.draw()
          

	    
	def press_key(event):
	  
	  if event.key == 'up':
	    dy = abs(self.xy_lim[3]-self.xy_lim[2])
	    self.xy_lim[3] += 0.15*dy
	    self.xy_lim[2] += 0.15*dy
	    self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	    self.fig.canvas.draw()
	    
	  elif event.key == 'down':
	    dy = abs(self.xy_lim[3]-self.xy_lim[2])
	    self.xy_lim[3] -= 0.15*dy
	    self.xy_lim[2] -= 0.15*dy
	    self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	    self.fig.canvas.draw()
	    
	  elif event.key == 'right':
	    dx = abs(self.xy_lim[0]-self.xy_lim[1])
	    self.xy_lim[0] -= 0.15*dx
	    self.xy_lim[1] -= 0.15*dx
	    self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	    self.fig.canvas.draw()
	  
	  elif event.key == 'left':
	    dx = abs(self.xy_lim[0]-self.xy_lim[1])
	    self.xy_lim[0] += 0.15*dx
	    self.xy_lim[1] += 0.15*dx
	    self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	    self.fig.canvas.draw()
	  
	  elif event.key == 'ctrl+z':   # zoom-in
	    self.xy_lim = zoom(self.xy_lim, ratio = 9/10.)
	    self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	    self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	    self.fig.canvas.draw()
	    
	  elif event.key == 'ctrl+x':   # zoom-out
            self.xy_lim = zoom(self.xy_lim, ratio = 10./9)
	    self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	    self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	    self.fig.canvas.draw()
	    
	###############################
	for p in range(0,len(col_pallet)):
	  resetax = axes([0.83, 0.9-p*0.05, 0.03, 0.03])
	  self.color_bts.append([Button(resetax, '', color=col_pallet[p], hovercolor=col_pallet[p]),True])
	  
	def cng_color(event):
	  bbox = self.fig.get_window_extent()
	  width, height = bbox.width, bbox.height
	  x, y = event.x, event.y
	  for p in range(0,len(col_pallet)):
	    if x>= 0.83*width and x<=0.86*width and y>=(0.9-p*0.05)*height and y<=(0.93-p*0.05)*height:
	      break 
	  
	  if self.color_bts[p][1]:
	    self.color_bts[p][0].color = 'white'
	    self.color_bts[p][0].label.set_text('off')
	    self.color_bts[p][1] = False
	    for obj in self.list_line:
	      if vel_limits(p, obj[3]):
	         obj[0].set_visible(False)
	    for obj in self.list_point:
	      if vel_limits(p, obj[3]):
	         obj[0].set_visible(False)	
	    for obj in self.lline:
	      if vel_limits(p, obj[2]):
	         obj[0].set_visible(False)	         
	    
	  else:
	    self.color_bts[p][0].color = col_pallet[p]
	    self.color_bts[p][1] = True
	    self.color_bts[p][0].label.set_text('')
	    for obj in self.list_line:
	      if vel_limits(p, obj[3]) and obj[4]:
	         obj[0].set_visible(True)
	    for obj in self.list_point:
	      if vel_limits(p, obj[3]) and obj[4]:
	         obj[0].set_visible(True)	
	    for obj in self.lline:
	      if vel_limits(p, obj[2]) and obj[3]:
	         obj[0].set_visible(True)	
	         
	for p in range(0,len(col_pallet)):
	   self.color_bts[p][0].on_clicked(cng_color)	
	
	
	###############################
	
	axcolor = 'lightgoldenrodyellow'
	resetax = axes([0.03, 0.77, 0.08, 0.03])
	self.resetax_button = Button(resetax, "Reset view", color=axcolor, hovercolor='navajowhite')
	self.resetax_button.label.set_fontsize(10)
	
	
	
	def reset_func(event):
	  self.xy_lim = self.xy_lim_0[:]
	  self.ax.set_xlim(self.xy_lim[0],self.xy_lim[1])
	  self.ax.set_ylim(self.xy_lim[2],self.xy_lim[3])
	  self.fig.canvas.draw()
	  
	  

	self.resetax_button.on_clicked(reset_func)
	###############################
	rax = axes([0.85, 0.1, 0.11, 0.13], axisbg='tan')
        self.radio = RadioButtons(rax, ('No', 'Yes'), active=0)
        annotate('Random Colors ?', (0.86,0.21), xycoords='figure fraction', size=10, color='black')
        
   
        def colorfunc(label):
           if label == 'Yes':
             for obj in self.color_bts:

	       obj[0].ax.set_visible(False)
	     for obj in self.color_labels:
	       obj.set_visible(False)
	       
	     for obj in self.list_line:
	       obj[0].set_color(obj[1])
	     for obj in self.list_point:
	       obj[0].set_color(obj[1]) 
	       obj[0].set_markeredgecolor(obj[1])
	     
	     for obj in self.lline:
	       obj[0].set_markeredgecolor(obj[4])
	       
	       
	   if label == 'No':
	     for obj in self.color_bts:
	       obj[0].ax.set_visible(True)
	     for obj in self.color_labels:
	       obj.set_visible(True)
	       
	     for obj in self.list_line:
	       obj[0].set_color(obj[2])
	     for obj in self.list_point:
	       obj[0].set_color(obj[2]) 
	       obj[0].set_markeredgecolor(obj[2])	       
	     for obj in self.lline:
	       obj[0].set_markeredgecolor(obj[1])
	       
           self.fig.canvas.draw()
	   
      
        self.radio.on_clicked(colorfunc)
        #######################          
        
	
	
	self.fig.canvas.mpl_connect('scroll_event', scroll_event)
	self.fig.canvas.mpl_connect('pick_event', self.browser.onpick)
	self.fig.canvas.mpl_connect('button_press_event', on_click)
	self.fig.canvas.mpl_connect('key_press_event', press_key)
	self.canvas.mpl_connect('button_press_event', lambda event:self.canvas._tkcanvas.focus_set())
	
	#plt.show()
	#plt.savefig('test.ps', dpi=600)
	
	self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)


        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
  
  
  
#################################################################

if __name__ == '__main__':
  
  	if (len(sys.argv) < 2): 
	  print "\nNot enough input arguments ..."
	  print >> sys.stderr, "Use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
	  exit(1)
	  
	opts, args =  arg_parser()
	print "\n------------------------------------"
	print "Input Arguments (provided by User)"
	print "------------------------------------"
	print "Input file:", opts.file
	print "Coordinate system:", opts.coordinate
	print "alpha:", opts.alpha
	print "delta:", opts.delta
	print "width:", opts.width
	print "pgc obj:", opts.pgc
	print "V_min:", opts.Vmin
	print "V_max:", opts.Vmax
	print "SGX: ", opts.SGX
	print "SGY: ", opts.SGY
	print "SGZ: ", opts.SGZ
	print "SGX Tickness: ", opts.XT	
	print "SGY Tickness: ", opts.YT	
	print "SGZ Tickness: ", opts.ZT	
	print "Projection: ", opts.projection
	print "Virgo on: ", opts.G
	virgo_on = opts.G
	print "------------------------------------"
	print "Use You can use \"python "+sys.argv[0]+" -h\""
	print "to see how you can set these values."
	print "------------------------------------"


	if opts.SGX!=None and opts.SGY!=None and opts.SGZ!=None and opts.projection==None:
	  if opts.XT!=None and opts.YT==None and opts.ZT==None: 
	    opts.projection=3
	    print "[Success] Chosen coordinate: SGY-SGZ "
	  if opts.XT==None and opts.YT!=None and opts.ZT==None: 
	    opts.projection=2
	    print "[Success] Chosen coordinate: SGX-SGZ "
	  if opts.XT==None and opts.YT==None and opts.ZT!=None: 
	    opts.projection=1
	    print "[Success] Chosen coordinate: SGX-SGY "
	  
	else: 
	    if opts.projection==None:
	      print "\n[Warning] the projection is not specified correctly ..."
	      print "[Success] Chosen coordinate: SGX-SGY "
	      opts.projection = 1
	      if opts.SGX!=None and opts.SGY!=None:
		opts.projection = 1
	      elif opts.SGX!=None and opts.SGZ!=None:
		opts.projection = 2
		print "[Success] Chosen coordinate: SGX-SGZ "
	      elif opts.SGY!=None and opts.SGZ!=None:
		opts.projection = 3	    
		print "[Success] Chosen coordinate: SGY-SGZ "
	    elif opts.projection % 3 > 2:
	      print "\n[Warning] Wrong projection number ..."
	      print "Still trying to figure out the best projection for your purpose ..."
	      opts.projection = 1
	      if opts.SGX!=None and opts.SGY!=None:
		opts.projection = 1
		print "[Success] Chosen coordinate: SGX-SGY "
	      elif opts.SGX!=None and opts.SGZ!=None:
		opts.projection = 2
		print "[Success] Chosen coordinate: SGX-SGZ "
	      elif opts.SGY!=None and opts.SGZ!=None:
		opts.projection = 3	 
		print "[Success] Chosen coordinate: SGY-SGZ "
	

	
	if opts.coordinate=='equatorial':
	  mode = "equatorial"
	elif opts.coordinate=='galactic':
	  mode = "galactic"
	elif opts.coordinate=='supergalactic':
	  mode = 'supergalactic'
	elif opts.coordinate=='cartesian':
	  mode = 'cartesian'	  
	else:
	  print "\n[Warning] the coordinate system not specified correctly ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
	  print "Using equatorial system by default ...\n"
	  mode = "equatorial"
	
	
	if opts.file != None:
	  filee = opts.file
	else:
	  print "Input filename:", filee
	  print "\n[Error] The input file was not specified correctly ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
	  exit(1)
	
	
	
	if opts.alpha == None:
	  print "\n[Warning] alpha not specified correctly ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
	  print "Using the default value ..."
	  print "alpha = 0\n"
	  alpha = 0
	else:
	  alpha = opts.alpha
	  if alpha >= 360:
	    alpha -= 360*floor(alpha/360)
	  if alpha < 0:
	    alpha *= -1
	    alpha -= 360*floor(alpha/360)
	    alpha = 360. - alpha 



	
	if opts.delta == None:
	  print "\n[Warning] delta not specified correctly ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
	  print "Using the default value ..."
	  print "delta = 0\n"
	  delta = 0
	else:
	  delta = opts.delta 
	  if delta>90 or delta<-90:
	    print "\n[Error] delta must be in [-90,90] range ..."
	    print "You entered: ", delta
	    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
	    if opts.pgc == None:
	      exit(1)
	  
	  
	  
	if opts.width == None:
	  print "\n[Warning] width is not specified correctly ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
	  print "Using the default value ..."
	  print "width = 20 deg/Mpc\n"
	  width = 10
	else:
	  width = opts.width/2.
  

	if opts.Vmin == None:
	  print "\n[Warning] V_min has not been specified ..."
	  print "[Warning] The Lower limit on radial velocity is -Inifinity ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
  

	if opts.Vmax == None:
	  print "\n[Warning] V_max has not been specified ..."
	  print "[Warning] The Upper limit on radial velocity is +Inifinity ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
	
	if opts.Vmin != None and opts.Vmax != None and opts.Vmin>opts.Vmax:
	  print "\n[Error] V_min must be smaller than V_max. Try again ..."
	  print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ...\n"
	  exit(1)
  
  # Main loop of the program
        app = MainWin()
        app.filee  = filee
        app.mode   = mode
        app.alpha  = alpha
        app.delta  = delta
        app.width  = width
        app.pgc    = opts.pgc
        
        app.Vmin   = opts.Vmin
        app.Vmax   = opts.Vmax
        
        app.SGX    = opts.SGX
        app.SGY    = opts.SGY
        app.SGZ    = opts.SGZ
        
        app.Projection    = opts.projection
        
        app.xt    = opts.XT
        app.yt    = opts.YT
        app.zt    = opts.ZT
        
        app.show_frame()
        app.mainloop()
  
  
  
  
  
  