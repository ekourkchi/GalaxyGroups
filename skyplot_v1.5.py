#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2016"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "1.5"
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
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib.widgets import Slider, Button, RadioButtons
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
#################################################
def zoom(i1, i2, j1, j2, ratio=1):

     xc = 0.5*(i1+i2)
     yc = 0.5*(j1+j2)
     delta_x = abs(i2-i1); delta_y =  abs(j2-j1)
     dx = 0.5*ratio*delta_x; dy=0.5*ratio*delta_y
     i1 = xc+dx ; i2=xc-dx
     j1 = yc-dy ; j2=yc+dy   
     return i1,i2,j1,j2
   
#################################################
def pan(i1, i2, j1, j2, xc, yc, ratio=1):

     delta_x = abs(i2-i1); delta_y =  abs(j2-j1)
     dx = 0.5*ratio*delta_x; dy=0.5*ratio*delta_y
     i1 = xc+dx ; i2=xc-dx
     j1 = yc-dy ; j2=yc+dy   
     return i1,i2,j1,j2   
#################################################
class GroupBrowser(object):
  
    def __init__(self):
      
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
	
	d_theta = 0.01
	u = np.arange(0,2*pi,d_theta)
	

	for p in range(len(self.id_gal)):
	   Dist = self.Vls_gal[p]/H0
	   if Dist<1: Dist = 1.
	   if self.mDist_gal[p] != 0: Dist = self.mDist_gal[p]
	   r = 180.*atan(self.R_2t2_gal[p]/Dist)/pi
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
	      r = 180.*atan(R_2t2[i]/Dist)/pi
	      
	      X = u*0.; Y = u*0.
	      
	      RA0 = self.ra[i]
	      DEC0 = self.dec[i]
	      
	      for q in range(0,len(u)):
		x = r*np.cos(u[q]) + RA0
		y = r*np.sin(u[q]) + DEC0
		X[q] = x
		Y[q] = y
		
	      if r <= 10000:
		  line, = plt.plot(X,Y, ':', markersize = 1, color='red', picker=False, visible=False)

	      i+=1  
	      X = []
	      Y = []
	      while i<N and flag[i]==1: 
		RA0 = self.ra[i]
		DEC0 = self.dec[i]
		X.append(RA0)
		Y.append(DEC0)
		i+=1
	      group_plt = plt.plot(X, Y, 'o', markersize = 3, mfc='none', markeredgecolor = 'red', visible=False)
	      
	      
	      if line != None and group_plt!= None:
		#line.set_visible(True)
		#for galplot in group_plt:
		   #galplot.set_visible(True)
	        plot_list = [nest_center, line, group_plt]
	        self.listplot.append(plot_list)
	      
	      
	      if i<N and flag[i]==2: 
		gr = 2
	      else:
		break
    
    
    def set_visible_gal(self, pgc_no):  # that's the nest pgc
      if self.listplot_gal != None:
	for plot_list in self.listplot_gal:
	  if plot_list[0] == pgc_no:
	    plot_list[1].set_visible(True)
	    break


    def set_invisible_gal(self):  # that's the nest pgc
      if self.listplot_gal != None:
	for plot_list in self.listplot_gal:
	    plot_list[1].set_visible(False)



    def set_visible(self, pgc_no):  # that's the nest pgc
      
      if self.listplot != None:
	for plot_list in self.listplot:
	  if plot_list[0] == pgc_no:
	    
	    plot_list[1].set_visible(True)
	    for galplot in plot_list[2]:
	      galplot.set_visible(True)
	    
	    break



    def set_invisible(self):  # that's the nest pgc
      
      if self.listplot != None:
	for plot_list in self.listplot:
	    plot_list[1].set_visible(False)
	    for galplot in plot_list[2]:
	      galplot.set_visible(False)

      
#################################################

class PointBrowser(object):
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """

    def __init__(self):
        self.lastind = 0

        self.selected, = ax.plot(A, B, '*', ms=6, alpha=0.8,
                                 color='red', markeredgecolor='red', visible=False)
	
	self.group_plot = GroupBrowser()

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
      
      if event.mouseevent.button == 1 :   # It's just active at left mouse click
        if event.artist != line:
            return True

        N = len(event.ind)
        if not N:
            return True

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        distances = np.hypot(x - A[event.ind], y - B[event.ind])
        indmin = distances.argmin()
        dataind = event.ind[indmin]
        
        if N>1: 
	  ind = event.ind
	  print "\n\n*****************************************************"
	  print "***** All Neighbour Galaxies ************************"
	  print "*****************************************************"

	  myTable = Table()
	  myTable.add_column(Column(data=id_gal[ind], name='pgc'))
	  myTable.add_column(Column(data=flag_gal[ind], name='flag'))
	  myTable.add_column(Column(data=ra_gal[ind],  name='ra'))  
	  myTable.add_column(Column(data=dec_gal[ind],  name='dec')) 
	  #myTable.add_column(Column(data=gl_gal[ind],  name='gl')) 
	  #myTable.add_column(Column(data=gb_gal[ind],  name='gb'))    
	  myTable.add_column(Column(data=sgl_gal[ind],  name='sgl')) 
	  myTable.add_column(Column(data=sgb_gal[ind],  name='sgb')) 
	  myTable.add_column(Column(data=Ks_gal[ind],  name='Ks'))
	  myTable.add_column(Column(data=Vls_gal[ind],  name='Vls')) 
	  myTable.add_column(Column(data=dcf2_gal[ind],  name='dcf2')) 
	  myTable.add_column(Column(data=ed_gal[ind],  name='ed')) 
	  myTable.add_column(Column(data=mDist_gal[ind],  name='mDist')) 
	  myTable.add_column(Column(data=mDistErr_gal[ind],  name='mDistErr')) 
	  myTable.add_column(Column(data=nest_gal[ind],  name='nest')) 
	  myTable.add_column(Column(data=objname_gal[ind],  name='objname')) 
	  print myTable
        
        ind = [dataind]
        print "\n\n------ Selected Galaxy -------------------------------"

        myTable = Table()
        myTable.add_column(Column(data=id_gal[ind], name='pgc'))
        myTable.add_column(Column(data=flag_gal[ind], name='flag'))
        myTable.add_column(Column(data=ra_gal[ind],  name='ra'))  
        myTable.add_column(Column(data=dec_gal[ind],  name='dec')) 
        #myTable.add_column(Column(data=gl_gal[ind],  name='gl')) 
        #myTable.add_column(Column(data=gb_gal[ind],  name='gb'))    
        myTable.add_column(Column(data=sgl_gal[ind],  name='sgl')) 
        myTable.add_column(Column(data=sgb_gal[ind],  name='sgb'))
        myTable.add_column(Column(data=Ks_gal[ind],  name='Ks'))
        myTable.add_column(Column(data=Vls_gal[ind],  name='Vls')) 
        myTable.add_column(Column(data=dcf2_gal[ind],  name='dcf2')) 
        myTable.add_column(Column(data=ed_gal[ind],  name='ed')) 
        myTable.add_column(Column(data=mDist_gal[ind],  name='mDist')) 
        myTable.add_column(Column(data=mDistErr_gal[ind],  name='mDistErr')) 
        myTable.add_column(Column(data=nest_gal[ind],  name='nest')) 
        myTable.add_column(Column(data=objname_gal[ind],  name='objname')) 
        print myTable        
        
        
        

        self.lastind = dataind
        self.update()

    def update(self):
        if self.lastind is None:
            return

        dataind = self.lastind
        
        self.group_plot.set_invisible()
        self.group_plot.set_invisible_gal()
        self.group_plot.set_visible(nest_gal[dataind])
        self.group_plot.set_visible_gal(id_gal[dataind])
        # set the desired point visible
        self.selected.set_visible(True)
        self.selected.set_data(A[dataind], B[dataind]) 

        fig.canvas.draw()

    def unselect(self, event):
        

        if self.lastind is None:
            return
	  
        if event.dblclick:
          self.group_plot.set_invisible()
          self.group_plot.set_invisible_gal()
	  dataind = self.lastind
	  # set the desired point invisible
	  self.selected.set_visible(False)
	  self.selected.set_data(A[dataind], B[dataind]) 
	  fig.canvas.draw()
       
#################################################
def table_deliver(file):
  
  

  try:
    mytable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
   
  except:
    print "\n[Error] The catalog \""+file+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)
  
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
  
  print "\n[Success] The catalog \""+file+"\" was sucessfully loaded ..."
  print "and has "+str(len(id))+" entries ...\n"
  
  
  out_table = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'dcf2':dcf2, 'ed':ed, 'objname':objname, \
		'mDist':mDist, 'mDistErr':mDistErr}
  return out_table
#################################################
def file_deliver(file):
  
  prefix = file.split('.')[0] 
  
  if prefix == 'north' or prefix == 'south':
    return table_deliver(file)
  
  Have_north = os.path.isfile('north.'+file)
  Have_south = os.path.isfile('south.'+file)
  
  if Have_north and not Have_south:
    print "\n[Warning] The catalog \""+'south.'+file+"\" is not available, ..."
    return table_deliver('north.'+file)
  
  elif Have_south and not Have_north:
    print "\n[Warning] The catalog \""+'north.'+file+"\" is not available, ..."
    return table_deliver('south.'+file)
  
  elif not Have_north and not Have_south:

    print "\n[Warning] The catalog \""+'south.'+file+"\" is not available, ..."
    print "[Warning] The catalog \""+'north.'+file+"\" is not available, ...\n"
    
    return table_deliver(file)
   
  
  else:
    north_table = table_deliver('north.'+file)
    south_table = table_deliver('south.'+file)
    
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


    out_table = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'dcf2':dcf2, 'ed':ed, 'objname':objname, \
		 'mDist':mDist, 'mDistErr':mDistErr}    
    
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
 - Version: v1.50
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
                      help="""Coordinate system (i.e. equatorial, galactic or supergalactic)""")
    parser.add_option("-p", "--pgc",
                      type='float', action='store',
                      help="""pgc number (optional, specify the plot center)""")    
    parser.add_option("-w", "--width",
                      type='float', action='store',
                      help="""The width of the plot in degrees""")    
    #parser.add_option('-n', "--north",
                  #action="store_true", dest="north", default=False,
                  #help="Selecting the north galactic sky")
      

    
    
    
    
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

def main(alpha, delta, fig, ax, width=10, mode='equatorial'):
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
  
  
  ax.annotate(x0_text+str(alpha), (0.016,0.86), xycoords='figure fraction', size=11)
  ax.annotate(y0_text+str(delta), (0.016,0.82), xycoords='figure fraction', size=11)
 
  ax.annotate("Coordinates ", (0.016,0.73), xycoords='figure fraction', size=14, color='maroon')
   
  Ux = ax.annotate(" ", (0.016,0.69), xycoords='figure fraction', size=12)
  Uy = ax.annotate(" ", (0.016,0.65), xycoords='figure fraction', size=12)
  
  ax.annotate("WCS ...", (0.016,0.60), xycoords='figure fraction', size=14, color='maroon')
   
  Ura = ax.annotate(" ", (0.016,0.56), xycoords='figure fraction', size=11)
  Ualf = ax.annotate(" ", (0.016,0.54), xycoords='figure fraction', size=10)
  Udec = ax.annotate(" ", (0.016,0.50), xycoords='figure fraction', size=11)
  Udelt = ax.annotate(" ", (0.016,0.47), xycoords='figure fraction', size=10)
  
  ax.annotate("Galactic", (0.016,0.42), xycoords='figure fraction', size=14, color='maroon')
  Ugl = ax.annotate(" ", (0.016,0.38), xycoords='figure fraction', size=11)
  Ugb = ax.annotate(" ", (0.016,0.34), xycoords='figure fraction', size=11)
  
  ax.annotate("Super-Galactic", (0.016,0.29), xycoords='figure fraction', size=14, color='maroon')
  Usgl = ax.annotate(" ", (0.016,0.25), xycoords='figure fraction', size=11)
  Usgb = ax.annotate(" ", (0.016,0.21), xycoords='figure fraction', size=11)
  
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

 
  # grid-Vertical
  for ra in ra_list:
     A=[]; B=[]
     for dec in np.arange(dec_min, dec_max+step, step):
        a, b = transform(alpha, delta, ra, dec)
        A.append(a)
        B.append(b)
     ax.plot(A, B, ':', color='#696969')

  

  # grid-Horiontal
  for dec in np.arange(dec_min, dec_max+5, 5):
     A=[]; B=[]
     for ra in np.arange(ra_min, ra_list_max+step, step):
        a, b = transform(alpha, delta, ra, dec)
        A.append(a)
        B.append(b)
     ax.plot(A, B, ':', color='#696969')  

  
  plt.xlim(width,-1*width)
  plt.ylim(-1*width,width)
  i1 = width
  i2 = -1*width
  j1 = -1*width
  j2 = width

  xvz = tangent(alpha, delta, xyz(alpha, delta))
  a, b = spherical(xvz)
  plt.plot([a], [b], '*', markersize = 10, color= 'black')
  
  
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
       
        Ux.set_text("X: "+'{:.2f}'.format(a))
        Uy.set_text("Y: "+'{:.2f}'.format(b))
        
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

        
        draw()
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
        draw()
   
   
  fig.canvas.mpl_connect('motion_notify_event', in_motion)
  
  return ra_list_min, ra_list_max, dec_min, dec_max
  
#################################################################
def ClusterPlot(id, ra, dec, Vls, R_2t2, mDist, nest, flag):
  
  N = len(id)
  NoGroups = len(id[np.where(flag==2)])
  print "Number of groups: ", NoGroups
  
  if NoGroups == 0:
    print "[Warning] No group found in data-base ..." 
    print "Check the input catalog and choose the right option ...\n" 
    return

  i = 0 
  if NoGroups!=0:
    while flag[i] != 2:
      i+=1
  
  gr = flag[i]
  while gr == 2:
    
    random.seed(nest[i])
    Red, Green, Blue = random.random(), random.random(), random.random()
    Dist = Vls[i]/H0
    if Dist<1: Dist = 1.
    
    #if mDist[i] != 0: Dist = mDist[i]
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
        line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Green, Blue), picker=False) 
	line.set_dashes([8, 3]) 

    i+=1  
    X = []
    Y = []
    while i<N and flag[i]==1: 
      RA0 = ra[i]
      DEC0 = dec[i]
      X.append(RA0)
      Y.append(DEC0)
      i+=1
    plt.plot(X, Y, 'o', markersize = 3, color=(Red, Green, Blue), markeredgecolor = (Red, Green, Blue))
    
        
    
    
    if i<N and flag[i]==2: 
      gr = 2
    else:
       break
  
  
  
  
#################################################################

if __name__ == '__main__':
  
  i1=0; i2=0;
  j1=0; j2=0;
  
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
  print "------------------------------------"
  
  
  if opts.coordinate=='equatorial':
    mode = "equatorial"
  elif opts.coordinate=='galactic':
    mode = "galactic"
  elif opts.coordinate=='supergalactic':
    mode = 'supergalactic'
  else:
    print "\n[Warning] the coordinate system not specified correctly ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ..."
    print "Using equatorial system by default ...\n"
    mode = "equatorial"
  
  
  if opts.file != None:
    file = opts.file
  else:
    print "\n[Error] The input file was not specified correctly ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)
  
  
   
  if opts.alpha == None:
    print "\n[Warning] alpha not specified correctly ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ..."
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
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ..."
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
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ..."
    print "Using the default value ..."
    print "width = 20 deg\n"
    width = 10
  else:
    width = opts.width/2.

  

  fig = plt.figure(figsize=(10,8 ), dpi=100)
  
  resetax = axes([0.01, 0.18, 0.16, 0.77])
  info_bts = Button(resetax, '', color='lightgray', hovercolor='lightgray')
  
  resetax = axes([0.01, 0.05, 0.16, 0.1])
  info_bts = Button(resetax, '', color='lightgray', hovercolor='lightgray')
  
  ax = fig.add_axes([0.25, 0.07, 0.72,  0.9]) 
  #subplots_adjust(left=0.25, bottom=0.25)
  fig.patch.set_facecolor('lightgray')

  
  
  
  mytable = file_deliver(file)
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
  
  N_galaxies = len(id)
  ra1    = np.zeros((N_galaxies,), dtype=np.int)
  ra2    = np.zeros((N_galaxies,), dtype=np.int)
  dec1   = np.zeros((N_galaxies,), dtype=np.int)
  dec2   = np.zeros((N_galaxies,), dtype=np.int)
  flag_p = np.zeros((N_galaxies,), dtype=np.int)
  
  

  
  
  if mode == "equatorial":
     x_cord = ra
     y_cord = dec
  elif mode == "galactic":
     x_cord = gl
     y_cord = gb 
  elif mode == "supergalactic":
     x_cord = sgl
     y_cord = sgb
  

  if opts.pgc != None:
    print "Using the specified object coordinates: pgc", opts.pgc 
    index = np.where(nest == opts.pgc)
    alpha = x_cord[index[0][0]]
    delta = y_cord[index[0][0]]
    print alpha, delta
  x_min, x_max, y_min, y_max = main(alpha, delta, fig, ax, width=width, mode = mode)


  
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
  flag_p[np.where(flag<2)] = 1
  
  
  #########
  cond = ra1 + ra2 + dec1 + dec2 + flag_p
  indices = np.where(cond==5)
  
  ra_gal      = ra[indices]
  dec_gal     = dec[indices]
  id_gal      = id[indices]
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
  flag_gal    = flag[indices]
  mDist_gal   = mDist[indices]
  mDistErr_gal   = mDistErr[indices]
  
  if mode == "equatorial":
     A, B = transform_list(alpha, delta, ra_gal, dec_gal)
  elif mode == "galactic":
     A, B = transform_list(alpha, delta, gl_gal, gb_gal)
  elif mode == "supergalactic":
     A, B = transform_list(alpha, delta, sgl_gal, sgb_gal)
 
  
  line, = ax.plot(A, B, 'o', markersize = 2.5, color='white', markeredgecolor='#696969', picker=5)   # color='#696969'
  #########

  #########
  cond = ra1 + ra2 + dec1 + dec2
  indices = np.where(cond==4)


  id    = id[indices]
  flag  = flag[indices]
  sgl   = sgl[indices]
  sgb   = sgb[indices]
  gl    = gl[indices]
  gb    = gb[indices]
  ra    = ra[indices]
  dec   = dec[indices]
  Vls   = Vls[indices]
  R_2t2 = R_2t2[indices]
  nest  = nest[indices]
  mDist = mDist[indices]
  #########  
  

  
  if mode=="equatorial":
     A_cluster, B_cluster = transform_list(alpha, delta, ra, dec)
  elif mode=="galactic":
     A_cluster, B_cluster = transform_list(alpha, delta, gl, gb)
  elif mode=="supergalactic":
     A_cluster, B_cluster = transform_list(alpha, delta, sgl, sgb)
  
  ClusterPlot(id, A_cluster, B_cluster, Vls, R_2t2, mDist, nest, flag)



  browser = PointBrowser()

  

  
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0) 
  
  
  plt.xlabel(r'$\Delta$' + 'X [deg]', fontsize=14)
  plt.ylabel(r'$\Delta$' + 'Y [deg]', fontsize=14)
  
  
  
  
  ax.annotate("Catalog Name:", (0.016,0.12), xycoords='figure fraction', size=12, color='maroon')
  ax.annotate(file, (0.016,0.08), xycoords='figure fraction', size=10, color='black')
  
    
  def scroll_event(event):
         #print 'you pressed', event.key, event.button, event.xdata, event.ydata, event.key
            
         global i1, i2, j1, j2
         
         
         if event.inaxes == ax: 
	   if event.key is None and event.button == 'up':
	     i1,i2,j1,j2 = zoom(i1, i2, j1, j2, ratio = 10./9)
	     ax.set_xlim(i1,i2)
             ax.set_ylim(j1,j2)
             draw()
           elif event.key is None and event.button == 'down':
	     i1,i2,j1,j2 = zoom(i1, i2, j1, j2, ratio = 9/10.)
	     ax.set_xlim(i1,i2)
             ax.set_ylim(j1,j2)
             draw()
         
  
  def on_click(event):
    global i1, i2, j1, j2
    if event.dblclick and event.button==1:
      browser.unselect(event)
    elif event.inaxes == ax and event.button == 2:
      i1, i2, j1, j2 = pan(i1, i2, j1, j2, event.xdata, event.ydata)
      ax.set_xlim(i1,i2)
      ax.set_ylim(j1,j2)
      draw()
      
  
      
    
  
  fig.canvas.mpl_connect('scroll_event', scroll_event)
  fig.canvas.mpl_connect('pick_event', browser.onpick)
  fig.canvas.mpl_connect('button_press_event', on_click)
  
  
  plt.show()
  #plt.savefig('test.ps', dpi=600)
  
  
  
  
  
  
  
  
  
  