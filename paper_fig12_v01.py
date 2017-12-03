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

event_active = True

virgo_on = False

LARGE_FONT= ("Verdana", 12)

col_pallet = ['darkviolet', 'blue', 'deepskyblue', 'forestgreen', 'y', 'gold', 'darkorange', 'red', 'magenta', 'maroon', 'sienna', 'slategrey', 'black']
vel_pallet = [-100000, 0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 100000]
#################################################

def ClusterPlot(id, ra, dec, Vls, R_2t2, mDist, nest, flag, alpha, delta, mode='none'):
  
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
    line, = plt.plot(X,Y, markersize = 2, color='black', picker=False) 
    #line.set_dashes([8, 3]) 
    
    
    if mode == "supergalactic":
      xA, yB = transform_list(alpha, delta, [110,105,105,96,0], [-3,-3,-4,-4,0])
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
        line, = plt.plot(X,Y, markersize = 2, color=colour, picker=False) 
	#line.set_dashes([8, 3]) 
	
	if nest[i] in [13505, 13434, 12412, 14084, 13620, 12626]:
            vertices = []
            for j in range(len(X)):
                vertices.append([X[j],Y[j]])
            vertices.append([X[0],Y[0]])
            ax.add_patch(Polygon(vertices, closed=True, fill=True, color='yellow'))
            
            
            
	

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
    
    if i<N and flag[i]==2: 
      gr = 2
    else:
       break
  
  

#################################################
#################################################################
################################################################  
#################################################################
def color_table(Vls):
  
  if Vls>=800 and Vls < 1000: p=1
  elif Vls>=1000 and Vls < 1200: p=3
  elif Vls>=1200 and Vls < 1400: p=6
  elif Vls>=1400 and Vls < 1600: p=7
  elif Vls>=1600 and Vls < 1800: p=7
  elif Vls>=1800 and Vls < 2000: p=8
  elif Vls>=2000 and Vls < 2200: p=10
  else: p=12
  
  
  return col_pallet[p]
  


 
#################################################################
################################################################   
def setcolor(x, color='black'):
     for m in x:
         for t in x[m][1]:
             t.set_color(color)
             s = t.get_text()
             n = len(s)
     
             if s[-1] != 'W' and s[-1] != 'E' and  s[-1] != 'N' and s[-1] != 'S':
               number = int(s[0:n-1])
             else:
               number = int(s[0:n-2])    
               
             if s[-1] == 'W' or s[-1] == 'S':
               number *= -1
               
             if s[-1] == 'W':
	       number += 360
       
             s = str(number)
             t.set_text(u''+s+u'\xb0')

########################################################################## 
########################################################################## 
########################################################################## 
########################################################################## 

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
#################################################################

def main(alpha, delta, fig, ax, width=10, mode='equatorial', V_range=None, SGX=None, SGY=None, SGZ=None, Projection=None):
  global i1, i2, j1, j2
  
  R = 2*width
  step = 0.5

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
      #plt.plot([a], [b], '*', markersize = 10, color= 'black')
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
      #plt.plot([0.5*(i1+i2)], 0.5*(j1+j2), '*', markersize = 10, color= 'black')   # central star 
  
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
  

  xy_lim = [i1, i2, j1, j2]
  
  return ra_list_min, ra_list_max, dec_min, dec_max, xy_lim 
  
#################################################################
def table_deliver_super(filee):
  
  

  try:
    mytable = np.genfromtxt(filee , delimiter='|', filling_values="-100000", names=True, dtype=None )
   
  except:
    print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot2 -h') \n"
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

  r1t_lum    = mytable['r1t_lum']
  dist       = mytable['dist']

  mDist      = mytable['mDist']
  mDistErr   = mytable['mDistErr']
  sigmaP_lum = mytable['sigmaP_lum']
  
  print "\n[Success] The catalog \""+filee+"\" was sucessfully loaded ..."
  print "and has "+str(len(id))+" entries ...\n"
  
  
  out_table = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'r1t_lum':r1t_lum, 'dist':dist, \
		'mDist':mDist, 'mDistErr':mDistErr, 'sigmaP_lum':sigmaP_lum}
  return out_table
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
#################################################
def file_deliver_super(filee):
  
  prefix = filee.split('.')[0] 
  
  if prefix == 'north' or prefix == 'south':
    return table_deliver_super(filee)
  
  Have_north = os.path.isfile('north.'+filee)
  Have_south = os.path.isfile('south.'+filee)
  
  if Have_north and not Have_south:
    print "\n[Warning] The catalog \""+'south.'+filee+"\" is not available, ..."
    return table_deliver_super('north.'+filee)
  
  elif Have_south and not Have_north:
    print "\n[Warning] The catalog \""+'north.'+filee+"\" is not available, ..."
    return table_deliver_super('south.'+filee)
  
  elif not Have_north and not Have_south:

    print "\n[Warning] The catalog \""+'south.'+filee+"\" is not available, ..."
    print "[Warning] The catalog \""+'north.'+filee+"\" is not available, ...\n"
    
    return table_deliver_super(filee)
   
  
  else:
    north_table = table_deliver_super('north.'+filee)
    south_table = table_deliver_super('south.'+filee)
    
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

    r1t_lum = np.concatenate((north_table['r1t_lum'][north_ind],south_table['r1t_lum'][south_ind]))
    dist = np.concatenate((north_table['dist'][north_ind],south_table['dist'][south_ind]))
    
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

    r1t_lum1 = np.concatenate((north_table['r1t_lum'][north_ind],south_table['r1t_lum'][south_ind]))
    dist1 = np.concatenate((north_table['dist'][north_ind],south_table['dist'][south_ind]))
    
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

    r1t_lum = np.concatenate((r1t_lum, r1t_lum1)) 
    dist = np.concatenate((dist, dist1)) 

    mDist = np.concatenate((mDist, mDist1)) 
    mDistErr = np.concatenate((mDistErr, mDistErr1))
    sigmaP_lum = np.concatenate((sigmaP_lum, sigmaP_lum1))

    out_table = {'pgc':id, 'flag':flag, 'sgl': sgl, 'sgb':sgb, 'gl':gl, 'gb':gb, 'ra':ra, 'dec':dec, 'Ks':Ks, 'Vls':Vls,  \
               'R2t_lum':R_2t2, 'nest':nest, 'r1t_lum':r1t_lum,  'dist':dist, \
		 'mDist':mDist, 'mDistErr':mDistErr, 'sigmaP_lum':sigmaP_lum}    
    
    return out_table


#################################################
def extractPGC(id, grp=False, supergrp=False):
  
  if not grp and not supergrp:
    return id
  
  
  if grp:
    pgc = int(id)%100000000
  
  
  if supergrp:
    grp = int(id)%10000000000
    pgc = int(grp)%100000000
  
  return pgc  
########################################################################## 
########################################################################## 
########################################################################## 
  
if __name__ == '__main__':
  
        filee = 'all.iter.2.v44.group'
        alpha = 261.660
        delta = -41.67
        width = 25
        mode = 'supergalactic'
        Vmin = 800
        Vmax = 2200
        Projection = None
        SGX = None
        SGY = None
        SGZ = None
        
        fig = plt.figure(figsize=(8,9.5), dpi=100)
	        
        ax = fig.add_axes([0.12, 0.23, 0.85,  0.85*8./9.5])
        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.set_xlabel('SGL [deg]', fontsize=18, labelpad=25)
        ax.set_ylabel('SGB [deg]', fontsize=18, labelpad=40)
           
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
	Vls_min  = np.zeros((N_galaxies,), dtype=np.int)
	Vls_max  = np.zeros((N_galaxies,), dtype=np.int)
	has_tick = False
	t1     = np.zeros((N_galaxies,), dtype=np.int)
	t2     = np.zeros((N_galaxies,), dtype=np.int)

        
        
        ## mode = supergalactic
	x_cord = sgl
	y_cord = sgb
	
	x_min, x_max, y_min, y_max, xy_lim = main(alpha, delta, fig, ax, width=width, mode = mode, V_range=[Vmin, Vmax], SGX=SGX, SGY=SGY, SGZ=SGZ, Projection=Projection)
	xy_lim_0 = xy_lim[:]
	
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


	#########
	flag_p[np.where(flag<2)] = 1
	Vls_min[np.where(Vls>=800)] = 1   # 
	Vls_max[np.where(Vls<=2200)] = 1   # 
        cond = ra1 + ra2 + dec1 + dec2 + flag_p +  t1 + t2 + Vls_min + Vls_max

        indices = np.where(cond==7)
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
	



        A, B = transform_list(alpha, delta, sgl_gal, sgb_gal)
        
	for i in range(len(id_gal)):
	  l, = plt.plot([A[i]], [B[i]], 'o', markersize = 2.5, color='white', markeredgecolor=color_table(Vls_gal[i]))



        cond = ra1 + ra2 + dec1 + dec2 + t1 + t2 + Vls_min + Vls_max
        indices = np.where(cond==6)


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
	
	
	A_cluster, B_cluster = transform_list(alpha, delta, sgl, sgb)
	
	
	ClusterPlot(id, A_cluster, B_cluster, Vls, R_2t2, mDist, nest, flag, alpha, delta, mode=mode)
	
         
         
        
        ########## Supergroups
        filee = 'all.iter.2.v44.supergroup'
        mytable = file_deliver_super(filee)
        
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

        r1t_lum    = mytable['r1t_lum']
        dist    = mytable['dist']
 
	mDist = mytable['mDist']
	mDistErr = mytable['mDistErr']
	sigmaP_lum = mytable['sigmaP_lum']
	
	superGroups_ID = [13418,13505]
	
	for i in range(len(id)):
            
            if extractPGC(id[i], supergrp=True) in superGroups_ID and flag[i] == 5:
                RA0, DEC0 = transform(alpha, delta, sgl[i], sgb[i])
                Dist = dist[i]
                Vls_grp = Vls[i]
                r = 180.*atan(r1t_lum[i]/Dist)/pi   # mode!='cartesian'
                
                d_theta = 0.01
                u = np.arange(0,2*pi,d_theta)
                X = u*0.
                Y = u*0.
                
                colour = color_table(Vls_grp)
            
                for q in range(0,len(u)):
                   x = r*np.cos(u[q]) + RA0
                   y = r*np.sin(u[q]) + DEC0
                   X[q] = x
                   Y[q] = y            
	
                if r <= 10000:
                   line, = plt.plot(X,Y, '-', markersize = 2, color=color_table(Vls_grp))
                   line.set_dashes([8, 3]) 
	

        
#####################################################
## DO NOT change the order
        ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
        ax2.set_axis_off()
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.xaxis.set_ticks_position('none')
        ax2.yaxis.set_ticks_position('none')        

	ax2.annotate(r'$ \/\/\/ 800 \/ \leq \/ V_{LS} \/ < \/ 1000$', (0.1-0.03,0.135), xycoords='figure fraction', size=15, color='black')     # p = 1 
	ax2.annotate(r'$\/ 1000 \/ \leq \/ V_{LS} \/ < \/ 1200$', (0.1-0.03,0.135-0.05), xycoords='figure fraction', size=15, color='black')        # p = 2 
	ax2.annotate(r'$\/ 1200 \/ \leq \/ V_{LS} \/ < \/ 1400$', (0.1-0.03,0.135-2*0.05), xycoords='figure fraction', size=15, color='black')        # p = 3 
	
	ax2.annotate(r'$1400 \/ \leq \/ V_{LS} \/ < \/ 1600$', (0.4-0.03,0.135), xycoords='figure fraction', size=15, color='black')         # p = 6 
	ax2.annotate(r'$1600 \/ \leq \/ V_{LS} \/ < \/ 1800$', (0.4-0.03,0.135-0.05), xycoords='figure fraction', size=15, color='black')         # p = 7 
	ax2.annotate(r'$1800 \/ \leq \/ V_{LS} \/ < \/ 2000$', (0.4-0.03,0.135-2*0.05), xycoords='figure fraction', size=15, color='black')         # p = 8 
	
	ax2.annotate(r'$2000 \/ \leq \/ V_{LS} \/ < \/ 2200\/\/ km\/ s^{-1}  $', (0.7-0.03,0.135), xycoords='figure fraction', size=15, color='black')        # p = 10 


	p_lst=[1,2,3,6,7,8,10]     
	p = 0
        for m in [0.1,0.4]:
          for n in [0.135,0.135-0.05,0.135-0.1]:
              ax2.add_patch(patches.Rectangle((m-0.055-0.03, n-0.005), 0.035, 0.025, color=col_pallet[p_lst[p]]))
              p+=1
        ax2.add_patch(patches.Rectangle((0.7-0.055-0.03, 0.135-0.005), 0.035, 0.025, color=col_pallet[p_lst[p]]))
        
        
        
        
        ax2.annotate(r'$240^o$', (0.92,0.95), xycoords='figure fraction', size=16, color='black')    
        ax2.annotate(r'$250^o$', (0.735,0.95), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$260^o$', (0.556,0.95), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$270^o$', (0.377,0.95), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$280^o$', (0.193,0.95), xycoords='figure fraction', size=16, color='black') 
        
        ax2.annotate(r'$300^o$', (0.209,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$290^o$', (0.302,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$280^o$', (0.386,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$270^o$', (0.460,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$260^o$', (0.530,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$250^o$', (0.606,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$240^o$', (0.688,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$230^o$', (0.769,0.208), xycoords='figure fraction', size=16, color='black') 
        ax2.annotate(r'$220^o$', (0.865,0.208), xycoords='figure fraction', size=16, color='black') 
        
        ax2.annotate(r'$-60^o$', (0.132,0.253), xycoords='figure fraction', size=16, color='blue', rotation=30) 
        ax2.annotate(r'$-50^o$', (0.059,0.379), xycoords='figure fraction', size=16, color='blue') 
        ax2.annotate(r'$-40^o$', (0.059,0.536), xycoords='figure fraction', size=16, color='blue') 
        ax2.annotate(r'$-30^o$', (0.059,0.687), xycoords='figure fraction', size=16, color='blue') 
        ax2.annotate(r'$-20^o$', (0.059,0.836), xycoords='figure fraction', size=16, color='blue') 
        
        #plt.savefig('fig12_v01.eps', dpi=300)  # the best view and size happens in eps 
	plt.show()
#####################################################






