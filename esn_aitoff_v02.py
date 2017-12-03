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
col_pallet = ['darkviolet', 'blue', 'deepskyblue', 'forestgreen', 'y', 'gold', 'darkorange', 'red', 'magenta', 'maroon', 'sienna', 'slategrey', 'black']
vel_pallet = [-100000, 0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 100000]

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
def file_deliver(filee):
  
  try:
    mytable = np.genfromtxt(filee , delimiter='|', filling_values="-100000", names=True, dtype=None )
   
  except:
    print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
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
  
  fig = plt.figure(figsize=(12, 8), dpi=100)
  ax = fig.add_subplot(111, projection="aitoff")
  plt.title("Supergalactic Aitoff Projection", y=1.08)
  ax.grid(True)
  ax.set_xticklabels([])
  plt.subplots_adjust(top=0.95, bottom=0.3, right=0.95, left=0.05)
  
  ax.annotate(r'$0^o$', (pi-0.1,pi/3.), size=11, color='black')
  ax.annotate(r'$90^o$', (pi/2.-0.2,pi/3.), size=11, color='black')
  ax.annotate(r'$180^o$', (-0.2,pi/3.), size=11, color='black')
  ax.annotate(r'$270^o$', (-pi/2.-0.2,pi/3.), size=11, color='black')
  
  
  ax.annotate(r'$0^o$', (pi-0.1,pi/3.), size=11, color='black')
  ax.annotate(r'$90^o$', (pi/2.-0.2,pi/3.), size=11, color='black')
  ax.annotate(r'$180^o$', (-0.2,pi/3.), size=11, color='black')
  ax.annotate(r'$270^o$', (-pi/2.-0.2,pi/3.), size=11, color='black')  
################################################################ 

  ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
  ax2.set_axis_off()
  ax2.set_xticks([])
  ax2.set_yticks([])
  ax2.xaxis.set_ticks_position('none')
  ax2.yaxis.set_ticks_position('none')
  ax2.annotate(r"$\/\/ V_{LS} \/ < \/ 0 \/\/\/\/ km\/ s^{-1}  $", (0.1,0.2), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$ \/\/\/\/ 0 \/ \leq \/ V_{LS} \/ < \/ 250$', (0.1,0.2-0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$\/ 250 \/ \leq \/ V_{LS} \/ < \/ 500$', (0.1,0.2-2*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$\/ 500 \/ \leq \/ V_{LS} \/ < \/ 750$', (0.1,0.2-3*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$\/ 750 \/ \leq \/ V_{LS} \/ < \/ 1000$', (0.30,0.2), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1000 \/ \leq \/ V_{LS} \/ < \/ 1250$', (0.30,0.2-0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1250 \/ \leq \/ V_{LS} \/ < \/ 1500$', (0.30,0.2-2*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1500 \/ \leq \/ V_{LS} \/ < \/ 1750$', (0.30,0.2-3*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1750 \/ \leq \/ V_{LS} \/ < \/ 2000$', (0.5,0.2), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$2000 \/ \leq \/ V_{LS} \/ < \/ 2500$', (0.5,0.2-0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$2500 \/ \leq \/ V_{LS} \/ < \/ 3000$', (0.5,0.2-2*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$3000 \/ \leq \/ V_{LS} \/ < \/ 3500$', (0.5,0.2-3*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$3500 \/ \leq \/ V_{LS}$', (0.7,0.2), xycoords='figure fraction', size=12, color='black')
  
  p = 0
  for m in [0.1,0.3,0.5]:
      for n in [0.2,0.2-0.05,0.2-2*0.05,0.2-3*0.05]:
          ax2.add_patch(patches.Rectangle((m-0.05, n), 0.03, 0.03, color=col_pallet[p]))
          p+=1
  ax2.add_patch(patches.Rectangle((0.7-0.05, 0.2), 0.03, 0.03, color=col_pallet[p]))
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
    
    mytable = file_deliver('all.iter.2.v43.group')
    
    
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
    
    
    indices = np.where(flag==0)
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
    
    for i in range(len(id_gal)):
      x, y = xymap_aitoff(sgl_gal[i], sgb_gal[i])
      l, = ax.plot([x], [y], '.', markersize = 1, color=color_table(Vls_gal[i]))
      #l, = ax.plot([x], [y], '.', markersize = 1, color=scalarMap.to_rgba(Vls_gal[i]))
    

    indices = np.where(flag>0)
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


    i = 0 
    N = len(id_gal)
    gr = flag[i]
    while gr == 2:
      

      colour = color_table(Vls_gal[i])
      #colour = scalarMap.to_rgba(Vls_gal[i])
      
      
      Dist = Vls_gal[i]/H0
      if Dist<1: Dist = 1.
      

      r = 180.*atan(R_2t2[i]/Dist)/pi

      
      d_theta = 0.01
      u = np.arange(0,2*pi,d_theta)
      X = u*0.
      Y = u*0.
      
      RA0  = sgl_gal[i]
      DEC0 = sgb_gal[i]
      
      for q in range(0,len(u)):
        X[q] = r*np.cos(u[q]) + RA0
        Y[q] = r*np.sin(u[q]) + DEC0

      OUT = xymap_aitoff_lst(X,Y)
      
      ### Plotting group borders  --- Brent: does not want to see the border of groups
      #if r <= 20 and Vls_gal[i]>0:
        #for XY in OUT:
            #line, = ax.plot(XY[0],XY[1], '-', markersize = 2, color=colour, picker=False) 
            #line.set_dashes([8, 3]) 


      i+=1  
      X = []
      Y = []
      while i<N and flag[i]==1: 

        x, y = xymap_aitoff(sgl_gal[i], sgb_gal[i])
        X.append(x)
        Y.append(y)
        i+=1

        point, = ax.plot(X, Y, 'o', markersize = 1, color=colour, markeredgecolor = colour)

      
      if i<N and flag[i]==2: 
	    gr = 2
      else:
	    break  
################################################################ 

  #plt.show()
  plt.savefig('All_sky_v05.png', dpi=600)  # the best view and size happens in eps format









