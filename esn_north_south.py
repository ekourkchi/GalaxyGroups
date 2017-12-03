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
from mpl_toolkits.basemap import Basemap

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
def draw_patch_poly( map, ax, x0, y0, GL_vertex=None, GB_vertex=None, color='blue'):
    
    if GL_vertex==None or GB_vertex==None:
      return
    
    lons = GL_vertex
    lats = GB_vertex
    
    x, y = map( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, closed=True, fill=True, color=color)
    ax.add_patch(poly)
########################################################################## 
def EBV_val(EBV):

    if EBV>2.0: val = 3
    elif EBV>1.0: val = 2.5
    elif EBV>0.5: val= 2.0
    elif EBV<=0.5 and EBV>0.1: val= 1.0
    elif EBV<0.1: val= 0
    else: val= 0
    
    return val
##########################################################################
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
    
# **************************************
def extractPGC(id, grp=False, supergrp=False):
  
  if not grp and not supergrp:
    return id
  
  
  if grp:
    pgc = int(id)%100000000
  
  
  if supergrp:
    grp = int(id)%10000000000
    pgc = int(grp)%100000000
  
  return pgc

# **************************************
##########################################################################   
if __name__ == '__main__':

  if len(sys.argv) < 2:
    print "\nEnter the cluster name as the input ..." 
    print "\nexample: \n > python " + str(sys.argv[0])+ " manual"
    print "\nPossible options: north, south" 
    print "Use north/south for whole sky.\n"
    sys.exit(1)

  cluster = str(sys.argv[1])
  
  
  fig = plt.figure(figsize=(10, 7), dpi=100)
  ax = fig.add_subplot(111)
  plt.subplots_adjust(top=0.95, bottom=0.05, right=0.65, left=0.07)
  
  if cluster == 'north' :
    map = Basemap(projection='npstere',boundinglat=10,lon_0=270,resolution='h', ax=ax)
  if cluster == 'south' :
    map = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='h', ax=ax)

  ## draw parallels and meridians.
  paral = map.drawparallels(np.arange(-80.,81.,20.))
  merid = map.drawmeridians(np.arange(-180.,181.,20.), labels=[1,1,1,1])
  setcolor(paral,'b')
  setcolor(merid,'black')
  
################################################################ 

  ax2 = plt.axes([0,0,1,1], axisbg=(1,1,1,0))
  ax2.set_axis_off()
  ax2.set_xticks([])
  ax2.set_yticks([])
  ax2.xaxis.set_ticks_position('none')
  ax2.yaxis.set_ticks_position('none')
  ax2.annotate(r"$\/\/ V_{ls} \/ < \/ 0 \/\/\/\/ km\/ s^{-1}  $", (0.80,0.800), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$ \/\/\/\/ 0 \/ \leq \/ V_{ls} \/ < \/ 250$', (0.80,0.800-0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$\/ 250 \/ \leq \/ V_{ls} \/ < \/ 500$', (0.80,0.800-2*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$\/ 500 \/ \leq \/ V_{ls} \/ < \/ 750$', (0.80,0.800-3*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$\/ 750 \/ \leq \/ V_{ls} \/ < \/ 1000$', (0.80,0.800-4*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1000 \/ \leq \/ V_{ls} \/ < \/ 1250$', (0.80,0.800-5*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1250 \/ \leq \/ V_{ls} \/ < \/ 1500$', (0.80,0.800-6*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1500 \/ \leq \/ V_{ls} \/ < \/ 1750$', (0.80,0.800-7*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$1750 \/ \leq \/ V_{ls} \/ < \/ 2000$', (0.80,0.800-8*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$2000 \/ \leq \/ V_{ls} \/ < \/ 2500$', (0.80,0.800-9*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$2500 \/ \leq \/ V_{ls} \/ < \/ 3000$', (0.80,0.800-10*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$3000 \/ \leq \/ V_{ls} \/ < \/ 3500$', (0.80,0.800-11*0.05), xycoords='figure fraction', size=12, color='black')
  ax2.annotate(r'$3500 \/ \leq \/ V_{ls}$', (0.80,0.800-12*0.05), xycoords='figure fraction', size=12, color='black')
  
  for p in range(0,len(col_pallet)):
     ax2.add_patch(patches.Rectangle((0.75, 0.800-p*0.05), 0.03, 0.03, color=col_pallet[p]))
  
  if cluster == 'north' :
    ax2.annotate("Galactic - NORTH", (0.75,0.9), xycoords='figure fraction', size=18, color='black')
  
  if cluster == 'south' :
    ax2.annotate("Galactic - SOUTH", (0.75,0.9), xycoords='figure fraction', size=18, color='black')
################################################################ 
################################################################ 
  jet = cm = plt.get_cmap('Greys')
  cNorm  = colors.Normalize(vmin=0, vmax=5)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
  ax3 = fig.add_axes([0.75, 0.1, 0.2, 0.03])
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
##########################################################################

  a = pyfits.open('EBV.nside.64.gal.fits')    # source: dust_test.pro
  d = a[1].data

  GL  = d['GL']
  GB  = d['GB']
  GL_v1  = d['GL_v1']
  GB_v1  = d['GB_v1']
  GL_v2  = d['GL_v2']
  GB_v2  = d['GB_v2']
  GL_v3  = d['GL_v3']
  GB_v3  = d['GB_v3']
  GL_v4  = d['GL_v4']
  GB_v4  = d['GB_v4']  
  EBV  = d['EBV']
  
  jet = cm = plt.get_cmap('Greys')
  cNorm  = colors.Normalize(vmin=0, vmax=5)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)  
  
  if cluster == 'north' :
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]>=0 and GL[i]<90 and GB[i]>-40 and GB[i]<90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)
      
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]<180 and GL[i]>=90 and GB[i]>-40 and GB[i]<90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)
    
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]>=180 and  GL[i]<270 and GB[i]>-40 and GB[i]<90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)  
      
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]>=270 and GB[i]>-40 and GB[i]<90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal) 



  if cluster == 'south' :
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]>=0 and GL[i]<90 and GB[i]<40 and GB[i]>-90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)
      
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]<180 and GL[i]>=90 and GB[i]<40 and GB[i]>-90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)
    
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]>=180 and  GL[i]<270 and GB[i]<40 and GB[i]>-90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)  
      
    for i in range(len(GB)):
      GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
      GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
      val = EBV_val(EBV[i])
      colorVal = scalarMap.to_rgba(val)
      if GL[i]>=270 and GB[i]<40 and GB[i]>-90:
        draw_patch_poly(map, ax, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal)
##########################################################################

  
  if True:
    
    mytable = file_deliver('all.iter.2.v39.group')
    
    
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
        x, y = map([gl_gal[i]], [gb_gal[i]])
        map.plot(x, y, '.', markersize = 2, color=color_table(Vls_gal[i]))
        #ax.plot([x], [y], '.', markersize = 1, color=scalarMap.to_rgba(Vls_gal[i]))
    


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
      
      RA0  = gl_gal[i]
      DEC0 = gb_gal[i]
      
      #for q in range(0,len(u)):
	#X[q] = r*np.cos(u[q]) + RA0
	#Y[q] = r*np.sin(u[q]) + DEC0

      #X, Y = map(X,Y)
      #if extractPGC(id_gal[i], grp=True) != 2557 and extractPGC(id_gal[i], grp=True) != 5064336 :
	  #line, = ax.plot(X,Y, '-', markersize = 2, color=colour, picker=False) 
	  #line.set_dashes([8, 3]) 


      i+=1  
      X = []
      Y = []
      while i<N and flag[i]==1: 
	X.append(gl_gal[i])
	Y.append(gb_gal[i])
	i+=1
      
      X, Y = map(X,Y)
      point, = ax.plot(X, Y, 'o', markersize = 2, color=colour, markeredgecolor = colour)

      
      if i<N and flag[i]==2: 
	gr = 2
      else:
	break 


##########################################################################

  ## draw tissot's indicatrix to show distortion.
  #ax = plt.gca()
  #for y in np.linspace(map.ymax/20,19*map.ymax/20,10):
      #for x in np.linspace(map.xmax/20,19*map.xmax/20,10):
          #lon, lat = m(x,y,inverse=True)
          #poly = map.tissot(lon,lat,2.5,100,\
                          #facecolor='green',zorder=10,alpha=0.5)
  #plt.title("North Polar Azimuthal Equidistant Projection")
  plt.show()
  
  