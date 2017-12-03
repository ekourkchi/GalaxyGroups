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
       
             s = str(number)
             t.set_text(u''+s+u'\xb0')

########################################################################## 
def draw_patch_poly( map, x0, y0, GL_vertex=None, GB_vertex=None, color='blue'):
    
    if GL_vertex==None or GB_vertex==None:
      return
    
    lons = GL_vertex
    lats = GB_vertex
    
    x, y = map( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, closed=True, fill=True, color=color)
    plt.gca().add_patch(poly)
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
    

##########################################################################   
if __name__ == '__main__':
  
  fig = plt.figure(figsize=(8, 6), dpi=100)
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
  map = Basemap(projection='npaeqd',boundinglat=10,lon_0=270,resolution='h', ax=ax)


  ## draw parallels and meridians.
  paral = map.drawparallels(np.arange(-80.,81.,20.))
  merid = map.drawmeridians(np.arange(-180.,181.,20.), labels=[1,1,1,1])
  setcolor(paral,'r')
  setcolor(merid,'b')
  

  
##########################################################################

  a = pyfits.open('EBV.nside.8.gal.fits')    # source: dust_test.pro
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
  
  for i in range(len(GB)):
    
    GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
    GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]
    
    val = EBV_val(EBV[i])
    colorVal = scalarMap.to_rgba(val)
    if GL[i]>=0 and GL[i]<90 and GB[i]>-40 and GB[i]<80:
     draw_patch_poly(map, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal )
     
  for i in range(len(GB)):
    GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
    GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
    val = EBV_val(EBV[i])
    colorVal = scalarMap.to_rgba(val)
    if GL[i]<180 and GL[i]>=90 and GB[i]>-40 and GB[i]<80:
     draw_patch_poly(map, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal )
  
  for i in range(len(GB)):
    GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
    GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
    val = EBV_val(EBV[i])
    colorVal = scalarMap.to_rgba(val)
    if GL[i]>=180 and  GL[i]<270 and GB[i]>-40 and GB[i]<80:
     draw_patch_poly(map, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal )  
     
  for i in range(len(GB)):
    GL_vertex = [GL_v1[i], GL_v2[i], GL_v3[i], GL_v4[i]]
    GB_vertex = [GB_v1[i], GB_v2[i], GB_v3[i], GB_v4[i]]    
    val = EBV_val(EBV[i])
    colorVal = scalarMap.to_rgba(val)
    if GL[i]>=270 and GB[i]>-40 and GB[i]<80:
     draw_patch_poly(map, GL[i], GB[i], GL_vertex=GL_vertex, GB_vertex=GB_vertex, color=colorVal ) 
     
##########################################################################
  
  lat = [50]
  lon = [90]
  x,y = map(lon, lat)
  map.plot(x,y,'o')
  
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
      x, y = map([sgl_gal[i]], [sgb_gal[i]])
      map.plot(x, y, '.', markersize = 1, color=color_table(Vls_gal[i]))
      #ax.plot([x], [y], '.', markersize = 1, color=scalarMap.to_rgba(Vls_gal[i]))
    





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
  
  