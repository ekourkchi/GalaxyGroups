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


#################################################
def extractPGC_lst(id, grp=False, supergrp=False):
  
  idl = []
  for i in id: idl.append(extractPGC(i, grp=grp, supergrp=supergrp))
  
  return np.asarray(idl)
      
      
  
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
  
if __name__ == '__main__':
     
  filee = 'all.iter.2.v43.group'
  
  try:
    intable = np.genfromtxt(filee , delimiter='|', filling_values="-100000", names=True, dtype=None )
  except:
    print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)  
  

  pgc          = intable['pgc']
  flag         = intable['flag']
  sgl          = intable['sgl']
  sgb          = intable['sgb']
  gl           = intable['gl']
  gb           = intable['gb']  
  ra           = intable['ra']
  dec          = intable['dec']
  Ty           = intable['Ty']
  B_mag        = intable['B_mag']
  Ks           = intable['Ks']
  logK         = intable['logK']
  Vls          = intable['Vls']
  Vhelio       = intable['Vhelio']
  dcf2         = intable['dcf2']
  ed           = intable['ed']
  mDist        = intable['mDist']
  mDistErr     = intable['mDistErr']
  R_2t2        = intable['R2t_lum']
  PGC1         = intable['nest']
  sigmaP_lum   = intable['sigmaP_lum']    
  sigmaP_dyn   = intable['sigmaP_dyn']    
  Mv_lum       = intable['Mv_lum']
  No_Galaxies  = intable['No_Galaxies']
  

  indices = np.where(flag!=1)
  pgc_grp          = pgc[indices]
  flag_grp         = flag[indices]
  sgl_grp          = sgl[indices]
  sgb_grp          = sgb[indices]
  gl_grp           = gl[indices]
  gb_grp           = gb[indices]
  ra_grp           = ra[indices]
  dec_grp          = dec[indices]
  Ks_grp           = Ks[indices]
  logK_grp         = logK[indices]
  Vls_grp          = Vls[indices]
  Vhelio_grp       = Vhelio[indices]
  dcf2_grp         = dcf2[indices]
  ed_grp           = ed[indices]
  mDist_grp        = mDist[indices]
  mDistErr_grp     = mDistErr[indices]
  R_2t2_grp        = R_2t2[indices]
  PGC1_grp         = PGC1[indices]
  sigmaP_lum_grp   = sigmaP_lum[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn[indices]  
  Mv_lum_grp       = Mv_lum[indices]
  No_Galaxies_grp  = No_Galaxies[indices] 
  
  indices = np.argsort(logK_grp)
  indices = indices[::-1]
  pgc_grp          = pgc_grp[indices]
  flag_grp         = flag_grp[indices]
  sgl_grp          = sgl_grp[indices]
  sgb_grp          = sgb_grp[indices]
  gl_grp           = gl_grp[indices]
  gb_grp           = gb_grp[indices] 
  ra_grp           = ra_grp[indices]
  dec_grp          = dec_grp[indices]
  Ks_grp           = Ks_grp[indices]
  logK_grp         = logK_grp[indices]
  Vls_grp          = Vls_grp[indices]
  Vhelio_grp       = Vhelio_grp[indices]
  dcf2_grp         = dcf2_grp[indices]
  ed_grp           = ed_grp[indices]
  mDist_grp        = mDist_grp[indices]
  mDistErr_grp     = mDistErr_grp[indices]
  R_2t2_grp        = R_2t2_grp[indices]
  PGC1_grp         = PGC1_grp[indices]
  sigmaP_lum_grp   = sigmaP_lum_grp[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn_grp[indices]  
  Mv_lum_grp       = Mv_lum_grp[indices]
  No_Galaxies_grp  = No_Galaxies_grp[indices]   
  
  for i in range(len(Ks_grp)): 
      if Ks_grp[i] == 0:
          Ks_grp[i] = None
          

  for i in range(len(logK_grp)): 
      if logK_grp[i] == 0:
          logK_grp[i] = None
                                  
  for i in range(len(mDist_grp)): 
      if mDist_grp[i] == 0:
          mDist_grp[i] = None
          mDistErr_grp[i] = None
  
      
  myTable = Table()
  myTable.add_column(Column(data=PGC1_grp, name='PGC1'))
  #myTable.add_column(Column(data=flag_grp, name='Flag'))
  myTable.add_column(Column(data=sgl_grp, name='SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb_grp, name='SGB', format='%0.4f'))
  myTable.add_column(Column(data=No_Galaxies_grp, name='Mem.'))
  myTable.add_column(Column(data=Ks_grp, name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK_grp, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vls_grp, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=mDist_grp, name='D', format='%0.2f'))
  myTable.add_column(Column(data=mDistErr_grp*100, name='eD', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_lum_grp, name='Sigma_l', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_dyn_grp, name='Sigma_d', format='%0.0f'))
  myTable.add_column(Column(data=R_2t2_grp, name='R2t', format='%0.3f'))
  myTable.add_column(Column(data=Mv_lum_grp, name='Mass', format='%1.2e'))

  myTable.write('group_table_v43.csv', format='ascii.fixed_width',delimiter='|', bookend=False)
  
  
  
  
  indices = np.where(flag<2)
  pgc_grp          = pgc[indices]
  flag_grp         = flag[indices]
  sgl_grp          = sgl[indices]
  sgb_grp          = sgb[indices]
  gl_grp           = gl[indices]
  gb_grp           = gb[indices]
  ra_grp           = ra[indices]
  dec_grp          = dec[indices]
  Ks_grp           = Ks[indices]
  Ty_grp           = Ty[indices]
  B_mag_grp        = B_mag[indices]
  logK_grp         = logK[indices]
  Vls_grp          = Vls[indices]
  Vhelio_grp       = Vhelio[indices]
  dcf2_grp         = dcf2[indices]
  ed_grp           = ed[indices]
  mDist_grp        = mDist[indices]
  mDistErr_grp     = mDistErr[indices]
  R_2t2_grp        = R_2t2[indices]
  PGC1_grp         = PGC1[indices]
  sigmaP_lum_grp   = sigmaP_lum[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn[indices]  
  Mv_lum_grp       = Mv_lum[indices]
  No_Galaxies_grp  = No_Galaxies[indices] 
 
 
  for i in range(len(Ty_grp)): 
      if Ty_grp[i] == -100000.0:
          Ty_grp[i] = None
  
  for i in range(len(B_mag_grp)): 
      if B_mag_grp[i] == -100000.0:
          B_mag_grp[i] = None          
   
  for i in range(len(Ks_grp)): 
      if Ks_grp[i] == 0:
          Ks_grp[i] = None
          

  for i in range(len(logK_grp)): 
      if logK_grp[i] == 0:
          logK_grp[i] = None
                                  
  for i in range(len(dcf2_grp)): 
      if dcf2_grp[i] == 0:
          dcf2_grp[i] = None
          ed_grp[i] = None
  
  
  myTable = Table()
  myTable.add_column(Column(data=PGC1_grp, name='PGC1'))
  myTable.add_column(Column(data=pgc_grp, name='PGC'))
  #myTable.add_column(Column(data=flag_grp, name='Flag'))
  myTable.add_column(Column(data=ra_grp, name='R.A.'))
  myTable.add_column(Column(data=dec_grp, name='Dec.'))
  myTable.add_column(Column(data=gl_grp, name='Glon'))
  myTable.add_column(Column(data=gb_grp, name='Glat'))
  
  myTable.add_column(Column(data=sgl_grp, name='SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb_grp, name='SGB', format='%0.4f'))
  myTable.add_column(Column(data=Ty_grp, name='Ty'))
  myTable.add_column(Column(data=B_mag_grp, name='B'))
  myTable.add_column(Column(data=Ks_grp, name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK_grp, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vls_grp, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=dcf2_grp, name='D', format='%0.2f'))
  myTable.add_column(Column(data=ed_grp*100, name='eD', format='%0.0f'))


  myTable.write('galaxy_table_v43.csv', format='ascii.fixed_width',delimiter='|', bookend=False)    
#########################################################################
  
  
  filee = 'all.iter.2.v43.supergroup'
    
  try:
    intable = np.genfromtxt(filee , delimiter='|', filling_values="-100000", names=True, dtype=None )
  except:
    print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)  
  

  pgc          = intable['pgc']
  flag         = intable['flag']
  sgl          = intable['sgl']
  sgb          = intable['sgb']
  gl           = intable['gl']
  gb           = intable['gb']  
  ra           = intable['ra']
  dec          = intable['dec']
  Ks           = intable['Ks']
  logK         = intable['logK']
  Vls          = intable['Vls']
  mDist        = intable['mDist']
  mDistErr     = intable['mDistErr']
  PGC1         = intable['nest']
  r1t_lum    = intable['r1t_lum']
  No_Galaxies  = intable['No_Galaxies']
  
  
  pgc_p = [] 
  pgc_tmp =  pgc[0]
  
  for i in range(len(pgc)):
      if flag[i] == 5 or flag[i] == 0:
         pgc_tmp  = pgc[i]
         
      if flag[i] == 5:
         PGC1[i] = PGC1[i+1]
      
      pgc_p.append(extractPGC(pgc_tmp, grp=False, supergrp=True))
  
  pgc_p = np.asarray(pgc_p)
  
  
  for i in range(len(Ks)): 
      if Ks[i] == 0:
          Ks[i] = None
          

  for i in range(len(logK)): 
      if logK[i] == 0:
          logK[i] = None
                                  
  for i in range(len(mDist)): 
      if mDist[i] == 0:
          mDist[i] = None
          mDistErr[i] = None  
  
  
  myTable = Table()
  myTable.add_column(Column(data=pgc_p, name='PGC+'))
  myTable.add_column(Column(data=flag, name='Flag'))
  myTable.add_column(Column(data=PGC1, name='PGC1'))
  myTable.add_column(Column(data=No_Galaxies, name='Mem.'))
  myTable.add_column(Column(data=sgl, name='SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb, name='SGB', format='%0.4f'))
  myTable.add_column(Column(data=Ks, name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vls, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=mDist, name='D', format='%0.2f'))
  myTable.add_column(Column(data=mDistErr*100, name='eD', format='%0.0f'))
  myTable.add_column(Column(data=r1t_lum, name='r1t', format='%0.3f'))


  myTable.write('supergroups_table_v43.csv', format='ascii.fixed_width',delimiter='|', bookend=False)   
  
  
  
  
  
  
  
    
    
    
    
    
    
