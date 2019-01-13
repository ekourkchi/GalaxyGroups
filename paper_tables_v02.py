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
def create_tables():
  
  filee_super = 'all.iter.2.v44.supergroup'
    
  try:
    intable_super = np.genfromtxt(filee_super , delimiter='|', filling_values="-100000", names=True, dtype=None )
  except:
    print "\n[Error] The catalog \""+filee_super+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)  
  

  pgc_super          = intable_super['pgc']
  flag_super         = intable_super['flag']
  sgl_super          = intable_super['sgl']
  sgb_super          = intable_super['sgb']
  gl_super           = intable_super['gl']
  gb_super           = intable_super['gb']  
  ra_super           = intable_super['ra']
  dec_super          = intable_super['dec']
  Ks_super           = intable_super['Ks']
  logK_super         = intable_super['logK']
  Vls_super          = intable_super['Vls']
  mDist_super        = intable_super['mDist']
  mDistErr_super     = intable_super['mDistErr']
  PGC1_super         = intable_super['nest']
  r1t_lum_super      = intable_super['r1t_lum']
  No_Galaxies_super  = intable_super['No_Galaxies']
  Mv_lum_super       = intable_super['Mv_lum']
  
  
  pgc_p_super = [] 
  mass_p_super = [] 
  pgc_tmp_super =  pgc_super[0]
  mass_tmp_super =  Mv_lum_super[0]
  
  for i in range(len(pgc_super)):
      if flag_super[i] == 5 or flag_super[i] == 0 or flag_super[i] == 2:
         pgc_tmp_super  = pgc_super[i]
         mass_tmp_super = Mv_lum_super[i]
         
      if flag_super[i] == 5:
         PGC1_super[i] = PGC1_super[i+1]
         
      
      pgc_p_super.append(extractPGC(pgc_tmp_super, grp=False, supergrp=True))
      mass_p_super.append(mass_tmp_super)
  
  pgc_p_super = np.asarray(pgc_p_super)    # PGC+   supergroup ID
  mass_p_super = np.asarray(mass_p_super)
  
  
  for i in range(len(Ks_super)): 
      if Ks_super[i] == 0:
          Ks_super[i] = None
          

  for i in range(len(logK_super)): 
      if logK_super[i] == 0:
          logK_super[i] = None
                                  
  for i in range(len(mDist_super)): 
      if mDist_super[i] == 0:
          mDist_super[i] = None
          mDistErr_super[i] = None  
  
    
  indices            = np.where(flag_super==5)
  pgc_p_super_t      = pgc_p_super[indices]
  PGC1_super_t       = PGC1_super[indices]
  sgl_super_t        = sgl_super[indices]
  sgb_super_t        = sgb_super[indices]
  logK_super_t       = logK_super[indices]
  Vls_super_t        = Vls_super[indices]
  mDist_super_t      = mDist_super[indices]
  mDistErr_super_t   = mDistErr_super[indices]
  r1t_lum_super_t    = r1t_lum_super[indices]
  Mv_lum_super_t     = Mv_lum_super[indices]
  
  
  indices            = np.argsort(Mv_lum_super_t)
  indices = indices[::-1]
  pgc_p_super_t      = pgc_p_super_t[indices]
  PGC1_super_t       = PGC1_super_t[indices]
  sgl_super_t        = sgl_super_t[indices]
  sgb_super_t        = sgb_super_t[indices]
  logK_super_t       = logK_super_t[indices]
  Vls_super_t        = Vls_super_t[indices]
  mDist_super_t      = mDist_super_t[indices]
  mDistErr_super_t   = mDistErr_super_t[indices]
  r1t_lum_super_t    = r1t_lum_super_t[indices]
  Mv_lum_super_t     = Mv_lum_super_t[indices]  
  
  
  myTable = Table()
  myTable.add_column(Column(data=pgc_p_super_t, name='PGC1+'))
  myTable.add_column(Column(data=sgl_super_t, name='SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb_super_t, name='SGB', format='%0.4f'))
  myTable.add_column(Column(data=logK_super_t, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vls_super_t, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=mDist_super_t, name='D', format='%0.2f'))
  myTable.add_column(Column(data=mDistErr_super_t*100, name='eD', format='%0.0f'))
  myTable.add_column(Column(data=r1t_lum_super_t, name='r1t', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_lum_super_t), name='Mass_lum', format='%0.3f'))


  myTable.write('SGroup_table_v44.csv', format='ascii.fixed_width',delimiter=',', bookend=False)   
  
  
########################################################################## 
########################################################################## 
     
  filee = 'all.iter.2.v44.group'
  
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
  R_2t2_dyn    = intable['R2t_dyn']
  Rg_dyn       = intable['Rg_dyn']
  PGC1         = intable['nest']
  sigmaP_lum   = intable['sigmaP_lum']    
  sigmaP_dyn   = intable['sigmaP_dyn']    
  Mv_lum       = intable['Mv_lum']
  Mv_dyn       = intable['Mv_dyn']
  No_Galaxies  = intable['No_Galaxies']
  objname      = intable['objname']
  

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
  R_2t2_dyn_grp    = R_2t2_dyn[indices]
  Rg_dyn_grp       = Rg_dyn[indices]
  PGC1_grp         = PGC1[indices]
  sigmaP_lum_grp   = sigmaP_lum[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn[indices]  
  Mv_lum_grp       = Mv_lum[indices]
  Mv_dyn_grp       = Mv_dyn[indices]
  No_Galaxies_grp  = No_Galaxies[indices] 
  objnames_grp     = objname[indices] 
  
  indices = np.argsort(logK_grp)
  #indices = indices[::-1]
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
  R_2t2_dyn_grp    = R_2t2_dyn_grp[indices]
  Rg_dyn_grp       = Rg_dyn_grp[indices]
  PGC1_grp         = PGC1_grp[indices]
  sigmaP_lum_grp   = sigmaP_lum_grp[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn_grp[indices]  
  Mv_lum_grp       = Mv_lum_grp[indices]
  Mv_dyn_grp       = Mv_dyn_grp[indices]
  No_Galaxies_grp  = No_Galaxies_grp[indices]   
  objnames_grp     = objnames_grp[indices]  
  
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

  for i in range(len(Mv_dyn_grp)): 
      if Mv_dyn_grp[i] == 0:
          Mv_dyn_grp[i] = None

  for i in range(len(R_2t2_dyn_grp)): 
      if R_2t2_dyn_grp[i] == 0:
          R_2t2_dyn_grp[i] = None

  for i in range(len(Rg_dyn_grp)): 
      if Rg_dyn_grp[i] == 0:
          Rg_dyn_grp[i] = None


  PGC_p  = []
  Mass_p = []
  for pgc1 in PGC1_grp:
      i = 0 
      while pgc1 != PGC1_super[i]: i+=1
      PGC_p.append(pgc_p_super[i])
      Mass_p.append(mass_p_super[i])
        
  PGC_p = np.asarray(PGC_p)
  Mass_p = np.asarray(Mass_p)
  
  
  indices = np.argsort(Mass_p, kind='mergesort')
  indices = indices[::-1]
  Mass_p           = Mass_p[indices]
  PGC_p            = PGC_p[indices]
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
  R_2t2_dyn_grp    = R_2t2_dyn_grp[indices]
  Rg_dyn_grp       = Rg_dyn_grp[indices]
  PGC1_grp         = PGC1_grp[indices]
  sigmaP_lum_grp   = sigmaP_lum_grp[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn_grp[indices]  
  Mv_lum_grp       = Mv_lum_grp[indices]
  Mv_dyn_grp       = Mv_dyn_grp[indices]
  No_Galaxies_grp  = No_Galaxies_grp[indices]    
  
      
  myTable = Table()
  myTable.add_column(Column(data=PGC1_grp, name='PGC1'))
  myTable.add_column(Column(data=PGC_p, name='PGC1+'))
  myTable.add_column(Column(data=No_Galaxies_grp, name='Mem'))
  myTable.add_column(Column(data=gl_grp, name='Glon', format='%0.4f'))
  myTable.add_column(Column(data=gb_grp, name='Glat', format='%0.4f'))  
  myTable.add_column(Column(data=sgl_grp, name='SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb_grp, name='SGB', format='%0.4f'))
  myTable.add_column(Column(data=Ks_grp, name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK_grp, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vhelio_grp, name='Vh', format='%0.0f'))
  myTable.add_column(Column(data=Vls_grp, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=mDist_grp, name='D', format='%0.2f'))
  myTable.add_column(Column(data=mDistErr_grp*100, name='eD', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_lum_grp, name='Sigma_L', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_dyn_grp, name='Sigma_V', format='%0.0f'))
  myTable.add_column(Column(data=R_2t2_grp, name='R2t_lum', format='%0.3f'))
  myTable.add_column(Column(data=Rg_dyn_grp, name='Rg_dyn', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_lum_grp), name='Mass_lum', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_dyn_grp), name='Mass_dyn', format='%0.3f'))

  myTable.write('Group_table_v44.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  
    #########################################################################

  
  indices          = np.where(flag<2)
  pgc_gal          = pgc[indices]
  flag_gal         = flag[indices]
  sgl_gal          = sgl[indices]
  sgb_gal          = sgb[indices]
  gl_gal           = gl[indices]
  gb_gal           = gb[indices]
  ra_gal           = ra[indices]
  dec_gal          = dec[indices]
  Ks_gal           = Ks[indices]
  Ty_gal           = Ty[indices]
  B_mag_gal        = B_mag[indices]
  logK_gal         = logK[indices]
  Vls_gal          = Vls[indices]
  Vhelio_gal       = Vhelio[indices]
  dcf2_gal         = dcf2[indices]
  ed_gal           = ed[indices]
  mDist_gal        = mDist[indices]
  mDistErr_gal     = mDistErr[indices]
  R_2t2_gal        = R_2t2[indices]
  PGC1_gal         = PGC1[indices]
  sigmaP_lum_gal   = sigmaP_lum[indices]   
  sigmaP_dyn_gal   = sigmaP_dyn[indices]  
  Mv_lum_gal       = Mv_lum[indices]
  No_Galaxies_gal  = No_Galaxies[indices] 
  objname_gal      = objname[indices] 
 
  for i in range(len(Ty_gal)): 
      if Ty_gal[i] == -100000.0:
          Ty_gal[i] = None
  
  for i in range(len(B_mag_gal)): 
      if B_mag_gal[i] == -100000.0:
          B_mag_gal[i] = None          
   
  for i in range(len(Ks_gal)): 
      if Ks_gal[i] == 0:
          Ks_gal[i] = None
          

  for i in range(len(logK_gal)): 
      if logK_gal[i] == 0:
          logK_gal[i] = None
                                  
  for i in range(len(dcf2_gal)): 
      if dcf2_gal[i] == 0:
          dcf2_gal[i] = None
          ed_gal[i] = None
  
  
  
  
  
  
  
  #########################################################################
  
  
  myTable = Table()
  
  myTable.add_column(Column(data=pgc_gal, name='PGC'))
  myTable.add_column(Column(data=objname_gal, name='Name', dtype='S35'))
  myTable.add_column(Column(data=ra_gal, name='R.A.', format='%0.4f'))
  myTable.add_column(Column(data=dec_gal, name='Dec.', format='%0.4f'))
  myTable.add_column(Column(data=gl_gal, name='Glon', format='%0.4f'))
  myTable.add_column(Column(data=gb_gal, name='Glat', format='%0.4f'))
  myTable.add_column(Column(data=sgl_gal, name='SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb_gal, name='SGB', format='%0.4f'))
  myTable.add_column(Column(data=Ty_gal, name='Ty', format='%0.1f'))
  myTable.add_column(Column(data=B_mag_gal, name='B', format='%0.2f'))
  myTable.add_column(Column(data=Ks_gal, name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK_gal, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vhelio_gal, name='Vh', format='%0.0f'))
  myTable.add_column(Column(data=Vls_gal, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=dcf2_gal, name='D', format='%0.2f'))
  myTable.add_column(Column(data=ed_gal*100, name='eD', format='%0.0f'))
  myTable.add_column(Column(data=PGC1_gal, name='PGC1'))

  myTable.write('Gal_table_v44.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  
  PGC1_galgrp = []
  PGC_p_gal = []
  No_Galaxies_galgrp = []
  gl_galgrp = []
  gb_galgrp = []
  sgl_galgrp = []
  sgb_galgrp = []
  Ks_galgrp = []
  logK_galgrp = []
  Vhelio_galgrp = []
  Vls_galgrp = []
  mDist_galgrp = []
  mDistErr_galgrp = []
  sigmaP_lum_galgrp = []
  sigmaP_dyn_galgrp = []
  R_2t2_galgrp = []
  Rg_dyn_galgrp = []
  Mv_lum_galgrp = []
  Mv_dyn_galgrp = []
  
  for pgc1 in PGC1_gal:
      i = 0 
      while PGC1_grp[i] != pgc1: i+=1
      PGC1_galgrp.append(PGC1_grp[i])
      PGC_p_gal.append(PGC_p[i])
      No_Galaxies_galgrp.append(No_Galaxies_grp[i])
      gl_galgrp.append(gl_grp[i])
      gb_galgrp.append(gb_grp[i])
      sgl_galgrp.append(sgl_grp[i])
      sgb_galgrp.append(sgb_grp[i])
      Ks_galgrp.append(Ks_grp[i])
      logK_galgrp.append(logK_grp[i])
      Vhelio_galgrp.append(Vhelio_grp[i])
      Vls_galgrp.append(Vls_grp[i])
      mDist_galgrp.append(mDist_grp[i])
      mDistErr_galgrp.append(mDistErr_grp[i])
      sigmaP_lum_galgrp.append(sigmaP_lum_grp[i])
      sigmaP_dyn_galgrp.append(sigmaP_dyn_grp[i])
      R_2t2_galgrp.append(R_2t2_grp[i])
      Rg_dyn_galgrp.append(Rg_dyn_grp[i])
      Mv_lum_galgrp.append(Mv_lum_grp[i])
      Mv_dyn_galgrp.append(Mv_dyn_grp[i])      


  PGC1_galgrp = np.asarray(PGC1_galgrp)
  PGC_p_gal = np.asarray(PGC_p_gal)
  No_Galaxies_galgrp = np.asarray(No_Galaxies_galgrp)
  gl_galgrp = np.asarray(gl_galgrp)
  gb_galgrp = np.asarray(gb_galgrp)
  sgl_galgrp = np.asarray(sgl_galgrp)
  sgb_galgrp = np.asarray(sgb_galgrp)
  Ks_galgrp = np.asarray(Ks_galgrp)
  logK_galgrp = np.asarray(logK_galgrp)
  Vhelio_galgrp = np.asarray(Vhelio_galgrp)
  Vls_galgrp = np.asarray(Vls_galgrp)
  mDist_galgrp = np.asarray(mDist_galgrp)
  mDistErr_galgrp = np.asarray(mDistErr_galgrp)
  sigmaP_lum_galgrp = np.asarray(sigmaP_lum_galgrp)
  sigmaP_dyn_galgrp = np.asarray(sigmaP_dyn_galgrp)
  R_2t2_galgrp = np.asarray(R_2t2_galgrp)
  Rg_dyn_galgrp = np.asarray(Rg_dyn_galgrp)
  Mv_lum_galgrp = np.asarray(Mv_lum_galgrp)
  Mv_dyn_galgrp = np.asarray(Mv_dyn_galgrp)

  
  #myTable.add_column(Column(data=PGC1_galgrp, name='g_PGC1'))
  myTable.add_column(Column(data=PGC_p_gal, name='PGC1+'))
  myTable.add_column(Column(data=No_Galaxies_galgrp, name='Mem.'))
  myTable.add_column(Column(data=gl_galgrp, name='g_Glat', format='%0.4f'))
  myTable.add_column(Column(data=gb_galgrp, name='g_Glon', format='%0.4f'))  
  myTable.add_column(Column(data=sgl_galgrp, name='g_SGL', format='%0.4f'))
  myTable.add_column(Column(data=sgb_galgrp, name='g_SGB', format='%0.4f'))
  myTable.add_column(Column(data=Ks_galgrp, name='g_Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK_galgrp, name='g_logK', format='%0.2f'))
  myTable.add_column(Column(data=Vhelio_galgrp, name='g_Vh', format='%0.0f'))
  myTable.add_column(Column(data=Vls_galgrp, name='g_Vls', format='%0.0f'))
  myTable.add_column(Column(data=mDist_galgrp, name='g_D', format='%0.2f'))
  myTable.add_column(Column(data=mDistErr_galgrp*100, name='g_eD', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_lum_galgrp, name='Sigma_L', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_dyn_galgrp, name='Sigma_V', format='%0.0f'))
  myTable.add_column(Column(data=R_2t2_galgrp, name='R2t_lum', format='%0.3f'))
  myTable.add_column(Column(data=Rg_dyn_galgrp, name='Rg_dyn', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_lum_galgrp), name='Mass_lum', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_dyn_galgrp), name='Mass_dyn', format='%0.3f'))



  myTable.write('EDD_table_v44.csv', format='ascii.fixed_width',delimiter='|', bookend=False)    
#########################################################################

########################################################################## 
  
if __name__ == '__main__':


  #create_tables()
#########################################################################

  intable = np.genfromtxt('Gal_table_v44.csv' , delimiter=',', filling_values=None, names=True, dtype=None)
  
  glat      = intable['Glat']
  dist      = intable['D']
  
  print 'total: ', len(glat)
  ind = np.where(dist>=0)
  dist_total = dist[ind]
  print 'total dist: ', len(dist_total)
  
  ind = np.where(glat>=0)
  glat_north = glat[ind]
  dist_north = dist[ind]
  print 'north all: ', len(glat_north)
  
  ind = np.where(dist_north>=0)
  glat_north_d = glat_north[ind]
  print 'north dist: ', len(glat_north_d)
    
  ind = np.where(glat<0)
  glat_south = glat[ind]
  dist_south = dist[ind]
  print 'south all: ', len(glat_south)
  
  ind = np.where(dist_south>=0)
  glat_south_d = glat_south[ind]
  print 'south dist: ', len(glat_south_d)
#########################################################################
    
  intable = np.genfromtxt('Group_table_v44.csv' , delimiter=',', filling_values=None, names=True, dtype=None)
    
  glat      = intable['Glat']   
  mem       = intable['Mem']  
  dist      = intable['D'] 
   
  ind = np.where(mem>=2)
  glat_tot = glat[ind]
  mem_tot  = mem[ind]
   
  print 'grp2 tot: ', len(glat_tot), sum(mem_tot)
   
  ind = np.where(glat_tot>=0)
  glat_tot_north = glat_tot[ind]
  mem_tot_north  = mem_tot[ind]  
  print 'grp2 north: ', len(glat_tot_north), sum(mem_tot_north)
   
  ind = np.where(glat_tot<0)
  glat_tot_south = glat_tot[ind]
  mem_tot_south  = mem_tot[ind]   
  print 'grp2 south: ', len(glat_tot_south), sum(mem_tot_south)
   

   
  ind = np.where(mem>=5)
  glat_tot = glat[ind]
  mem_tot  = mem[ind]
   
  print 'grp5 tot: ', len(glat_tot), sum(mem_tot)
   
  ind = np.where(glat_tot>=0)
  glat_tot_north = glat_tot[ind]
  mem_tot_north  = mem_tot[ind]  
  print 'grp5 north: ', len(glat_tot_north), sum(mem_tot_north)
   
  ind = np.where(glat_tot<0)
  glat_tot_south = glat_tot[ind]
  mem_tot_south  = mem_tot[ind]   
  print 'grp5 south: ', len(glat_tot_south), sum(mem_tot_south)   
   
  glat      = intable['Glat']   
  mem       = intable['Mem']  
  ind = np.where(dist>=0)
  glat = glat[ind]
  mem  = mem[ind]
   
   
  ind = np.where(mem>=2)
  glat_tot = glat[ind]
  mem_tot  = mem[ind]
   
  print 'd grp tot: ', len(glat_tot), sum(mem_tot)
   
  ind = np.where(glat_tot>=0)
  glat_tot_north = glat_tot[ind]
  mem_tot_north  = mem_tot[ind]  
  print 'd grp north: ', len(glat_tot_north), sum(mem_tot_north)
   
  ind = np.where(glat_tot<0)
  glat_tot_south = glat_tot[ind]
  mem_tot_south  = mem_tot[ind]   
  print 'd grp south: ', len(glat_tot_south), sum(mem_tot_south)   
   
  glat      = intable['Glat'] 
  mem       = intable['Mem']  
  ind = np.where(mem==1)
  glat_tot = glat[ind]
   
  print 'grp=1 tot: ', len(glat_tot)   
