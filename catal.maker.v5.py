#!/home/ehsan/Ureka/Ureka/variants/common/bin/python


import sys
import os
import random
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from numpy import cos, sin
import astropy.coordinates as coord
import astropy.units as u
from math import *
from time import time
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py
from itertools import chain
from astropy import coordinates as coord
from astropy import units as unit


# This is the default circular velocity and LSR peculiar velocity of the Sun
# TODO: make this a config item?
VCIRC = 220. # u.km/u.s
VLSR = [10., 5.25, 7.17] # *u.km/u.s
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

# **************************************
# Given the apparent magnitude of a galaxy, m
# it returns the absolute luminosity 
# Vls is the radial velocity of the galaxy
# So, the distance to a galaxy is estimated 
# based on its velocity
def m_logK(m, Vls, d=0):
    
    Mk_sun = 3.28   # K-band
    distance = Vls / H0  # Mpc
    
    if distance < 1:
      distance = 1 # Mpc
    
    if d!=0:
      distance = d
    
    M = m - 5*log10(distance) - 30 + 5
    logK = -0.4 * (M - Mk_sun)
    
    if Vls == 0:
      logK = 0
    
    
    return logK
# **************************************



def vgsr_to_vhel(gl, gb, vgsr, vcirc=VCIRC, vlsr=VLSR):
    """ Convert a radial velocity in the Galactic standard of rest (GSR) to
        a barycentric radial velocity.

        Parameters
        ----------
        coordinate : l and b are galactic coordinates (gl, gb)
        vgsr : :class:`~astropy.units.Quantity`
            GSR line-of-sight velocity.
        vcirc : :class:`~astropy.units.Quantity`
            Circular velocity of the Sun.
        vlsr : :class:`~astropy.units.Quantity`
            Velocity of the Sun relative to the local standard
            of rest (LSR).

        Returns
        -------
        vhel : :class:`~astropy.units.Quantity`
            Radial velocity in a barycentric rest frame.

    """
    l = gl*pi/180.
    b = gb*pi/180.


    #if not isinstance(vgsr, u.Quantity):
        #raise TypeError("vgsr must be a Quantity subclass")

    # compute the velocity relative to the LSR
    lsr = vgsr - vcirc*sin(l)*cos(b)

    # velocity correction for Sun relative to LSR
    v_correct = vlsr[0]*cos(b)*cos(l) + \
        vlsr[1]*cos(b)*sin(l) + \
        vlsr[2]*sin(b)
    vhel = lsr - v_correct

    return vhel

def vhel_to_vgsr(gl, gb, vhel, vcirc=VCIRC, vlsr=VLSR):
    """ Convert a velocity from a heliocentric radial velocity to
        the Galactic standard of rest (GSR).

        Parameters
        ----------
        coordinate : :class:`~astropy.coordinates.SkyCoord`
            An Astropy SkyCoord object or anything object that can be passed
            to the SkyCoord initializer.
        vhel : :class:`~astropy.units.Quantity`
            Barycentric line-of-sight velocity.
        vcirc : :class:`~astropy.units.Quantity`
            Circular velocity of the Sun.
        vlsr : :class:`~astropy.units.Quantity`
            Velocity of the Sun relative to the local standard
            of rest (LSR).

        Returns
        -------
        vgsr : :class:`~astropy.units.Quantity`
            Radial velocity in a galactocentric rest frame.

    """
    l = gl*pi/180.
    b = gb*pi/180.

    if not isinstance(vhel, u.Quantity):
        raise TypeError("vhel must be a Quantity subclass")

    lsr = vhel + vcirc*sin(l)*cos(b)

    # velocity correction for Sun relative to LSR
    v_correct = vlsr[0]*cos(b)*cos(l) + \
        vlsr[1]*cos(b)*sin(l) + \
        vlsr[2]*sin(b)
    vgsr = lsr + v_correct

    return vgsr





def isNaN(num):
    return num != num


def Vh2Vls(el,b, Vh):
  
    alpha = pi / 180.
    cosb = cos(b*alpha)
    sinb = sin(b*alpha)
    cosl = cos(el*alpha)
    sinl = sin(el*alpha)
    
    vls = float(Vh)-26.*cosl*cosb+317.*sinl*cosb-8.*sinb

    
    return vls

### (another "Vlg" has been given by Courteau and van den Bergh; another by Yahil et al.)
### The Vlg used in MKgroups is their own version. 
### The following function just works fine for MK-groups
def Vlg2Vls(el,b, Vlg):
  
    alpha = pi / 180.
    
    cosb = cos(b*alpha)
    sinb = sin(b*alpha)
    cosl = cos(el*alpha)
    sinl = sin(el*alpha)
    
    
    #v = 316  # km/s
    #ll = 93  # deg
    #bb = -4  # deg
    #x = v * cos(ll*alpha) * cos(bb*alpha)
    #y = v * sin(ll*alpha) * cos(bb*alpha)
    #z = v * sin(bb*alpha)
    
    #Vh = float(Vlg)-(x*cosl*cosb+y*sinl*cosb+z*sinb)
    Vh=float(Vlg)+16.*cosl*cosb-315.*sinl*cosb+22.*sinb
    vls = float(Vh)-26.*cosl*cosb+317.*sinl*cosb-8.*sinb
    
    
    return vls
  
  
### The Vlg used in MKgroups is their own version. 
def Vlg2Vh(el,b, Vlg):
  
    alpha = pi / 180.
    
    cosb = cos(b*alpha)
    sinb = sin(b*alpha)
    cosl = cos(el*alpha)
    sinl = sin(el*alpha)
    
    

    Vh=float(Vlg)+16.*cosl*cosb-315.*sinl*cosb+22.*sinb
    
    
    
    return Vh
  
################################################################# 

if __name__ == '__main__':
  
  #if len(sys.argv) < 2:
    #print "\nEnter the sky patch as input ..." 
    #print "\nexample: \n > python " + str(sys.argv[0])+ " north"
    #print "\nPossible options: north, south" 
    #print "Use north/south for whole sky.\n"
    #sys.exit(1)
  
  #sky = str(sys.argv[1])
  
  VirgoW_g = np.genfromtxt('brent_virgoW.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  VirgoW_id  = VirgoW_g['PGC']
  VirgoW_d   = VirgoW_g['D']
  VirgoW_ed  = VirgoW_g['eD']  
  VirgoW_T   = VirgoW_g['T'] 
  VirgoW_Ks  = VirgoW_g['Ks'] 
  
  MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  MW_id  = MW_g['PGC']
  MW_dm  = MW_g['dm']
  MW_edm = MW_g['edm']
  
  M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  M31_id  = M31_g['PGC']  
  M31_dm  = M31_g['dm']
  M31_edm = M31_g['edm'] 
  
  LEDA	= np.genfromtxt('LEDA_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 1
  TMRS	= np.genfromtxt('2MRS_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 2
  TMPP	= np.genfromtxt('2M++_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 3
  CF3D	= np.genfromtxt('CF3D_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 4
  MKgr	= np.genfromtxt('MKgroups_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )	# 5 
  Upda	= np.genfromtxt('Updated_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )	# 6
  V8KK	= np.genfromtxt('V8K_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 7
  
  
  LEDAp	= np.genfromtxt('LEDA++_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 8
  # LEDA++_allsky.csv
  # http://leda.univ-lyon1.fr/leda/fullsql.html
  # 0<v<5000
  # btc>16
  # NED REF code: 1989PGC...C...0000P
  ### NOTE all catals should be sorted in terms of pgc
  ### NOTE all catals should be sorted in terms of pgc 
  NEWdist	= np.genfromtxt('new_distances.txt' , delimiter=',', filling_values="-100000", names=True, dtype=None )	
  pgc_newdist = NEWdist['PGC']
  dist_newdist = NEWdist['d']

  
  NorthFlag01 = np.genfromtxt('NorthFlag_0_300.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
  NorthFlag02 = np.genfromtxt('NorthFlag_300_10000.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
  SouthFlag01 = np.genfromtxt('SouthFlag_0_1500.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
  
  ##NOTE
  ## They are made by "Brent_flag_maker.py"
  ###############
  ###############
  BrentNorthFlag = np.genfromtxt('brent_north_bad_stars.v01.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
  BrentSouthFlag = np.genfromtxt('brent_south_bad_stars.v01.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
  pgcNBFlag = BrentNorthFlag['pgc']
  pgcSBFlag = BrentSouthFlag['pgc']
  pgcBFlag = np.concatenate((pgcNBFlag, pgcSBFlag))
  indices = np.argsort(pgcBFlag)
  pgcBFlag = pgcBFlag[indices]
  
  pgcNFlag01 = NorthFlag01['pgc']
  NFlag01 = NorthFlag01['flag']
  
  pgcNFlag02 = NorthFlag02['pgc']
  NFlag02 =   NorthFlag02['flag']
  
  pgcSFlag01 = SouthFlag01['pgc']
  SFlag01 =  SouthFlag01['flag']
  

  pgcFlag = np.concatenate((pgcNFlag01, pgcNFlag02, pgcSFlag01))
  StarGalFlags =  np.concatenate((NFlag01, NFlag02, SFlag01))
  indices = np.argsort(pgcFlag)
  pgcFlag   =  pgcFlag[indices]
  StarGalFlags  =    StarGalFlags[indices]
  
  #for i in range(0,20):
    #print pgcFlag[i], '  ...   ', StarGalFlags[i]
  
   
  #sys.exit("EHSAN STOPPing Point")
  
  
  pgc1 = LEDA['pgc']  # LEDA
  LEDA_sgl  = LEDA['sgl']
  LEDA_sgb  = LEDA['sgb']
  LEDA_gl   = LEDA['l2']
  LEDA_gb   = LEDA['b2']
  LEDA_b    = LEDA['btc']
  LEDA_vhelio = LEDA['v']
  LEDA_Ty = LEDA['t']
  LEDA_RA = LEDA['al2000']
  LEDA_DEC = LEDA['de2000']
  LEDA_OBJNAME = LEDA['objname']
  
  pgc1p = LEDAp['pgc']  # LEDAp
  LEDA_sglp  = LEDAp['sgl']
  LEDA_sgbp  = LEDAp['sgb']
  LEDA_glp  = LEDAp['l2']
  LEDA_gbp   = LEDAp['b2']
  LEDA_bp    = LEDAp['btc']
  LEDA_vheliop = LEDAp['v']
  LEDA_Typ = LEDAp['t']
  LEDA_RAp = LEDAp['al2000']
  LEDA_DECp = LEDAp['de2000']  
  LEDA_OBJNAMEp = LEDAp['objname']
  
  pgc2 = TMRS['pgc']  # 2MRS
  TMRS_sgl = TMRS['SGL']
  TMRS_sgb = TMRS['SGB']
  TMRS_gl = TMRS['Glon']
  TMRS_gb = TMRS['Glat']
  TMRS_Kc = TMRS['K_c']
  TMRS_vh = TMRS['Vhel']
  TMRS_RA = TMRS['RAJ']
  TMRS_DEC = TMRS['DEJ']
  
  pgc3 = TMPP['pgc']  # 2M++
  TMPP_sgl = TMPP['SGL']
  TMPP_sgb = TMPP['SGB']
  TMPP_gl = TMPP['Glon']
  TMPP_gb = TMPP['Glat']
  TMPP_Ks = TMPP['Ks']  
  TMPP_vhell = TMPP['Vhel']  
    
  
  pgc4 = CF3D['pgc']  # CF3D
  CF3D_sgl = CF3D['SGL']
  CF3D_sgb = CF3D['SGB']
  CF3D_gl = CF3D['Glon']
  CF3D_gb = CF3D['Glat']
  CF3D_Ks = CF3D['Ks']  
  CF3D_vhell = CF3D['Vhel']    
  CF3D_Ty = CF3D['Ty']
  CF3D_Dist =  CF3D['Dist']  # Cosmic flows
  CF3D_DM =  CF3D['DM']  # Cosmic flows
  CF3D_eDM =  CF3D['eDM']  # Cosmic flows
  CF3D_eD   =  (0.2*log(10.))*CF3D_eDM
  
  
  
  
  
  pgc5 = MKgr['pgc']  # MKgroups
  MKgr_Ks = MKgr['Ks']
  MKgr_Vlg = MKgr['Vlg']
  MKgr_Ty = MKgr['T']
  
  
  pgc6 = Upda['pgc']  # Updated
  Upda_Ks = Upda['Ks']
  Upda_Vh = Upda['Vh']

  pgc7 = V8KK['pgc']  # V8K
  V8K_sgl = V8KK['sgL']
  V8K_sgb = V8KK['sgB']
  V8K_gl = V8KK['l']
  V8K_gb = V8KK['b']
  V8K_absB = V8KK['absB']  
  V8K_vsgr = V8KK['V_GSR']    
  V8K_Ty = V8KK['T']
   
  
  BB8 =  np.genfromtxt('brent_added_gals.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
  BB8_PGC     = BB8['PGC']
  BB8_sgl     = BB8['sgl']
  BB8_sgb     = BB8['sgb']
  BB8_gl      = BB8['l']
  BB8_gb      = BB8['b']
  BB8_DM      = BB8['dm']
  BB8_eDM     = BB8['edm']
  BB8_Vmag    = BB8['Vmag']
  BB8_OBJNAME = BB8['Name']
  BB8_Vhel    = BB8['Vhel']
  BB8_Vls     = BB8['Vls']
  BB8_dist     = BB8['dist']
  ra_h  = BB8['ra_h']
  ra_m  = BB8['ra_m']
  ra_s  = BB8['ra_s']
  dec_d = BB8['dec_d']
  dec_m = BB8['dec_m']
  dec_s = BB8['dec_s']
  
  BB8_RA  = 15.*(ra_h+ra_m/60.+ra_s/3600.)
  BB8_DEC = np.sign(dec_d)*(np.abs(dec_d)+dec_m/60.+dec_s/3600.)
  
  indices = np.argsort(BB8_PGC)
  BB8_PGC     = BB8_PGC[indices]
  BB8_sgl     = BB8_sgl[indices]
  BB8_sgb     = BB8_sgb[indices]
  BB8_gl      = BB8_gl[indices]
  BB8_gb      = BB8_gb[indices]
  BB8_DM      = BB8_DM[indices]
  BB8_eDM     = BB8_eDM[indices]
  BB8_Vmag    = BB8_Vmag[indices]
  BB8_OBJNAME = BB8_OBJNAME[indices]
  BB8_Vhel    = BB8_Vhel[indices]
  BB8_Vls     = BB8_Vls[indices] 
  BB8_RA      = BB8_RA[indices]
  BB8_DEC     = BB8_DEC[indices]
  BB8_dist    = BB8_dist[indices]
  
  
  n1 = len(pgc1)  # LEDA
  n1p = len(pgc1p)  # LEDAp
  n2 = len(pgc2)
  n3 = len(pgc3)
  n4 = len(pgc4)
  n40 = len(pgc_newdist)
  n5 = len(pgc5)
  n6 = len(pgc6)
  n7 = len(pgc7)
  n8 = len(BB8_PGC)
  
  
  n_flag = len(pgcFlag)
  n_bflag = len(pgcBFlag) # Brent star flags
  
  pgcMax = max(max(pgc1), max(pgc1p), max(pgc2), max(pgc3), max(pgc4), max(pgc5), max(pgc6), max(pgc7), max(BB8_PGC))
  
  print "MAx PGC:   ",  pgcMax
  pgc = []
  
  SGL1 = [] # 2MRS
  SGB1 = [] # 2MRS 
  RA1 =  [] # 2MRS  
  DEC1 = [] # 2MRS 
  SGL3 = [] # CF3D
  SGB3 = [] # CF3D   
  SGL2 = [] # 2M++
  SGB2 = [] # 2M++  
  SGL0 = [] # LEDA
  SGB0 = [] # LEDA
  RA0  = [] # LEDA
  DEC0 = [] # LEDA
  RA0p = [] # LEDAp
  DEC0p= [] # LEDAp
  SGL0p= [] # LEDAp
  SGB0p= [] # LEDAp
  SGL8 = [] # V8K
  SGB8 = [] # V8K
  GL1 = [] # 2MRS
  GB1 = [] # 2MRS 
  GL3 = [] # CF3D
  GB3 = [] # CF3D   
  GL2 = [] # 2M++
  GB2 = [] # 2M++  
  GL8 = [] # V8K
  GB8 = [] # V8K 
  GL0  = []   # LEDA
  GB0  = []   # LEDA
  Ty0  = []   # LEDA
  GL0p  = []   # LEDAp
  GB0p  = []   # LEDAp
  Ty0p  = []   # LEDAp  
  Ty3  = []  # Ty3 
  Ty4  = []  # Ty4
  Ty8  = []  # V8K
  Kc  = []    # 2MRS
  Ks3 = []   # CF3D
  Ks5 = []   # Updated Nearby Galaxy Catalog
  Ks4 = []   # MKgroups
  Ks2 = []   # 2M++
  absB = []   # V8K
  B   = []   # LEDA (B-corrected magnitude)
  Bp   = []   # LEDAp (B-corrected magnitude)
  Dist =  []  # Cosmic flows
  eD   =  []    # Cosmic flows
  DM =  []  # Cosmic flows
  eDM   =  []    # Cosmic flows
  
  Vhel3 =   []      # CF3D
  Vh =   []            # 2MRS
  Vhel2 =   []      # 2M++
  Vh5 =   []          # Updated Nearby Galaxy Catalog
  Vlg4 =   []        # MKgroups
  v_helio =   []  # LEDA
  v_heliop =   []  # LEDAp  
  v_sgr =   []  # V8K
  
  bb8_sgl  = []
  bb8_sgb  = []
  bb8_gl   = []
  bb8_gb   = []
  bb8_vmag = []
  bb8_obj  = []
  bb8_vhel = []
  bb8_vls  = []
  bb8_ra   = []
  bb8_dec  = []
  
  
  OBJ = []
  SelectionFlag = []
  
  print "\nI read the data, ready for the catalog compilation ...\nPlease wait ...\n"
  i1=0
  i1p=0
  i2=0
  i3=0
  i4=0
  i40 = 0  # distance extenstion April,26,2015
  i5=0
  i6=0
  i7=0
  i8 = 0 # bb8
  
  i_flag = 0
  i_bflag = 0 
  print 'LEDA #',n1,'\n2MRS #', n2,'\n2M++ #',n3,'\nCF3D #',n4,'\nMKgroups #',n5,'\nUpdated #',n6,'\nV8K #',n7,'\nLEDAp #',n1p,'\n'
  
  sys.stdout.flush()
  print '\n'
  print ' Main Table :     ',
  #print '\b'*52,
  sys.stdout.flush()
  steps = 2
  d = steps
  
  
    

  for pg in range(1,pgcMax+10):
      
   
      
      
      pp = (100.*pg/pgcMax)
      if pp > d:
        if pp<10:
          print '\b'*4+'. '+str(int(floor(pp)))+'%',
        else:
	  print '\b'*5+'. '+str(int(floor(pp)))+'%',
        sys.stdout.flush()
        d += steps

      Parcham = 0
      sgl1  = -100000 # 2MRS
      sgb1  = -100000 # 2MRS 
      ra1   = -100000 # 2MRS 
      dec1  = -100000 # 2MRS 
      sgl3  = -100000 # CF3D
      sgb3  = -100000 # CF3D   
      sgl2  = -100000 # 2M++
      sgb2  = -100000 # 2M++  
      sgl0  = -100000 # LEDA
      sgb0  = -100000 # LEDA
      ra0   = -100000 # LEDA
      dec0  = -100000 # LEDA
      ra0p  = -100000 # LEDAp
      dec0p = -100000 # LEDAp
      sgl0p = -100000 # LEDAp
      sgb0p = -100000 # LEDAp      
      sgl8  = -100000 # V8K
      sgb8  = -100000 # V8K     
      gl1   = -100000 # 2MRS
      gb1   = -100000 # 2MRS 
      gl3   = -100000 # CF3D
      gb3   = -100000 # CF3D   
      gl2   = -100000 # 2M++
      gb2   = -100000 # 2M++  
      gl0   = -100000   # LEDA
      gb0   = -100000   # LEDA
      gl0p  = -100000   # LEDAp
      gb0p  = -100000   # LEDAp      
      gl8   = -100000   # V8K
      gb8   = -100000   # V8K      
      ty0   = -100000   # LEDA
      ty0p  = -100000   # LEDAp      
      ty3   = -100000  # CF3D 
      ty4   = -100000  # MKgroups
      ty8   = -100000  # V8K
      kc    = -100000    # 2MRS
      ks3   = -100000   # CF3D
      ks5   = -100000   # Updated Nearby Galaxy Catalog
      ks4   = -100000   # MKgroups
      ks2   = -100000   # 2M++
      b     = -100000   # LEDA (B-corrected magnitude)
      bp    = -100000   # LEDAp (B-corrected magnitude)      
      absb  = -100000  # V8K
      dist  =  -100000  # CF3D
      ed    =  -100000    # CF3D
      dm    =  -100000    # CF3D
      edm   =  -100000    # CF3D
      
      v1 =   -100000  # LEDA
      v1p =   -100000  # LEDAp    
      v2 =   -100000  # 2MRS
      v3 =   -100000  # 2M++  
      v4 =   -100000  # CF3D
      v5 =   -100000  # MKgroups
      v6 =   -100000  # Updated Nearby Galaxy Catalog
      v7 =   -100000  # V8K
      objname = '-'
      sflag = 1
      
      # BB8 varibales
      bb8_sgl_var  = -100000
      bb8_sgb_var  = -100000
      bb8_gl_var   = -100000
      bb8_gb_var   = -100000
      bb8_vmag_var = -100000
      bb8_obj_var  = '-'
      bb8_vhel_var = -100000
      bb8_vls_var  = -100000
      bb8_ra_var   = -100000
      bb8_dec_var  = -100000      
      
      if i1<n1 and pgc1[i1] == pg:  # LEDA
	Parcham+=1
	ra0  =  LEDA_RA[i1]
	dec0 =  LEDA_DEC[i1]
	sgl0 =  LEDA_sgl[i1]
	sgb0 =  LEDA_sgb[i1]
        gl0  =  LEDA_gl[i1]
        gb0  =  LEDA_gb[i1]
        ty0  =  LEDA_Ty[i1]
        b    =  LEDA_b[i1]
        v1 = LEDA_vhelio[i1]
        objname = LEDA_OBJNAME[i1]
        i1+=1
 
      if i1p<n1p and pgc1p[i1p] == pg:  # LEDAp
	Parcham+=1
	#if pg == 4075680: print "CHEEECKCKCK"
	ra0p  =  LEDA_RAp[i1p]
	dec0p =  LEDA_DECp[i1p]	
	sgl0p =  LEDA_sglp[i1p]
	sgb0p =  LEDA_sgbp[i1p]
        gl0p  =  LEDA_glp[i1p]
        gb0p  =  LEDA_gbp[i1p]
        ty0p  =  LEDA_Typ[i1p]
        bp    =  LEDA_bp[i1p]
        v1p = LEDA_vheliop[i1p]
        objname = LEDA_OBJNAMEp[i1p]
        i1p+=1 
 
      if i2<n2 and pgc2[i2] == pg:  # 2MRS
	#print pg, pgc1[i2], i2
	Parcham+=1
	ra1  = TMRS_RA[i2]
	dec1 = TMRS_DEC[i2]
	sgl1 = TMRS_sgl[i2]
	sgb1 = TMRS_sgb[i2]
        gl1  = TMRS_gl[i2]
        gb1  = TMRS_gb[i2]
        kc   = TMRS_Kc[i2]
        v2 = TMRS_vh[i2]
	i2+=1

      if i3<n3 and pgc3[i3] == pg:  # 2M++
	Parcham+=1
	sgl2 = TMPP_sgl[i3]
	sgb2 = TMPP_sgb[i3]
	gl2 =  TMPP_gl[i3]
	gb2 =  TMPP_gb[i3]
	ks2 =  TMPP_Ks[i3]
	v3 =   TMPP_vhell[i3]
	i3+=1  
  
      if i4<n4 and pgc4[i4] == pg:  # CF3D
	Parcham+=1
	sgl3 = CF3D_sgl[i4]
	sgb3 = CF3D_sgb[i4]
	gl3 = CF3D_gl[i4]
	gb3 = CF3D_gb[i4]
	ks3 = CF3D_Ks[i4]
	ty3 = CF3D_Ty[i4]
	v4 = CF3D_vhell[i4]
	dist = CF3D_Dist[i4]
	ed = CF3D_eD[i4]
	dm = CF3D_DM[i4]
	edm = CF3D_eDM[i4]
	i4+=1  
      
      
      if i40<n40 and pgc_newdist[i40] == pg and dist==0:  # distance extension
	#dist = dist_newdist[i40]
	#ed = 0.1
	i40+=1  
	
      
      
      if i5<n5 and pgc5[i5] == pg:  # MKgroups
	Parcham+=1
	ty4 = MKgr_Ty[i5]
	ks4 = MKgr_Ks[i5]
	v5 = MKgr_Vlg[i5]
	i5+=1  
	
      if i6<n6 and pgc6[i6] == pg:  # Updated
	Parcham+=1
	ks5 = Upda_Ks[i6]
	v6 = Upda_Vh[i6]
	i6+=1  	


      if i7<n7 and pgc7[i7] == pg:  # V8K
	Parcham+=1
        sgl8  = V8K_sgl[i7]  # V8K
        sgb8  = V8K_sgb[i7]  # V8K  
        gl8  = V8K_gl[i7]    # V8K
        gb8  = V8K_gb[i7]    # V8K
        ty8  = V8K_Ty[i7]    # V8K
        absb = V8K_absB[i7]  # V8K
        v7 =   V8K_vsgr[i7]  # V8K
	i7+=1  	



      if i8<n8 and BB8_PGC[i8] == pg:  # BB8
	Parcham+=1
	bb8_sgl_var  = BB8_sgl[i8]
	bb8_sgb_var  = BB8_sgb[i8]
	bb8_gl_var   = BB8_gl[i8]
	bb8_gb_var   = BB8_gb[i8]
	dm   = BB8_DM[i8]
	edm  = BB8_eDM[i8]
	dist = BB8_dist[i8]
	ed   = edm*(0.2*log(10.))
	bb8_vmag_var = BB8_Vmag[i8]
	bb8_obj_var  = BB8_OBJNAME[i8]
	bb8_vhel_var = BB8_Vhel[i8]
	bb8_vls_var  = BB8_Vls[i8]
	bb8_ra_var   = BB8_RA[i8]
	bb8_dec_var  = BB8_DEC[i8]
	i8+=1

      if i_flag<n_flag and pgcFlag[i_flag] == pg:
	sflag = StarGalFlags[i_flag]  # if sflag = 0 --> then the object is a star 
	i_flag += 1
      if i_bflag<n_bflag and pgcBFlag[i_bflag] == pg:
	sflag = 0 # Brent Stars
	i_bflag += 1

      
      if Parcham > 0 :
	
	pgc.append(pg)
	
	RA0.append(ra0)     # LEDA
	DEC0.append(dec0)     # LEDA
	SGL0.append(sgl0)     # LEDA
	SGB0.append(sgb0)     # LEDA
	GL0.append(gl0)      # LEDA
	GB0.append(gb0)      # LEDA
	Ty0.append(ty0)      # LEDA
	B.append(b)        # LEDA (B-corrected magnitude)
	v_helio.append(v1)  # LEDA

	RA0p.append(ra0p)     # LEDAp
	DEC0p.append(dec0p)     # LEDAp
	SGL0p.append(sgl0p)     # LEDAp
	SGB0p.append(sgb0p)     # LEDAp
	GL0p.append(gl0p)      # LEDAp
	GB0p.append(gb0p)      # LEDAp
	Ty0p.append(ty0p)      # LEDAp
	Bp.append(bp)        # LEDAp (B-corrected magnitude)
	v_heliop.append(v1p)  # LEDAp


	RA1.append(ra1)  # 2MRS
	DEC1.append(dec1) # 2MRS
	SGL1.append(sgl1) # 2MRS
	SGB1.append(sgb1) # 2MRS 
	GL1.append(gl1)  # 2MRS
	GB1.append(gb1)  # 2MRS 
	Kc.append(kc)   # 2MRS
	Vh.append(v2)   # 2MRS

	SGL2.append(sgl2)  # 2M++
	SGB2.append(sgb2)  # 2M++ 
	GL2.append(gl2)    # 2M++
	GB2.append(gb2)    # 2M++  
	Ks2.append(ks2)    # 2M++
	Vhel2.append(v3)   # 2M++

	SGL3.append(sgl3)  # CF3D
	SGB3.append(sgb3)  # CF3D   
	GL3.append(gl3)    # CF3D
	GB3.append(gb3)    # CF3D   
	Ks3.append(ks3)    # CF3D
	Vhel3.append(v4)   # CF3D
	Dist.append(dist)  # CF3D
	eD.append(ed)      # CF3D
	DM.append(dm)  # CF3D
	eDM.append(edm)      # CF3D	
	Ty3.append(ty3)    # CF3D

	Ks4.append(ks4)   # MKgroups
	Vlg4.append(v5)   # MKgroups
	Ty4.append(ty4)   # MKgroups

	Ks5.append(ks5)   # Updated Nearby Galaxy Catalog
	Vh5.append(v6)    # Updated Nearby Galaxy Catalog
  
        SGL8.append(sgl8) # V8K
        SGB8.append(sgb8) # V8K
        GL8.append(gl8) # V8K
        GB8.append(gb8) # V8K 
        Ty8.append(ty8)  # V8K
        absB.append(absb)   # V8K
        v_sgr.append(v7)  # V8K

        OBJ.append(objname)
        SelectionFlag.append(sflag)
        
	bb8_sgl.append(bb8_sgl_var)
	bb8_sgb.append(bb8_sgb_var)
	bb8_gl.append(bb8_gl_var)
	bb8_gb.append(bb8_gb_var)
	bb8_vmag.append(bb8_vmag_var)
	bb8_obj.append(bb8_obj_var)
	bb8_vhel.append(bb8_vhel_var)
	bb8_vls.append(bb8_vls_var)
	bb8_ra.append(bb8_ra_var)
	bb8_dec.append(bb8_dec_var)
  

  
  pgc   = np.asarray(pgc)
  N_galaxies = len(pgc)
  print "\nN_gal", N_galaxies
  
  SGL1 = np.asarray(SGL1) # 2MRS
  SGB1 = np.asarray(SGB1) # 2MRS 
  
  RA1 = np.asarray(RA1) # 2MRS
  DEC1 = np.asarray(DEC1) # 2MRS   
  
  SGL3 = np.asarray(SGL3) # CF3D
  SGB3 = np.asarray(SGB3) # CF3D   
  
  SGL2 = np.asarray(SGL2) # 2M++
  SGB2 = np.asarray(SGB2) # 2M++  

  SGL0  = np.asarray(SGL0) # LEDA
  SGB0  = np.asarray(SGB0) # LEDA
  RA0   = np.asarray(RA0)  # LEDA
  DEC0  = np.asarray(DEC0) # LEDA 
  
  SGL0p  = np.asarray(SGL0p) # LEDAp
  SGB0p  = np.asarray(SGB0p) # LEDAp  
  RA0p   = np.asarray(RA0p)  # LEDAp
  DEC0p  = np.asarray(DEC0p) # LEDAp   
  
  GL1 = np.asarray(GL1) # 2MRS
  GB1 = np.asarray(GB1) # 2MRS 
  
  GL3 = np.asarray(GL3) # CF3D
  GB3 = np.asarray(GB3) # CF3D   
  
  GL2 = np.asarray(GL2) # 2M++
  GB2 = np.asarray(GB2) # 2M++  

  GL0  = np.asarray(GL0)   # LEDA
  GB0  = np.asarray(GB0)   # LEDA
  GL0p  = np.asarray(GL0p)   # LEDAp
  GB0p  = np.asarray(GB0p)   # LEDAp
  
  Ty0  = np.asarray(Ty0)   # LEDA
  Ty0p  = np.asarray(Ty0p)   # LEDAp  
  Ty3  = np.asarray(Ty3)  # Ty3 
  Ty4  = np.asarray(Ty4)  # Ty4
  
  Kc  = np.asarray(Kc)    # 2MRS
  Ks3 = np.asarray(Ks3)   # CF3D
  Ks5 = np.asarray(Ks5)   # Updated Nearby Galaxy Catalog
  Ks4 = np.asarray(Ks4)   # MKgroups
  Ks2 = np.asarray(Ks2)   # 2M++
  B   = np.asarray(B)   # LEDA (B-corrected magnitude)
  Bp   = np.asarray(Bp)   # LEDAp (B-corrected magnitude)
  
  Dist =  np.asarray(Dist)  # Cosmic flows
  eD   =  np.asarray(eD)    # Cosmic flows
  DM =  np.asarray(DM)  # Cosmic flows
  eDM =  np.asarray(eDM)  # Cosmic flows
  
  SGL8 = np.asarray(SGL8)    # V8K
  SGB8 = np.asarray(SGB8)    # V8K
  GL8 = np.asarray(GL8)      # V8K
  GB8 = np.asarray(GB8)      # V8K
  Ty8 = np.asarray(Ty8)      # V8K
  absB = np.asarray(absB)    # V8K
  v_sgr = np.asarray(v_sgr)  # V8K
 
  
  bb8_sgl  = np.asarray(bb8_sgl)    # BB8
  bb8_sgb  = np.asarray(bb8_sgb)    # BB8
  bb8_gl   = np.asarray(bb8_gl)     # BB8
  bb8_gb   = np.asarray(bb8_gb)     # BB8
  bb8_vmag = np.asarray(bb8_vmag)   # BB8
  bb8_obj  = np.asarray(bb8_obj)    # BB8
  bb8_vhel = np.asarray(bb8_vhel)   # BB8
  bb8_vls  = np.asarray(bb8_vls)    # BB8
  bb8_ra   = np.asarray(bb8_ra)     # BB8
  bb8_dec  = np.asarray(bb8_dec)    # BB8
  
  
  
  
  KK = np.zeros((N_galaxies,), dtype=float)
  BB = np.zeros((N_galaxies,), dtype=float)
  for i in range(0, N_galaxies):
    if Kc[i]!= -100000 and Kc[i]!= 0: 
       KK[i] = Kc[i]
    elif Ks3[i]!= -100000 and Ks3[i]!= 0:
       KK[i] = Ks3[i]
    elif Ks5[i]!= -100000 and Ks5[i]!= 0:
       KK[i] = Ks5[i]
    elif Ks4[i]!= -100000 and Ks4[i]!= 0:
       KK[i] = Ks4[i]
    elif Ks2[i]!= -100000 and Ks2[i]!= 0:
       KK[i] = Ks2[i]
    else:
       KK[i] = -100000
     
    if  B[i] != -100000 and B[i]!= 0:
       BB[i] = B[i]   
    elif  Bp[i] != -100000 and Bp[i]!= 0:
       BB[i] = Bp[i]
    elif absB[i] != -100000 and absB[i]!= 0 and v_sgr[i]>0:
       BB[i] = absB[i] + 15 + 5*log10(v_sgr[i])
    else: 
       BB[i] = -100000  
  
  

  ra = np.zeros((N_galaxies,), dtype=float)
  dec = np.zeros((N_galaxies,), dtype=float)  
  sgl = np.zeros((N_galaxies,), dtype=float)
  sgb = np.zeros((N_galaxies,), dtype=float)
  gl = np.zeros((N_galaxies,), dtype=float)
  gb = np.zeros((N_galaxies,), dtype=float)  
  Ty  = np.zeros((N_galaxies,), dtype=float)
  Ks  = np.zeros((N_galaxies,), dtype=float)
  coordinate = np.zeros((N_galaxies,), dtype='a10')
  Ty_source = np.zeros((N_galaxies,), dtype='a10')
  Ks_source = np.zeros((N_galaxies,), dtype='a15')
  flag  = np.zeros((N_galaxies,), dtype=float)
  LEDAflag  = np.zeros((N_galaxies,), dtype=float)
  

  Vls =  np.zeros((N_galaxies,), dtype=float)
  Vhelio =  np.zeros((N_galaxies,), dtype=float)
  Vls_source = np.zeros((N_galaxies,), dtype='a10')
  OBJECT  = np.zeros((N_galaxies,), dtype='a35')
  
  

  print '\n'
  #print ' Starting [                                                  ]',
  #print '\b'*52,
  print ' Catalogue :     ',
  #print '\b'*52,
  sys.stdout.flush()
  steps = 2
  d = steps
  
  for i in range(0, N_galaxies):
    
     pp = (100.*i/N_galaxies)
     if pp > d:
       if pp<10:
          print '\b'*4+'. '+str(int(floor(pp)))+'%',
       else:
	  print '\b'*5+'. '+str(int(floor(pp)))+'%',
       sys.stdout.flush()
       d += steps
     
     LEDAflag[i] = 1
     V_heliocenter = -100000
     if bb8_sgl[i] != -100000:  # BB8
       sgl[i] = bb8_sgl[i]
       sgb[i] = bb8_sgb[i]
       gl[i]  = bb8_gl[i]
       gb[i]  = bb8_gb[i]
       ra[i]  = bb8_ra[i]
       dec[i] = bb8_dec[i]
       coordinate[i] = 'Brent'
       LEDAflag[i] = 0
     elif SGL0[i] != -100000:  # LEDA
       sgl[i] = SGL0[i]
       sgb[i] = SGB0[i]  
       gl[i]  = GL0[i]
       gb[i]  = GB0[i] 
       ra[i]  = 15.*RA0[i]  # hour to degree
       dec[i] = DEC0[i]
       coordinate[i] = 'LEDA'
       #LEDAflag[i] = 1
     elif SGL0p[i] != -100000:  # LEDAp
       sgl[i] = SGL0p[i]
       sgb[i] = SGB0p[i]  
       gl[i]  = GL0p[i]
       gb[i]  = GB0p[i] 
       ra[i]  = 15.*RA0p[i]  # hour to degree
       dec[i] = DEC0p[i]  
       coordinate[i] = 'LEDA++'  
       #LEDAflag[i] = 1
     elif SGL1[i] != -100000:    # 2MRS
       sgl[i] = SGL1[i]
       sgb[i] = SGB1[i]
       gl[i]  = GL1[i]
       gb[i]  = GB1[i]
       ra[i]  = RA1[i]
       dec[i] = DEC1[i]       
       coordinate[i] = '2MRS'
       LEDAflag[i] = 0
     elif SGL3[i] != -100000:  # CF3D
       sgl[i] = SGL3[i]
       sgb[i] = SGB3[i]
       gl[i]  = GL3[i]
       gb[i]  = GB3[i]
       point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
       ra[i]  = point.fk5.ra.degree
       dec[i] = point.fk5.dec.degree
       coordinate[i] = 'CF3D'
       LEDAflag[i] = 0
     elif SGL8[i] != -100000 and SGL2[i] != -100000:  # 2M++
       sgl[i] = SGL2[i]
       sgb[i] = SGB2[i]
       gl[i]  = GL2[i]
       gb[i]  = GB2[i]
       point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
       ra[i]  = point.fk5.ra.degree
       dec[i] = point.fk5.dec.degree
       coordinate[i] = '2M++'
       LEDAflag[i] = 0
  
     elif SGL8[i] != -100000:  # V8K
       sgl[i] = SGL8[i]
       sgb[i] = SGB8[i]  
       gl[i]  = GL8[i]
       gb[i]  = GB8[i] 
       point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
       ra[i]  = point.fk5.ra.degree
       dec[i] = point.fk5.dec.degree
       coordinate[i] = 'V8K'  
       LEDAflag[i] = 0
     else:
       LEDAflag[i] = 0
       flag[i] += 1
     
     
     
     if Ty0[i]!= -100000: 
       Ty[i] =  Ty0[i]   # LEDA
       Ty_source[i] = 'LEDA'
       #LEDAflag[i] = 0
     elif Ty0p[i]!= -100000: 
       Ty[i] =  Ty0p[i]   # LEDAp
       Ty_source[i] = 'LEDA++'  
       #LEDAflag[i] = 0
     elif Ty3[i]!= -100000:  
       Ty[i] =  Ty3[i]   # CF3D
       Ty_source[i] = 'CF3D'
       LEDAflag[i] = 0
     elif Ty4[i] != -100000:  
       Ty[i] =  Ty4[i]   # MKgroups
       Ty_source[i] = 'MKgroups'
       LEDAflag[i] = 0
     elif Ty8[i] != -100000: 
       Ty[i] =  Ty8[i]   # V8K
       Ty_source[i] = 'V8K'   
       LEDAflag[i] = 0
     else:
       Ty[i] = -100000
       Ty_source[i] = 'Nan'
       
    
     
     #print i, Vhel3[i], Vh[i], Vhel2[i], Vh5[i], Vlg4[i], v_helio[i]
     
     
     # Converting the heliocentric/localGroup velocity to V_ls (local sheet velocity)
     if Vhel3[i] != -100000 and Vhel3[i]!=0:
       Vls[i] = Vh2Vls(gl[i], gb[i], Vhel3[i])
       Vhelio[i] =  Vhel3[i]
       Vls_source[i] = 'CF3D'
       V_heliocenter = Vhel3[i]
       LEDAflag[i] = 0
       
     if Vhel3[i] == 0:
       Vls[i] = 0
       Vhelio[i] =  0
       V_heliocenter = 0
       Vls_source[i] = 'CF3D'
       LEDAflag[i] = 0  
     elif bb8_sgl[i] != -100000:  # BB8
       Vls[i] = bb8_vls[i]
       Vhelio[i] =  bb8_vhel[i]
       Vls_source[i] = 'Brent'
       LEDAflag[i] = 0
     elif Vh[i] != -100000:
       V_heliocenter = Vh[i]
       Vhelio[i] =  Vh[i]
       Vls[i] = Vh2Vls(gl[i], gb[i], Vh[i])
       
       Vls_source[i] = '2MRS'
       LEDAflag[i] = 0
     elif Vhel2[i] != -100000:
       V_heliocenter = Vhel2[i]
       Vhelio[i] =  Vhel2[i]
       Vls[i] = Vh2Vls(gl[i], gb[i], Vhel2[i])
       Vls_source[i] = '2M++'
       LEDAflag[i] = 0
     elif Vh5[i] != -100000:
       V_heliocenter = Vh5[i]
       Vhelio[i] =  Vh5[i]
       Vls[i] = Vh2Vls(gl[i], gb[i], Vh5[i])
       Vls_source[i] = 'Updated'
       LEDAflag[i] = 0
     elif Vlg4[i] != -100000:
       Vhelio[i] = Vlg2Vh(gl[i], gb[i], Vlg4[i])
       Vls[i] = Vlg2Vls(gl[i], gb[i], Vlg4[i])
       Vls_source[i] = 'MKgroups'
       LEDAflag[i] = 0
     elif v_helio[i] != -100000:
       V_heliocenter = v_helio[i]
       Vhelio[i] =  v_helio[i]
       Vls[i] = Vh2Vls(gl[i], gb[i], v_helio[i])
       Vls_source[i] = 'LEDA'
     elif v_heliop[i] != -100000:
       V_heliocenter = v_heliop[i] 
       Vhelio[i] =  v_heliop[i] 
       Vls[i] = Vh2Vls(gl[i], gb[i], v_heliop[i])
       Vls_source[i] = 'LEDA++'       
     elif v_sgr[i] != -100000:
       vh = vgsr_to_vhel(gl[i], gb[i], v_sgr[i])
       Vhelio[i] =  vh
       Vls[i] = Vh2Vls(gl[i], gb[i], vh)
       Vls_source[i] = 'V8K'
       LEDAflag[i] = 0
     else:
       #print pgc[i], "Error: Vls"
       Vls_source[i] = 'NaN'
       Vls[i] = 0.
       Vhelio[i] = 0
       #flag[i] += 1
       LEDAflag[i] = 2
     if pgc[i] in [166167, 2815820, 2815822, 2815823, 4689187] :
       Vhelio[i] = 0
       Vls[i] = 0
       Vls_source[i] = 'TRGB-d'
       LEDAflag[i] = 0
     
     if pgc[i] == 48515:
       Vhelio[i] = 600
       Vls[i] = 356
       Vls_source[i] = 'Brent'
     if pgc[i] == 166179:
       Vhelio[i] = 507
       Vls[i] = 264
       Vls_source[i] = 'Brent'      
       

     #### To test Vls by V8K ...
     ##if pgc[i] == 2800438 or pgc[i] == 2807070:
      ##vh_t = vgsr_to_vhel(gl[i], gb[i], v_sgr[i])
      ##Vls_t = Vh2Vls(gl[i], gb[i], vh_t)
      ##print pgc[i] , Vls_t, vh_t
     
     
     if bb8_sgl[i] != -100000:
       OBJECT[i] = bb8_obj[i]
     else:
       OBJECT[i] = OBJ[i] 
      

     if Kc[i]!= -100000 and Kc[i]!= 0: 
       Ks[i] = Kc[i]
       Ks_source[i] = '2MRS'
       LEDAflag[i] = 0
     elif Ks3[i]!= -100000 and Ks3[i]!= 0:
       Ks[i] = Ks3[i]
       Ks_source[i] = 'CF3D'
       LEDAflag[i] = 0
     elif Ks5[i]!= -100000 and Ks5[i]!= 0:
       Ks[i] = Ks5[i]
       Ks_source[i] = 'Updated'
       LEDAflag[i] = 0
     elif Ks4[i]!= -100000 and Ks4[i]!= 0:
       Ks[i] = Ks4[i]
       Ks_source[i] = 'MKgroups'
       LEDAflag[i] = 0
     elif Ks2[i]!= -100000 and Ks2[i]!= 0:
       Ks[i] = Ks2[i]
       Ks_source[i] = '2M++'
       LEDAflag[i] = 0
     elif  B[i] != -100000 and B[i]!= 0 and Ty[i] > -20:
       Ks_source[i] = 'B2Ks_LEDA'   # Converting B magnitude to Ks magnitude
       if Ty[i] < 2:
         Ks[i] = B[i]-4.10
       elif  Ty[i] >= 2 and Ty[i] <=9:
         Ks[i] = B[i]-4.60+0.25*Ty[i]
       else:
         Ks[i] = B[i]-2.35      
     
     elif  Bp[i] != -100000 and Bp[i]!= 0 and Ty[i] > -20:
       Ks_source[i] = 'B2Ks_LEDA++'   # Converting B magnitude to Ks magnitude
       if Ty[i] < 2:
         Ks[i] = Bp[i]-4.10
       elif  Ty[i] >= 2 and Ty[i] <=9:
         Ks[i] = Bp[i]-4.60+0.25*Ty[i]
       else:
         Ks[i] = Bp[i]-2.35        
     
     elif absB[i] != -100000 and absB[i]!= 0 and Ty[i] > -20 and v_sgr[i]>0:
       Ks_source[i] = 'V8K'
       LEDAflag[i] = 0
       B_app_mag = absB[i] + 15 + 5*log10(v_sgr[i])
       if Ty[i] < 2:
         Ks[i] = B_app_mag-4.10
       elif  Ty[i] >= 2 and Ty[i] <=9:
         Ks[i] = B_app_mag-4.60+0.25*Ty[i]
       else:
         Ks[i] = B_app_mag-2.35  


     elif  B[i] != -100000 and B[i]!= 0 and Ty[i] == -100000: 

       Ks_source[i] = 'B2Ks_LED?'

       q1 = np.zeros((N_galaxies,), dtype=np.int)
       q2 = np.zeros((N_galaxies,), dtype=np.int)
       q1[np.where(BB<=B[i]+0.25)] = 1
       q2[np.where(BB>=B[i]-0.25)] = 1
       qq = q1+q2
       indices = np.where(qq==2)
       all_K =[]
       for index in indices[0]: 
	 if KK[index]!= -100000:
	    all_K.append(KK[index])
	    
       if len(all_K) != 0:
          all_K = np.asarray(all_K)
          Ks[i] = np.median(all_K)
          if Ks[i]<=-20: Ks[i] = B[i]-2.5
       else:
          Ks_source[i] = 'NaN'
          Ks[i] = 0.
          flag[i] += 1
          LEDAflag[i] = 2  # remove it from the catalog	 
       
       
     elif  Bp[i] != -100000 and Bp[i]!= 0 and Ty[i] == -100000:
       Ks_source[i] = 'B2Ks_LE+?'
       q1 = np.zeros((N_galaxies,), dtype=np.int)
       q2 = np.zeros((N_galaxies,), dtype=np.int)
       q1[np.where(BB<=Bp[i]+0.25)] = 1
       q2[np.where(BB>=Bp[i]-0.25)] = 1
       qq = q1+q2
       indices = np.where(qq==2)
       all_K =[]
       for index in indices[0]: 
	 if KK[index]!= -100000:
	    all_K.append(KK[index])
       
       if len(all_K) != 0:
          all_K = np.asarray(all_K)
          Ks[i] = np.median(all_K)
          if Ks[i]<=-20: Ks[i] = Bp[i]-2.5
       else:
          Ks_source[i] = 'NaN'
	  Ks[i] = 0.
          flag[i] += 1
          LEDAflag[i] = 2  # remove it from the catalog	  
    
     elif bb8_sgl[i]!= -100000:
       Ks[i] = bb8_vmag[i] - 3.4  #  I found it by fitting a straight line
       Ks_source[i] = 'Brent.Vmag'
       LEDAflag[i] = 0	  
       

       
     elif absB[i] != -100000 and absB[i]!= 0 and Ty[i] == -100000 and v_sgr[i]>0:

       Ks_source[i] = 'V8K?'
       LEDAflag[i] = 0
       B_app_mag = absB[i] + 15 + 5*log10(v_sgr[i])
       q1 = np.zeros((N_galaxies,), dtype=np.int)
       q2 = np.zeros((N_galaxies,), dtype=np.int)
       q1[np.where(BB<=B_app_mag+0.25)] = 1
       q2[np.where(BB>=B_app_mag-0.25)] = 1
       qq = q1+q2
       indices = np.where(qq==2)
       all_K =[]
       for index in indices[0]: 
	 if KK[index]!= -100000:
	    all_K.append(KK[index])
	    
       if len(all_K) != 0:
          all_K = np.asarray(all_K)
          Ks[i] = np.median(all_K)
          if Ks[i]<=-20: Ks[i] = B_app_mag-2.5
       else:
          Ks_source[i] = 'NaN'
	  Ks[i] = 0.
          flag[i] += 1
          LEDAflag[i] = 2  # remove it from the catalog	  
	  
       
     else:
         Ks_source[i] = 'NaN'
	 Ks[i] = 0.
         flag[i] += 1
         LEDAflag[i] = 2  # remove it from the catalog
     
     if LEDAflag[i] == 1 and V_heliocenter > 900.:
       LEDAflag[i] = 0
     
     if Ks[i]<=-20:
         flag[i] += 1

     if Dist[i] == -100000:
       Dist[i] = 0
       eD[i] = 0
       DM[i] = 0
       eDM[i] = 0
       
     
     # Brent set these to zero, probably because of bad distances they have
     if pgc[i] in [2142, 70027, 46938, 43330, 14897]:
       Dist[i] = 0
       eD[i] = 0
       DM[i] = 0
       eDM[i] = 0
     
     ## from SN-type II (it's not the official distance)
     if pgc[i] in [40001]:
       Dist[i] = 16
       eD[i] = 0.20
       DM[i] = (log10(Dist[i])+5.)*5.
       eDM[i] = (5.*eD[i])/log(10.)

     if pgc[i] in [166096]:
       Dist[i] = 9.27
       eDM[i] = 0.20
       eD[i] = eDM[i]*log(10.)*0.2
       DM[i] = (log10(Dist[i])+5.)*5.

#########################################
# accurate distances for local groups MW and M31
     for j in range(len(MW_g)):
       if pgc[i] == MW_id[j]:
	  DM[i]   = MW_dm[j]
	  eDM[i]  = MW_edm[j]
	  Dist[i] = 10**(DM[i]/5.-5.)
	  eD[i]   = eDM[i]*log(10.)*0.2
	  break
  
     for j in range(len(M31_g)):
       if pgc[i] == M31_id[j]:
	  DM[i]   = M31_dm[j]
	  eDM[i]  = M31_edm[j]
	  Dist[i] = 10**(DM[i]/5.-5.)
	  eD[i]   = eDM[i]*log(10.)*0.2
	  break
  

# Milky Way Group (5064336) assumes half of M31 flux Ng=16
     
     if pgc[i] in [5064336]:
        Mk_sun  = 3.28   # K-band
	M31_d    =  0.783  # Mpc
	M31_Vls  =  -32.
	M31_Ks   =  0.797
	M31_logK = m_logK(M31_Ks, M31_Vls, d=M31_d)
	
	M31_MAG = -2.5*M31_logK+Mk_sun     # M31 Absolut Magnitude
	MW_d    = 0.008
	M       = M31_MAG + 2.5*log10(2.)  # half of the M31 flux
	MW_Ks   = M + 25 + 5. * log10(MW_d)
	Ks[i]   = MW_Ks
        Ks_source[i] = 'Brent' 

        eD[i] = 0.01
        eDM[i] = (5.*eD[i])/log(10.)
        


                
#########################################
# Virgo W (taken from brent_virgoW.csv)
     for j in range(len(VirgoW_id)):
       if pgc[i] == VirgoW_id[j]:
	  if VirgoW_d[j]>0:
	    Dist[i] = VirgoW_d[j]
            eD[i]   = VirgoW_ed[j]
            DM[i]   = (log10(Dist[i])+5.)*5.
            eDM[i]  = (5.*eD[i])/log(10.)
          if VirgoW_Ks[j]>0 : Ks[i]   = VirgoW_Ks[j]
          Ty[i]   = VirgoW_T[j]
	  break

# Virgo Foreground
     if pgc[i] == 40045: 
          Dist[i] = 9.29
          eD[i]   = 0.10
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)

     if pgc[i] == 42081: 
          Dist[i] = 9.52
          eD[i]   = 0.10
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)
          
     if pgc[i] == 43072: 
          Dist[i] = 9.63
          eD[i]   = 0.10
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)          

     if pgc[i] == 44491: 
          Dist[i] = 2.18
          eD[i]   = 0.06
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 
          
# Virgo W'
     if pgc[i] == 40375: 
          Dist[i] = 22.80
          eD[i]   = 0.12
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 

     if pgc[i] == 40122: 
          Dist[i] = 21.88
          eD[i]   = 0.12
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 
          
     if pgc[i] == 40886: 
          Dist[i] = 21.48
          eD[i]   = 0.12
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 
          
     if pgc[i] == 40033: 
          Dist[i] = 20.23
          eD[i]   = 0.15
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 
          
     if pgc[i] == 40119: 
          Dist[i] = 27.93
          eD[i]   = 0.15
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)       
          
     if pgc[i] == 41050: 
          Dist[i] = 23.01
          eD[i]   = 0.15
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)           

# Virgo W"
     if pgc[i] == 42741: 
          Dist[i] = 22.39
          eD[i]   = 0.07
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)  

     if pgc[i] == 42619: 
          Dist[i] = 21.38
          eD[i]   = 0.14
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 
          
     if pgc[i] == 43001: 
          Dist[i] = 21.38
          eD[i]   = 0.15
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)           
          
     if pgc[i] == 43254: 
          Dist[i] = 25.47
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 
          
     if pgc[i] == 42833: 
          Dist[i] = 28.44
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)           


# Virgo M

     if pgc[i] == 38890: 
          Dist[i] = 32.21
          eD[i]   = 0.12
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 

     if pgc[i] == 38916: 
          Dist[i] = 37.50
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 

     if pgc[i] == 39002: 
          Dist[i] = 34.99
          eD[i]   = 0.20
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.) 

     if pgc[i] == 39025: 
          Dist[i] = 25.82
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)           

     if pgc[i] == 39040: 
          Dist[i] = 38.19
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)           

     if pgc[i] == 39152: 
          Dist[i] = 40.55
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)             

     if pgc[i] == 39390: 
          Dist[i] = 33.42
          eD[i]   = 0.16
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)     
          
     if pgc[i] == 38749: 
          Dist[i] = 34.80
          eD[i]   = 0.15
          DM[i]   = (log10(Dist[i])+5.)*5.
          eDM[i]  = (5.*eD[i])/log(10.)           
#########################################
     # Brent name replacements 
     if pgc[i] == 9002227: pgc[i] = 3994648
     if pgc[i] == 9002244: pgc[i] = 538615
     if pgc[i] == 9002322: pgc[i] = 2801623
     
     ### 5064293 === 14241 already available in the catalog (i.e. UGCA 86)
     if pgc[i] in [2017209, 4713565, 5056949, 5056950, 5057224, 5064293]:
       flag[i] += 1
     
     if pgc[i] == 2801026:
       DM[i]   = 26.42
       eDM[i]  = 0.07 
       Dist[i] = 10**(DM[i]/5.-5.)
       eD[i]   = eDM[i]*log(10.)*0.2

     
     if pgc[i] == 2807150:
       DM[i]   = 27.08
       eDM[i]  = 0.10 
       Dist[i] = 10**(DM[i]/5.-5.)
       eD[i]   = eDM[i]*log(10.)*0.2

     
     if pgc[i] == 4689223:
       DM[i]   = 21.69
       eDM[i]  = 0.10 
       Dist[i] = 10**(DM[i]/5.-5.)
       eD[i]   = eDM[i]*log(10.)*0.2

     
     if pgc[i] == 4713553:
       DM[i]   = 19.11
       eDM[i]  = 0.08 
       Dist[i] = 10**(DM[i]/5.-5.)
       eD[i]   = eDM[i]*log(10.)*0.2

     
     if pgc[i] == 4713560:
       DM[i]   = 20.60
       eDM[i]  = 0.20 
       Dist[i] = 10**(DM[i]/5.-5.)
       eD[i]   = eDM[i]*log(10.)*0.2

     
     if pgc[i] == 4713564:
       DM[i]   = 23.10
       eDM[i]  = 0.10 
       Dist[i] = 10**(DM[i]/5.-5.)
       eD[i]   = eDM[i]*log(10.)*0.2

     if pgc[i] == 2807107:
       Vhelio[i] =  38.
       Vls[i] = Vh2Vls(gl[i], gb[i], Vhelio[i])
       Vls_source[i] = 'Modified'
       
     if pgc[i] == 9892:
       Ks[i] = 4.73
       Ks_source[i] = 'Maffei_1'       
       
     if pgc[i] == 10217:
       Ks[i] = 4.82
       Ks_source[i] = 'Maffei_2'    
 
     if pgc[i] == 100170:
       Ks[i] = 8.37
       Ks_source[i] = 'Dwing_1'    
################################################################################# 
## NAMES
################################################################################# 
     if pgc[i] == 2801026:
       OBJECT[i] = 'KKR25'     
     elif pgc[i] == 2807150:
       OBJECT[i] = 'KKH86'         
     elif pgc[i] == 4689223:
       OBJECT[i] = 'CVn I'         
     elif pgc[i] == 4713553:
       OBJECT[i] = 'BootesI'         
     elif pgc[i] == 4713564:
       OBJECT[i] = 'Leo T'         
     elif pgc[i] == 2807107:
       OBJECT[i] = 'KKH006'         
     elif pgc[i] == 4689212:
       OBJECT[i] = 'Sagittarius dSph'         
     elif pgc[i] == 5060430:
       OBJECT[i] = 'Andromeda 29'         
     elif pgc[i] == 5064336:
       OBJECT[i] = 'Milky Way'     
     elif pgc[i] == 5065056:
       OBJECT[i] = 'Andromeda 30'         
     elif pgc[i] == 5065058:
       OBJECT[i] = 'Leo P'         
     elif pgc[i] == 5065677:
       OBJECT[i] = 'Andromeda 31'         
     elif pgc[i] == 5065678:
       OBJECT[i] = 'Andromeda 32'   
     elif pgc[i] == 5067061:
       OBJECT[i] = 'Andromeda 33'         
     elif pgc[i] == 5742923:
       OBJECT[i] = 'Crater 2'         
     elif pgc[i] == 5952936:
       OBJECT[i] = 'UGC 8508'         
     elif pgc[i] == 17223:
       OBJECT[i] = 'LMC'     
     elif pgc[i] == 3085:
       OBJECT[i] = 'SMC'         
     elif pgc[i] == 4713553:
       OBJECT[i] = 'Bootes (I)'         
     elif pgc[i] == 60095:
       OBJECT[i] = 'Draco'         
     elif pgc[i] == 54074:
       OBJECT[i] = 'Ursa Minor'    
     elif pgc[i] == 3589:
       OBJECT[i] = 'Sculptor'     
     elif pgc[i] == 88608:
       OBJECT[i] = 'Sextans (I)'         
     elif pgc[i] == 19441:
       OBJECT[i] = 'Carina'         
     elif pgc[i] == 4713560:
       OBJECT[i] = 'Hercules'         
     elif pgc[i] == 10074:
       OBJECT[i] = 'Fornax'    
     elif pgc[i] == 4689223:
       OBJECT[i] = 'Canes Venatici (I)'     
     elif pgc[i] == 34176:
       OBJECT[i] = 'Leo II'         
     elif pgc[i] == 29488:
       OBJECT[i] = 'Leo I'         
     elif pgc[i] == 2557:
       OBJECT[i] = 'Andromeda'         
     elif pgc[i] == 2555:
       OBJECT[i] = 'M32'    
     elif pgc[i] == 2429:
       OBJECT[i] = 'NGC 205'     
     elif pgc[i] == 4689222:
       OBJECT[i] = 'Andromeda IX'         
     elif pgc[i] == 2666:
       OBJECT[i] = 'Andromeda I'         
     elif pgc[i] == 4608690:
       OBJECT[i] = 'Andromeda XVII'         
     elif pgc[i] == 5057230:
       OBJECT[i] = 'Andromeda XXVII' 
     elif pgc[i] == 2121:
       OBJECT[i] = 'Andromeda II'     
     elif pgc[i] == 5057228:
       OBJECT[i] = 'Andromeda XXV'         
     elif pgc[i] == 5057229:
       OBJECT[i] = 'Andromeda XXVI'         
     elif pgc[i] == 3097824:
       OBJECT[i] = 'Andromeda V'         
     elif pgc[i] == 5056923:
       OBJECT[i] = 'Andromeda XI'    
     elif pgc[i] == 5056919:
       OBJECT[i] = 'Andromeda XIX'     
     elif pgc[i] == 5057226:
       OBJECT[i] = 'Andromeda XXIII'         
     elif pgc[i] == 5056920:
       OBJECT[i] = 'Andromeda XX'         
     elif pgc[i] == 5056925:
       OBJECT[i] = 'Andromeda XIII'         
     elif pgc[i] == 5057231:
       OBJECT[i] = 'Andromeda XXI' 
     elif pgc[i] == 5056921:
       OBJECT[i] = 'Andromeda X'         
     elif pgc[i] == 5065678:
       OBJECT[i] = 'And XXXII'         
     elif pgc[i] == 2004:
       OBJECT[i] = 'NGC 147'         
     elif pgc[i] == 5065056:
       OBJECT[i] = 'And XXX'    
     elif pgc[i] == 5056922:
       OBJECT[i] = 'Andromeda XIV'     
     elif pgc[i] == 5056924:
       OBJECT[i] = 'Andromeda XII'         
     elif pgc[i] == 5056926:
       OBJECT[i] = 'Andromeda XV'         
     elif pgc[i] == 4601:
       OBJECT[i] = 'Andromeda II'         
     elif pgc[i] == 2329:
       OBJECT[i] = 'NGC 185' 
     elif pgc[i] == 5060430:
       OBJECT[i] = 'Andromeda XXIX'     
     elif pgc[i] == 5818:
       OBJECT[i] = 'Triangulum'         
     elif pgc[i] == 5057227:
       OBJECT[i] = 'Andromeda XXIV'         
     elif pgc[i] == 2807155:
       OBJECT[i] = 'Andromeda VII'         
     elif pgc[i] == 1305:
       OBJECT[i] = 'IC 10'    
     elif pgc[i] == 5065677:
       OBJECT[i] = 'And XXXI'     
     elif pgc[i] == 3792:
       OBJECT[i] = 'LGS 3'         
     elif pgc[i] == 2807158:
       OBJECT[i] = 'Andromeda VI'         
     elif pgc[i] == 5057232:
       OBJECT[i] = 'Andromeda XXII'         
     elif pgc[i] == 5056927:
       OBJECT[i] = 'Andromeda XVI'  
     elif pgc[i] == 143:
       OBJECT[i] = 'WLM'         
     elif pgc[i] == 3844:
       OBJECT[i] = 'IC 1613'         
     elif pgc[i] == 6830:
       OBJECT[i] = 'Phoenix'         
     elif pgc[i] == 26142:
       OBJECT[i] = 'UGC 4879'    
     elif pgc[i] == 28868:
       OBJECT[i] = 'Leo A'     
     elif pgc[i] == 63616:
       OBJECT[i] = 'NGC 6822'         
     elif pgc[i] == 63287:
       OBJECT[i] = 'SagDIG'         
     elif pgc[i] == 65367:
       OBJECT[i] = 'Aquarius'         
     elif pgc[i] == 69519:
       OBJECT[i] = 'Tucana' 
     elif pgc[i] == 71538:
       OBJECT[i] = 'Pegasus dIrr'     
     elif pgc[i] == 3097691:
       OBJECT[i] = 'Cetus'         
     elif pgc[i] == 4713564:
       OBJECT[i] = 'Leo T'         
     elif pgc[i] == 5056918:
       OBJECT[i] = 'Andromeda 18'         
     elif pgc[i] == 5060429:
       OBJECT[i] = 'Andromeda 28'    
     elif pgc[i] == 5067061:
       OBJECT[i] = 'Andromeda 33'     
     elif pgc[i] == 9892:
       OBJECT[i] = 'Maffei 1'         
     elif pgc[i] == 10217:
       OBJECT[i] = 'Maffei 2'         
     elif pgc[i] == 100170:
       OBJECT[i] = 'Dwingeloo 1'         
     elif pgc[i] == 101304:
       OBJECT[i] = 'Dwingeloo 2'        
     elif pgc[i] == 166068:
       OBJECT[i] = 'MB1'         
     elif pgc[i] == 166069:
       OBJECT[i] = 'KK 22'         
     elif pgc[i] == 168300:
       OBJECT[i] = 'KKH 11' 
     elif pgc[i] == 2807107:
       OBJECT[i] = 'KKH 12'     
     elif pgc[i] == 2807102:
       OBJECT[i] = 'KKH 5'         
     elif pgc[i] == 2807103:
       OBJECT[i] = 'KKH 6'         
        
################################################################################# 
##### FLAG HANDLING .... ##### 
#################################################################################  
     #sky = "south"
#################################################################################  
     if Vls[i] > 4000: flag[i] += 1
     #if sky == 'south' and gb[i] > 0: flag[i] += 1    # north < 0  or south > 0
     #if sky == 'north' and gb[i] < 0: flag[i] += 1    # north < 0  or south > 0
     if LEDAflag[i] > 0: flag[i] += 1
     if SelectionFlag[i] == 0 : flag[i] += 1

#################################################################################  
#################################################################################  
#################################################################################  
  # Select those galaxies for which the flag = 0
  
  flag = np.where(flag==0)

  if len(flag[0]) != 0 :

      pgc_file = pgc[flag]
      gl_file = gl[flag]
      gb_file = gb[flag]
      sgl_file = sgl[flag]
      sgb_file = sgb[flag]
      ra_file = ra[flag]
      dec_file = dec[flag]
      coordinate_file = coordinate[flag]
      Ty_file = Ty[flag]
      Ty_source_file = Ty_source[flag]
      Ks_file = Ks[flag]
      Ks_source_file = Ks_source[flag]
      Vls_file = Vls[flag]
      Vhelio_file = Vhelio[flag]
      Vls_source_file = Vls_source[flag]
      Dist_file = Dist[flag]
      eD_file = eD[flag]
      DM_file = DM[flag]
      eDM_file = eDM[flag]
      OBJECT_file = OBJECT[flag]
      BB_file = BB[flag]
      
      
      myTable = Table()
      myTable.add_column(Column(data=pgc_file, name='pgc'))
      myTable.add_column(Column(data=ra_file,  name='ra'))  
      myTable.add_column(Column(data=dec_file,  name='dec'))  
      myTable.add_column(Column(data=gl_file,  name='gl'))
      myTable.add_column(Column(data=gb_file,  name='gb'))
      myTable.add_column(Column(data=sgl_file, name='sgl'))
      myTable.add_column(Column(data=sgb_file, name='sgb'))
      myTable.add_column(Column(data=coordinate_file, name='coordinate_src'))
      myTable.add_column(Column(data=Ty_file,name='Ty'))
      myTable.add_column(Column(data=Ty_source_file,name='Ty_src'))
      myTable.add_column(Column(data=BB_file,name='B_mag'))
      myTable.add_column(Column(data=Ks_file,name='Ks'))
      myTable.add_column(Column(data=Ks_source_file,name='Ks_src'))
      myTable.add_column(Column(data=Vls_file,name='Vls'))
      myTable.add_column(Column(data=Vhelio_file,name='Vhelio'))
      myTable.add_column(Column(data=Vls_source_file,name='Vls_src'))
      myTable.add_column(Column(data=Dist_file,name='dcf2'))
      myTable.add_column(Column(data=eD_file,name='ed'))
      myTable.add_column(Column(data=DM_file,name='DM'))
      myTable.add_column(Column(data=eDM_file,name='eDM'))
      myTable.add_column(Column(data=OBJECT_file,name='objname'))
      
      
#      pgc | flag |       ra |      dec |       gl |       gb |      sgl |      sgb |        Ty |      B_mag |    Ks |    logK |  Vls | Vhelio |  dcf2 |   ed | mDist | mDistErr | R_theta | sigmaP_dyn | sigmaP_lum |   Mv_dyn |   Mv_lum | Rg_angular | Rg_dyn | R2t_dyn | R2t_lum |   tX_dyn |   tX_lum | No_Galaxies |    nest |  coordinate_src |     Ty_src |        Ks_src |    Vls_src |                         objname
#    10180 |    0 |  40.3146 |  38.7433 | 145.3293 | -19.3053 | 341.4192 |  -9.5717 |      -5.0 |      15.65 | 11.55 |  8.8816 |  933 |    740 |  0.00 | 0.00 |  0.00 |     0.00 | 0.00000 |        0.0 |            |          |          |            |        |         |         |          |          |           1 |   10180 |            LEDA |      LEDA  |    B2Ks_LEDA  |      LEDA  |                        UGC02165
      if True:
	pgc_new = 10180
	gl_new = 145.3293
	gb_new = -19.3053
	sgl_new = 341.4192
	sgb_new = -9.5717
	ra_new = 40.3146
	dec_new = 38.7433
	coordinate_new = 'LEDA'
	Ty_new = -5.0
	Ty_source_new = 'LEDA'
	Ks_new = 11.55
	Ks_source_new = 'BRENT'
	Vls_new = 933
	Vhelio_new = 740
	Vls_source_new = 'LEDA'
	Dist_new = 0.0
	eD_new = 0.0
	DM_new = 0.0
	eDM_new = 0.0
	OBJECT_new = 'UGC02165'
	BB_new = 15.65
	myTable.add_row([pgc_new, ra_new, dec_new, gl_new, gb_new, sgl_new, sgb_new, coordinate_new, \
	  Ty_new, Ty_source_new, BB_new, Ks_new, Ks_source_new, Vls_new, Vhelio_new, Vls_source_new, Dist_new, eD_new, DM_new, eDM_new, OBJECT_new])





      myTable.write('AllSky.v20.cf3.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
      
#################################################################################  
#################################################################################  
