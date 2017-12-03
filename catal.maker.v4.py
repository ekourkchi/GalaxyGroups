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
  
  if len(sys.argv) < 2:
    print "\nEnter the sky patch as input ..." 
    print "\nexample: \n > python " + str(sys.argv[0])+ " north"
    print "\nPossible options: north, south" 
    print "Use north/south for whole sky.\n"
    sys.exit(1)
  
  sky = str(sys.argv[1])
  
  
  LEDA	= np.genfromtxt('LEDA_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 1
  TMRS	= np.genfromtxt('2MRS_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 2
  TMPP	= np.genfromtxt('2M++_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )		# 3
  CF2D	= np.genfromtxt('CF2D2.1_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )	# 4
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
    
  
  pgc4 = CF2D['pgc']  # CF2D
  CF2D_sgl = CF2D['SGL']
  CF2D_sgb = CF2D['SGB']
  CF2D_gl = CF2D['Glon']
  CF2D_gb = CF2D['Glat']
  CF2D_Ks = CF2D['Ks']  
  CF2D_vhell = CF2D['Vhel']    
  CF2D_Dist =  CF2D['Dist']  # Cosmic flows
  CF2D_eD   =  CF2D['eD']    # Cosmic flows
  CF2D_Ty = CF2D['Ty']
  
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
   
  
  
  n1 = len(pgc1)  # LEDA
  n1p = len(pgc1p)  # LEDAp
  n2 = len(pgc2)
  n3 = len(pgc3)
  n4 = len(pgc4)
  n40 = len(pgc_newdist)
  n5 = len(pgc5)
  n6 = len(pgc6)
  n7 = len(pgc7)
  
  n_flag = len(pgcFlag)
  n_bflag = len(pgcBFlag) # Brent star flags
  
  pgcMax = max(max(pgc1), max(pgc1p), max(pgc2), max(pgc3), max(pgc4), max(pgc5), max(pgc6), max(pgc7))
  
  print "MAx PGC:   ",  pgcMax
  pgc = []
  
  SGL1 = [] # 2MRS
  SGB1 = [] # 2MRS 
  RA1 =  [] # 2MRS  
  DEC1 = [] # 2MRS 
  SGL3 = [] # CF2D
  SGB3 = [] # CF2D   
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
  GL3 = [] # CF2D
  GB3 = [] # CF2D   
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
  Ks3 = []   # CF2D
  Ks5 = []   # Updated Nearby Galaxy Catalog
  Ks4 = []   # MKgroups
  Ks2 = []   # 2M++
  absB = []   # V8K
  B   = []   # LEDA (B-corrected magnitude)
  Bp   = []   # LEDAp (B-corrected magnitude)
  Dist =  []  # Cosmic flows
  eD   =  []    # Cosmic flows
  
  Vhel3 =   []      # CF2D
  Vh =   []            # 2MRS
  Vhel2 =   []      # 2M++
  Vh5 =   []          # Updated Nearby Galaxy Catalog
  Vlg4 =   []        # MKgroups
  v_helio =   []  # LEDA
  v_heliop =   []  # LEDAp  
  v_sgr =   []  # V8K
  
  OBJ =[]
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
  
  i_flag = 0
  i_bflag = 0 
  print 'LEDA #',n1,'\n2MRS #', n2,'\n2M++ #',n3,'\nCF2D #',n4,'\nMKgroups #',n5,'\nUpdated #',n6,'\nV8K #',n7,'\nLEDAp #',n1p,'\n'
  
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
      sgl3  = -100000 # CF2D
      sgb3  = -100000 # CF2D   
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
      gl3   = -100000 # CF2D
      gb3   = -100000 # CF2D   
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
      ty3   = -100000  # CF2D 
      ty4   = -100000  # MKgroups
      ty8   = -100000  # V8K
      kc    = -100000    # 2MRS
      ks3   = -100000   # CF2D
      ks5   = -100000   # Updated Nearby Galaxy Catalog
      ks4   = -100000   # MKgroups
      ks2   = -100000   # 2M++
      b     = -100000   # LEDA (B-corrected magnitude)
      bp    = -100000   # LEDAp (B-corrected magnitude)      
      absb  = -100000  # V8K
      dist  =  -100000  # CF2D
      ed    =  -100000    # CF2D
      
      v1 =   -100000  # LEDA
      v1p =   -100000  # LEDAp    
      v2 =   -100000  # 2MRS
      v3 =   -100000  # 2M++  
      v4 =   -100000  # CF2D
      v5 =   -100000  # MKgroups
      v6 =   -100000  # Updated Nearby Galaxy Catalog
      v7 =   -100000  # V8K
      objname = '-'
      sflag = 1
      
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
  
      if i4<n4 and pgc4[i4] == pg:  # CF2D
	Parcham+=1
	sgl3 = CF2D_sgl[i4]
	sgb3 = CF2D_sgb[i4]
	gl3 = CF2D_gl[i4]
	gb3 = CF2D_gb[i4]
	ks3 = CF2D_Ks[i4]
	ty3 = CF2D_Ty[i4]
	v4 = CF2D_vhell[i4]
	dist = CF2D_Dist[i4]
	ed = CF2D_eD[i4]
	i4+=1  
      
      
      if i40<n40 and pgc_newdist[i40] == pg and dist==0:  # distance extension
	dist = dist_newdist[i40]
	ed = 0.1
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

	SGL3.append(sgl3)  # CF2D
	SGB3.append(sgb3)  # CF2D   
	GL3.append(gl3)    # CF2D
	GB3.append(gb3)    # CF2D   
	Ks3.append(ks3)    # CF2D
	Vhel3.append(v4)   # CF2D
	Dist.append(dist)  # CF2D
	eD.append(ed)      # CF2D
	Ty3.append(ty3)    # CF2D

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
  

  
  pgc   = np.asarray(pgc)
  N_galaxies = len(pgc)
  print "\nN_gal", N_galaxies
  
  SGL1 = np.asarray(SGL1) # 2MRS
  SGB1 = np.asarray(SGB1) # 2MRS 
  
  RA1 = np.asarray(RA1) # 2MRS
  DEC1 = np.asarray(DEC1) # 2MRS   
  
  SGL3 = np.asarray(SGL3) # CF2D
  SGB3 = np.asarray(SGB3) # CF2D   
  
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
  
  GL3 = np.asarray(GL3) # CF2D
  GB3 = np.asarray(GB3) # CF2D   
  
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
  Ks3 = np.asarray(Ks3)   # CF2D
  Ks5 = np.asarray(Ks5)   # Updated Nearby Galaxy Catalog
  Ks4 = np.asarray(Ks4)   # MKgroups
  Ks2 = np.asarray(Ks2)   # 2M++
  B   = np.asarray(B)   # LEDA (B-corrected magnitude)
  Bp   = np.asarray(Bp)   # LEDAp (B-corrected magnitude)
  
  Dist =  np.asarray(Dist)  # Cosmic flows
  eD   =  np.asarray(eD)    # Cosmic flows
  
  SGL8 = np.asarray(SGL8)    # V8K
  SGB8 = np.asarray(SGB8)    # V8K
  GL8 = np.asarray(GL8)      # V8K
  GB8 = np.asarray(GB8)      # V8K
  Ty8 = np.asarray(Ty8)      # V8K
  absB = np.asarray(absB)    # V8K
  v_sgr = np.asarray(v_sgr)  # V8K
 
  
  
  
  
  
  
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
     if SGL0[i] != -100000:  # LEDA
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
     elif SGL3[i] != -100000:  # CF2D
       sgl[i] = SGL3[i]
       sgb[i] = SGB3[i]
       gl[i]  = GL3[i]
       gb[i]  = GB3[i]
       point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
       ra[i]  = point.fk5.ra.degree
       dec[i] = point.fk5.dec.degree
       coordinate[i] = 'CF2D'
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
       Ty[i] =  Ty3[i]   # CF2D
       Ty_source[i] = 'CF2D'
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
       Vls_source[i] = 'CF2D'
       V_heliocenter = Vhel3[i]
       LEDAflag[i] = 0
     if Vhel3[i] == 0:
       Vls[i] = 0
       Vhelio[i] =  0
       V_heliocenter = 0
       Vls_source[i] = 'CF2D'
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
     elif pgc[i] in [166167, 2815820, 2815822, 2815823, 4689187] :
       Vhelio[i] = 0
       Vls[i] = 0
       Vls_source[i] = 'TRGB-d'
       LEDAflag[i] = 0
     else:
       #print pgc[i], "Error: Vls"
       Vls_source[i] = 'NaN'
       Vls[i] = 0.
       Vhelio[i] = 0
       #flag[i] += 1
       LEDAflag[i] = 2
     
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
     
     
     
     OBJECT[i] = OBJ[i] 

     if Kc[i]!= -100000 and Kc[i]!= 0: 
       Ks[i] = Kc[i]
       Ks_source[i] = '2MRS'
       LEDAflag[i] = 0
     elif Ks3[i]!= -100000 and Ks3[i]!= 0:
       Ks[i] = Ks3[i]
       Ks_source[i] = 'CF2D'
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
     
     # Brent set these to zero, probably because of bad distances they have
     if pgc[i] in [2142, 70027, 46938, 43330]:
       Dist[i] = 0
       eD[i] = 0
     
     # from SN-type II (it's not the official distance)
     if pgc[i] in [40001]:
       Dist[i] = 16
       eD[i] = 0.20
     
     if pgc[i] in [14897]:  # bad distance info
       Dist[i] = 0.
       eD[i] = 0.
          
     
     # Brent name replacements 
     # 9002227 -> 3994648
     # 9002244 ->  538615
     # 9002322 -> 2801623
     if pgc[i] == 9002227: pgc[i] = 3994648
     if pgc[i] == 9002244: pgc[i] = 538615
     if pgc[i] == 9002322: pgc[i] = 2801623
     
     
     
     #if  pgc[i]==166133: print "\n   checkpoint166133:    ", Dist[i], eD[i], flag[i], LEDAflag[i], SelectionFlag[i], gb[i]
     #if  pgc[i]==4326021: print "\n   checkpoint4326021:    ", Dist[i], eD[i], flag[i], LEDAflag[i],  SelectionFlag[i], gb[i]
################################################################################# 
##### FLAG HANDLING .... ##### 
#################################################################################  
     #sky = "south"
#################################################################################  
     if Vls[i] > 4000: flag[i] += 1
     if sky == 'south' and gb[i] > 0: flag[i] += 1    # north < 0  or south > 0
     if sky == 'north' and gb[i] < 0: flag[i] += 1    # north < 0  or south > 0
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
      myTable.add_column(Column(data=OBJECT_file,name='objname'))
      
      
#      pgc | flag |       ra |      dec |       gl |       gb |      sgl |      sgb |        Ty |      B_mag |    Ks |    logK |  Vls | Vhelio |  dcf2 |   ed | mDist | mDistErr | R_theta | sigmaP_dyn | sigmaP_lum |   Mv_dyn |   Mv_lum | Rg_angular | Rg_dyn | R2t_dyn | R2t_lum |   tX_dyn |   tX_lum | No_Galaxies |    nest |  coordinate_src |     Ty_src |        Ks_src |    Vls_src |                         objname
#    10180 |    0 |  40.3146 |  38.7433 | 145.3293 | -19.3053 | 341.4192 |  -9.5717 |      -5.0 |      15.65 | 11.55 |  8.8816 |  933 |    740 |  0.00 | 0.00 |  0.00 |     0.00 | 0.00000 |        0.0 |            |          |          |            |        |         |         |          |          |           1 |   10180 |            LEDA |      LEDA  |    B2Ks_LEDA  |      LEDA  |                        UGC02165
      if sky == 'south':
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
	OBJECT_new = 'UGC02165'
	BB_new = 15.65
	myTable.add_row([pgc_new, ra_new, dec_new, gl_new, gb_new, sgl_new, sgb_new, coordinate_new, \
	  Ty_new, Ty_source_new, BB_new, Ks_new, Ks_source_new, Vls_new, Vhelio_new, Vls_source_new, Dist_new, eD_new, OBJECT_new])





      myTable.write('AllSky.'+sky+'.v17.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
      
#################################################################################  
#################################################################################  
