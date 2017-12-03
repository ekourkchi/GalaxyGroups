from time import time
import sys
import os
import numpy as np
from math import *

# spherematch related libraries
from astrometry.util.starutil_numpy import * 
import spherematch as spherematch


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

VCIRC = 220. # u.km/u.s
VLSR = [10., 5.25, 7.17] # *u.km/u.s
#################################################################

def angleto3Dradius(angle, isDegree=True):
  
  if isDegree:
    angle = angle*pi/180.
  
  
  return sqrt((sin(angle))**2 + (1-cos(angle))**2)

#################################################################


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
#################################################################

### export OMP_NUM_THREADS=2

if __name__ == '__main__':
  
   time0 = time()

   
   a = pyfits.open('LEDA.head.long.fits')
   d = a[1].data
   pgc = d.field('pgc')
   #gl_ = d.field('sgl')
   #gb_ = d.field('sgb')
   gl_ = d.field('l2')
   gb_ = d.field('b2')
   logd25 = d.field('logd25')
   logr25 = d.field('logr25')
   v = d.field('v')
   vgsr = d.field('vgsr')
   
   

   
   

   
   
   EDD = np.genfromtxt('2MASS_K_11.25.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
   EDD_pgc = EDD['PGC']  # EDD
   EDD_gl_ = EDD['Glon']
   EDD_gb_ = EDD['Glat'] 
   EDD_Vhel = EDD['Vhel']   
   
   
   indices = np.where(EDD_pgc >= 9000000)
   print "all >9xxxx:", len(indices[0])
   EDD_pgc = EDD_pgc[indices]
   EDD_gl_ = EDD_gl_[indices]
   EDD_gb_ = EDD_gb_[indices]
   
   EDD_Vhel = EDD_Vhel[indices]

   

   LEDA	= np.genfromtxt('/home/ehsan/GalGroups/start/AllSky/LEDA_allsky.csv' , delimiter=',', filling_values="-100000", names=True, dtype=None )
   pgc_p = LEDA['pgc']   

   
   print "preparation time Ehsan = ", time()-time0
   print "... GO match ..."
   time0 = time()
   

   
   (m1, m2, d)  = spherematch.match_radec(EDD_gl_, EDD_gb_, gl_, gb_, 200./3600, nearest=False)
   
   p1 = 0.
   M1 = []
   M2 = []
   Del = []
   PGC1 = []
   PGC2 = []
   Angle = []
   
   for i in range(len(d)):
     d25 = 0.1*(10**logd25[m2[i]])

     
     delta = abs(v[m2[i]]- EDD_Vhel[m1[i]])
     if d[i]*3600<1 or (delta<100 and v[m2[i]]>0): #  pgc[m2[i]] == EDD_pgc[m1[i]]
       p1 += 1
       M1.append(m1[i])
       M2.append(m2[i])
       Del.append(delta)
       PGC1.append(EDD_pgc[m1[i]])
       PGC2.append(pgc[m2[i]])
       Angle.append(d[i]*3600)
       
       
       #print pgc[m2[i]], EDD_pgc[m1[i]], d[i]*3600, ' <> ', v[m2[i]], EDD_Vhel[m1[i]]
       #print pgc[m2[i]], EDD_pgc[m1[i]], d[i]*3600, ' <> ', v[m2[i]], vh, ' <> ', d25*60
   
   
   
   M1 = np.asarray(M1)
   M2 = np.asarray(M2)
   Del = np.asarray(Del)
   PGC1 = np.asarray(PGC1)
   PGC2 = np.asarray(PGC2)
   Angle = np.asarray(Angle)
   
   indices = np.argsort(PGC1)
   PGC1 = PGC1[indices]
   PGC2 = PGC2[indices]
   Del  = Del[indices]
   M1   = M1[indices]
   M2   = M2[indices]
   Angle = Angle[indices]
   
   old_id = -1
   del_list = []
   pgc1_list = []
   pgc2_list = []
   m1_list = []
   m2_list = []
   angle_list = []
   
   
   
   del_list_f = []
   pgc1_list_f = []
   pgc2_list_f = []
   m1_list_f = []
   m2_list_f = []   
   angle_list_f = []
   
   for i in range(len(PGC1)):
     
     if PGC1[i] != old_id:
       
       if len(del_list) != 0:
	 
	 indices = np.argsort(angle_list)
	 index = indices[0]
	 

	 
	 del_list_f.append(del_list[index])
	 pgc1_list_f.append(pgc1_list[index])
	 pgc2_list_f.append(pgc2_list[index])
	 m1_list_f.append(m1_list[index])
	 m2_list_f.append(m2_list[index])
	 angle_list_f.append(angle_list[index])
	 
         #if len(del_list) > 1:
	   #print len(del_list)
	   #print del_list
	   #print pgc1_list
	   #print pgc2_list
	   #print angle_list
	   #print pgc1_list[index], pgc2_list[index]
	   #print
       
       del_list = []
       pgc1_list = []
       pgc2_list = []
       m1_list = []
       m2_list = []
       angle_list = []
       del_list.append(Del[i])
       pgc1_list.append(PGC1[i])
       pgc2_list.append(PGC2[i])
       m1_list.append(M1[i])
       m2_list.append(M2[i])
       angle_list.append(Angle[i])
       old_id = PGC1[i]
     else:
       
       del_list.append(Del[i])
       pgc1_list.append(PGC1[i])
       pgc2_list.append(PGC2[i])
       m1_list.append(M1[i])
       m2_list.append(M2[i])    
       angle_list.append(Angle[i])
       
       
   myTable = Table()
   empty = []
   myTable.add_column(Column(data=empty,name='flag', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='pgc_2MASS', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='Glon_2MASS', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='Glat_2MASS', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='Vhel_2MASS', format='%0.1f'))
   myTable.add_column(Column(data=empty,name='pgc_LEDA', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='l2_LEDA', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='b2_LEDA', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='Vhel_LEDA', format='%0.1f'))
   myTable.add_column(Column(data=empty,name='d25_arcsec_LEDA', format='%0.1f'))
   myTable.add_column(Column(data=empty,name='angle_arcsec', format='%0.2f'))
   
   t = 0
   for i in range(len(m1_list_f)):
          
       flag = 0 
       for id in pgc_p:
	 if pgc2_list_f[i] == id:
	   flag = 1
	   break   
       myTable.add_row([flag, pgc1_list_f[i], EDD_gl_[m1_list_f[i]], EDD_gb_[m1_list_f[i]], EDD_Vhel[m1_list_f[i]], pgc2_list_f[i], gl_[m2_list_f[i]], gb_[m2_list_f[i]], v[m2_list_f[i]], 6.*(10**logd25[m2_list_f[i]]), angle_list_f[i]])
       t +=1

     
   myTable.write('2MASSK11.25_LEDA_match.txt', format='ascii.fixed_width',delimiter='|', bookend=False)
   print myTable














