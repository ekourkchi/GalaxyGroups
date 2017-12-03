#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2016"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"
# As of April, 8, 2016
###
# Written by Ehsan Kourkchi (September 2015)
# email: ehsan@ifa.hawaii.edu
# This code, identifies groups of galaxies, given a 
# a cataloge of galaxies. The radial velocity of galaxies would 
# be used to find galaxies with relatively the same radial velocities
# and almost the same position on the sky
###

## Importing Important Python Libraries
import sys
import os
import random
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
import subprocess



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

# Returns mass in solar unit, given the absolute K-band luminosity
def Mass(L_k):

  L = L_k / 1.E10
  
  if L < 1.:
    MtoL = 43.0#*(L**-0.3)
  elif L > 1000.:
    MtoL = 121.19
  else:
    MtoL = 43*(L**0.15)
  
  
  Mass_out = h * L_k * MtoL
  
  return Mass_out


# **************************************
# returns angular separation of 
# two vectors in radian
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
   #X2 = sqrt(x1**2 + y1**2 + z1**2)
   #Y2 = sqrt(x2**2 + y2**2 + z2**2)
   
   if XdotY > 1 :
     theta12 = 0.
   elif XdotY < -1 :
     theta12 = -1.*pi
   else:
     theta12 = acos(XdotY)  
   return theta12   # radian
# **************************************
# returns sign of a number 
def sign(x):
  if x<0 : return -1
  if x>=0 : return 1
# **************************************
# L is luminosity
# l,b are galaxtic coordinates, and d is distance
# returns the barycentric coordiantes of a pair 
def barycenter(L1, l1, b1, d1, L2, l2, b2, d2):
  
   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)
   
   x1 = d1 * cl1 * cb1
   y1 = d1 * sl1 * cb1
   z1 = d1 * sb1
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = d2 * cl2 * cb2
   y2 = d2 * sl2 * cb2
   z2 = d2 * sb2   


   L_tot = L1 + L2
   x = (L1*x1+L2*x2)/L_tot
   y = (L1*y1+L2*y2)/L_tot
   z = (L1*z1+L2*z2)/L_tot


   r = sqrt(x**2+y**2)
   b = atan(z/r)*180/pi
   l = atan(y/x)*180/pi
   

   if sign(x) < 0 : l+=180
   if sign(x) > 0 and sign(y) < 0: l+=360   
   
   if d1==0 or d2==0:
      d = 0.
   else:
      d = sqrt(x**2+y**2+z**2)


   
   return l, b, d

# **************************************
def touch(gr1, gr2, coeff=1.0):
   
   l1  = gr1.sgl
   b1  = gr1.sgb
   d1    = gr1.dist
   r1t_1 = gr1.r1t
   
   l2  = gr2.sgl
   b2  = gr2.sgb
   d2    = gr2.dist
   r1t_2 = gr2.r1t
   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)   
   
   x1 = d1 * cl1 * cb1
   y1 = d1 * sl1 * cb1
   z1 = d1 * sb1 
   
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = d2 * cl2 * cb2
   y2 = d2 * sl2 * cb2
   z2 = d2 * sb2    
   
   R12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
   
   if R12 < coeff*(r1t_1 + r1t_2):
     return True
   else: 
     return False
   

# **************************************
def readgrouplist(filename):
  
   
   try:
      mytable = np.genfromtxt(filename , delimiter='|', filling_values="-100000", names=True, dtype=None )
   
   except:
      print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
      print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
      exit(1)
  
   id          = mytable['pgc']
   flag        = mytable['flag']
   sgl         = mytable['sgl']
   sgb         = mytable['sgb']
   gl          = mytable['gl']
   gb          = mytable['gb']  
   ra          = mytable['ra']
   dec         = mytable['dec']
   Ks          = mytable['Ks']
   Vls         = mytable['Vls']
   nest        = mytable['nest']
   dcf2        = mytable['dcf2']
   ed          = mytable['ed']
   mDist       = mytable['mDist']
   mDistErr    = mytable['mDistErr']
   No_Galaxies = mytable['No_Galaxies']

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
   
   GroupList = []

   
   while gr == 2:
     
     id_grp       = [id[i]]
     flag_grp     = [flag[i]]
     sgl_grp      = [sgl[i]]
     sgb_grp      = [sgb[i]]
     gl_grp       = [gl[i]]
     gb_grp       = [gb[i]]
     ra_grp       = [ra[i]]
     dec_grp      = [dec[i]]
     Ks_grp       = [Ks[i]]
     Vls_grp      = [Vls[i]]
     nest_grp     = [nest[i]]
     dcf2_grp     = [dcf2[i]]
     ed_grp       = [ed[i]]
     mDist_grp    = [mDist[i]]
     mDistErr_grp = [mDistErr[i]]
     No_Galaxies_grp = [No_Galaxies[i]]
     
     i += 1
     while i<N and flag[i]==1: 
         
         id_grp.append(id[i])
         flag_grp.append(flag[i])
         sgl_grp.append(sgl[i])
         sgb_grp.append(sgb[i])
         gl_grp.append(gl[i])
         gb_grp.append(gb[i])
         ra_grp.append(ra[i])
         dec_grp.append(dec[i])
         Ks_grp.append(Ks[i])
         Vls_grp.append(Vls[i])
         nest_grp.append(nest[i])
         dcf2_grp.append(dcf2[i])
         ed_grp.append(ed[i])
         mDist_grp.append(mDist[i])
         mDistErr_grp.append(mDistErr[i])
         No_Galaxies_grp.append(No_Galaxies[i])
         
         i+=1
     
     mygroup = {'pgc':id_grp, 'flag':flag_grp, 'sgl': sgl_grp, 'sgb':sgb_grp, 'gl':gl_grp, 'gb':gb_grp, 'ra':ra_grp, 'dec':dec_grp, 'Ks':Ks_grp, 'Vls':Vls_grp,  \
               'nest':nest_grp, 'dcf2':dcf2_grp, 'ed':ed_grp, 'mDist':mDist_grp, 'mDistErr':mDistErr_grp, 'No_Galaxies':No_Galaxies_grp} 
     
     GroupList.append(GroupNode(data=mygroup))
     

     
     if i<N and flag[i]==2: 
          gr = 2
     else:
	  break
   
   return GroupList
# **************************************
class GroupNode:
   
   def __init__(self, data=None):
     
     
     
      self.dist = 0
      ########self.mDist = 0
      
      self.sgl        = 0
      self.sgb        = 0
      self.gl         = 0
      self.gb         = 0 
      self.ra         = 0
      self.dec        = 0
      self.Ks         = 0
      self.logK       = 0
      
      self.M_v2       = 0
      self.R_2t2      = 0
      self.r2t        = 0 
      self.r1t        = 0
      
      self.flag   = 0   #     
      
      if data == None:
	self.id        = 0
	self.Vls       = 0
	self.sigma     = 0 
	self.nest      = 0
	self.mDist     = 0
	self.mDistErr  = 0
	self.Vls_list  = None
	self.ngal      = 0
      
        return
      
      #self.data = data   ### ????
      
      id              = data['pgc']
      flag            = data['flag']
      sgl             = data['sgl']
      sgb             = data['sgb']
      gl              = data['gl']
      gb              = data['gb']  
      ra              = data['ra']
      dec             = data['dec']
      Ks              = data['Ks']
      nest            = data['nest']
      dcf2            = data['dcf2']
      ed              = data['ed']
      mDist           = data['mDist']
      mDistErr        = data['mDistErr']
      self.Vls_list   = data['Vls']
      ngal            = data['No_Galaxies']
      
      
      #### ??? Radius = r1t, Mass=M_v2, R_2t2
      self.id         = id[0]
      self.ngal       = ngal[0]
      self.Vls, self.sigma   = self.v_ave(self.Vls_list[1:])
      self.nest       = nest[0]
      self.mDist, self.mDistErr = self.dist_av(dcf2, ed)

      self.initialize(sgl, sgb, gl, gb, ra, dec, Ks)
   
   
   # ************
   def initialize(self, sgl, sgb, gl, gb, ra, dec, Ks):
     
     Ltot = 0
     Msgl=0; Msgb=0
     Mgl=0; Mgb=0
     Mra=0; Mdec=0
     for i in range(1, len(Ks)):
       
       L = 10**m_logK(Ks[i], self.Vls, distance=self.mDist)
       Msgl += L*sgl[i]
       Msgb += L*sgb[i]
       Mgl += L*gl[i]
       Mgb += L*gb[i]
       Mra += L*ra[i]
       Mdec += L*dec[i]       
       Ltot += L
     
     self.sgl  = Msgl/Ltot
     self.sgb  = Msgb/Ltot
     self.gl   = Mgl/Ltot
     self.gb   = Mgb/Ltot
     self.ra   = Mra/Ltot
     self.dec  = Mdec/Ltot
     
     self.logK = log10(Ltot)
     Mk_sun = 3.28   # K-band
     M = Mk_sun - (self.logK / 0.4)
     Dist_v = self.mDist
     if Dist_v==0: 
       Dist_v =  self.Vls / H0
       if Dist_v<1. : Dist_v=1
     
     self.dist =   Dist_v
     self.Ks = M + 5*log10(Dist_v) + 30 - 5
     
     self.M_v2 = Mass(10**self.logK)
     self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
     
     self.r2t = sqrt(1.5)*self.R_2t2
     self.r1t = 3.5*self.r2t
   
   # ************
   def dist_av(self, dcf2, ed):
     
     sum_dist = 0
     sum_edist = 0
     mD=0; emD=0
     n = 0
     for i in range(1, len(dcf2)):
       
       if dcf2[i] != 0 and ed[i] != 0:
	 
	 err =  dcf2[i]*ed[i]
	 sum_dist += dcf2[i]/(err**2)
	 sum_edist += 1/(err**2)
	 n += 1
     
     if sum_edist!=0:
       mD = sum_dist/sum_edist
       if mD!=0:
         emD = sqrt(1/sum_edist)/mD
     
     return mD, emD
     
   
   
   
   def v_ave(self, v_list):
     
     m = 0 
     sum_v = 0
     sum_v2 = 0 
     
     for v in v_list:
       if v != 0:
	 sum_v += v
	 sum_v2 += v**2
	 m +=1
     
     if m != 0 : 
        mean_v = sum_v/m
        sigma_v = sqrt((sum_v2/m)-mean_v**2)
     else:
        mean_v = 0
        sigma_v = 0
     
     return mean_v, sigma_v
     
#################################################################
def m_logK(m, Vls, distance=0):
    
    Mk_sun = 3.28   # K-band
    
    if distance == 0: 
       distance = Vls / H0  # Mpc    
       if distance < 1:
          distance = 1 # Mpc
    
    M = m - 5*log10(distance) - 30 + 5
    logK = -0.4 * (M - Mk_sun)
    
    if Vls==0 and distance==0:
      logK = 0
    
    
    return logK
# **************************************

def SuperGroup(GroupList, old_SuperList=None):
  
  
  if old_SuperList != None:
    
    inGame = []
    Lum_mass = []
    for SuperGroup in old_SuperList:
      Lum_mass.append(SuperGroup[0].M_v2)
    indices = np.argsort(Lum_mass)
    for i in indices[::-1]: 
      for group in old_SuperList[i][1]:
          inGame.append(group)
          
    for group in GroupList:
      if group.flag == 0:
	inGame.append(group)
	
    for group in GroupList:
      group.flag = 0

  else:
  # -------------------------

    inGame = GroupList[:]
    
    Lum_mass = []
    for i in range(0, len(inGame)): Lum_mass.append(inGame[i].M_v2)
    Lum_mass = np.asarray(Lum_mass)
    indices = np.argsort(Lum_mass)
    
    tmp = []
    for i in indices[::-1]: tmp.append(inGame[i]) 
    inGame = tmp
    
  # -------------------------
  
  
  SuperList = []
  N = len(inGame)
  print N
  
  p = 0 
  while p < N-1:
    
    Players = []
    Players.append(inGame[p])
    
    q = (p+1)%N
    pivot = p
    head = Players[0]
    Bol = False
    while q!=pivot:  
      
      if touch(head, inGame[q], coeff=1.0):
	Bol = True
	Players.append(inGame[q])
	head =  Gjoint(Players)
	
	inGame.pop(q)
	N-=1
	if p>=q:  p-=1
  
	if len(Players) == 2:
	   inGame.pop(p)
	   N-=1	
	   if q>p: q = (q-1)%N    
	
	q = (q-1)%N
	pivot = q  
  
      q = (q+1)%N
    if Bol: 
      SuperList.append([head, Players])
    p+=1 
    
  return SuperList
#################################################################
def mergSGroup(SuperList):
  
    NoGroups = len(SuperList)
 
    i = 0
    while i < NoGroups:
      j = i+1
      while j < NoGroups-1:
	
	if touch(SuperList[i][0], SuperList[j][0], coeff=1.0):
	  
	  Players = []
	  
	  for group in SuperList[i][1]:
	    Players.append(group)
	  
	  for group in SuperList[j][1]:
	    Players.append(group)
	    
	  head =  Gjoint(Players)	  
	  
	  SuperList[i] = [head, Players]
	  SuperList.pop(j)
	  NoGroups-=1
	  i=0
	  j=0
        j+=1
      i+=1
    
    return SuperList  
  

#################################################################

def Gjoint(Players):
  
   if len(Players)==0: return None

   if len(Players)==1: return Players[0]

   Lum_mass = []
   for i in range(0, len(Players)):
        Lum_mass.append(Players[i].Ks)
   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass)
   
   id = Players[indices[0]].id
   
   head = None
   for i in indices:
       head = Joint(head, Players[i], ID = id+10000000000)
       
   return head
 
#################################################################


def Joint(gr1, gr2, ID = 0):
  
   if gr1 == None and gr2 != None:
     return gr2
   if gr1 != None and gr2 == None:
     return gr1
   if gr1 == None and gr2 == None:
     return None  
   
   gr = GroupNode()
   gr.id = ID
   
   L1 = 10**gr1.logK
   L2 = 10**gr2.logK
   
   if L1 == 1: L1 = 0
   if L2 == 1: L2 = 0
   
   L_tot = L1 + L2
   
   gr.logK = log10(L_tot)
   
   
   d1 = gr1.dist
   d2 = gr2.dist
   
   sgl1 = gr1.sgl
   sgl2 = gr2.sgl
   sgb1 = gr1.sgb
   sgb2 = gr2.sgb
   
   gl1 = gr1.gl
   gl2 = gr2.gl
   gb1 = gr1.gb
   gb2 = gr2.gb 
   
   ra1 = gr1.ra
   ra2 = gr2.ra
   dec1 = gr1.dec
   dec2 = gr2.dec   

   v1 = gr1.Vls
   v2 = gr2.Vls   
   gr.Vls  = (L1*v1 + L2*v2)/L_tot
   
   gr.sgl, gr.sgb, d     = barycenter(L1,  sgl1,  sgb1, d1, L2,  sgl2,  sgb2, d2)
   gr.gl,  gr.gb,  d0    = barycenter(L1, gl1, gb1, d1, L2, gl2, gb2, d2)
   gr.ra,  gr.dec, d00   = barycenter(L1, ra1, dec1, d1, L2, ra2, dec2, d2)
   ###if abs(d- d0)>1 :
     ###print gr1.id, gr2.id
     ###print "  ", d, d0, d00
     
   gr.dist = np.median([d,d0,d00])
   Mk_sun = 3.28   # K-band
   M = Mk_sun - (gr.logK / 0.4)
   gr.Ks = M + 5*log10(gr.dist) + 30 - 5
   
   gr.M_v2 = Mass(10**gr.logK)
   gr.R_2t2 = 0.215*((gr.M_v2/1.E12)**(1./3))  # Mpc
   gr.r2t = sqrt(1.5)*gr.R_2t2
   gr.r1t = 3.5*gr.r2t
   

   gr.flag = 2
   gr1.flag = 1
   gr2.flag = 1
   
   gr.ngal = gr1.ngal + gr2.ngal
   gr.Vls_list = np.asarray([gr.Vls])
   gr.Vls_list = np.concatenate((gr.Vls_list, gr1.Vls_list[1:], gr2.Vls_list[1:]))
   tmp, gr.sigma = gr.v_ave(gr.Vls_list[1:])
   
   
   return gr
  
#################################################################
#################################################################   
def SGroupwrite(outfile, SuperList, GroupList):
  
   NoSGroups = len(SuperList)
   
   myTable = Table()
   
   empty = []
   myTable.add_column(Column(data=empty,name='pgc', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='flag', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='ra', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='dec', format='%0.4f', length=10))
   myTable.add_column(Column(data=empty,name='gl', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='gb', format='%0.4f', length=10))    
   myTable.add_column(Column(data=empty,name='sgl', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='sgb', format='%0.4f', length=10))

   myTable.add_column(Column(data=empty,name='Ks', format='%0.2f'))
   myTable.add_column(Column(data=empty,name='logK', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='Vls', format='%0.0f'))
   myTable.add_column(Column(data=empty,name='dist', format='%0.2f'))
   myTable.add_column(Column(data=empty,name='mDist', format='%0.2f'))
   myTable.add_column(Column(data=empty,name='mDistErr', format='%0.2f'))

   myTable.add_column(Column(data=empty,name='sigmaP_dyn', format='%0.1f'))
   myTable.add_column(Column(data=empty,name='sigmaP_lum', format='%0.1f'))
  
   myTable.add_column(Column(data=empty,name='Mv_lum', format='%1.2e'))
   myTable.add_column(Column(data=empty,name='R2t_lum', format='%0.3f'))
   myTable.add_column(Column(data=empty,name='r1t_lum', format='%0.3f'))  
   myTable.add_column(Column(data=empty,name='tX_lum', format='%1.2e'))  
   
   myTable.add_column(Column(data=empty,name='No_Galaxies',dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='nest', dtype=np.dtype(int)))
   
   for i in range(0, NoSGroups):  # for all groups
        table_row(myTable, SuperList[i][0])
        for Group in SuperList[i][1]:
	  table_row(myTable, Group)
   for Group in GroupList:
     if Group.flag<=0: 
       table_row(myTable, Group)
   
   
   pgc = 999999999; 
   ra = 999.9999; dec=-99.99;
   gl = ra; gb = dec
   sgl = ra; sgb=dec
   Ty = -100000.00; B_mag=Ty
   Ks= 99.99
   logK = 99.9999
   Vls = 9999
   dcf2 = 99.99
   ed = 9.99
   Mv_dyn = 9.99E99; Mv_lum = Mv_dyn
   tX_dyn = Mv_lum; tX_lum=Mv_lum
   nest = 9999999   
   
   flag = 0
   mDist = 0
   mDistErr = 0
   dist = 0
   sigmaP_dyn = 0
   sigmaP_lum = 0 
   R2t_lum = 0
   r1t_lum = 0
   subGalaxies = 0
   
   myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ks,logK,Vls, dist, \
	       mDist, mDistErr, sigmaP_dyn, sigmaP_lum, \
	          Mv_lum, R2t_lum, r1t_lum, tX_lum, subGalaxies, nest])
   
   
   myTable.write(outfile, format='ascii.fixed_width',delimiter='|', bookend=False)
   
   ### removing the last line, (it sits o adjust the column wodths)
   command =  ["csh", "remove_lastline.csh", outfile]
   subprocess.call(command)  
   #print myTable
#################################################################

def table_row(myTable, Group):
  
        
        pgc    = Group.id  
        flag   = Group.flag
        ra = Group.ra
        dec = Group.dec
        gl     = Group.gl  
        gb     = Group.gb
        sgl    = Group.sgl  
        sgb = Group.sgb  
        
        Ks = Group.Ks
        logK = Group.logK  
        
        Vls = Group.Vls  
        dist = Group.dist
        
        mDist = Group.mDist
        mDistErr = Group.mDistErr
        
        subGalaxies = Group.ngal  
        sigmaP_dyn = Group.sigma  
        nest = Group.nest  
        
        Mv_lum = Group.M_v2
        R2t_lum = Group.R_2t2
        r1t_lum = Group.r1t
        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqrt(2.5)

        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ks,logK,Vls, dist, \
	       mDist, mDistErr, sigmaP_dyn, sigmaP_lum, \
	          Mv_lum, R2t_lum, r1t_lum, tX_lum, subGalaxies, nest])

#################################################################

if __name__ == '__main__':
  

   filename = "south.iter.2.v31.group" 
   print filename
   GroupList = readgrouplist(filename)   
   
   #print GroupList[0].Vls, GroupList[0].sigma
   #print GroupList[0].mDist
   #print GroupList[0].mDistErr
   #print GroupList[0].logK, GroupList[0].Ks
   #print GroupList[0].ra, GroupList[0].dec
   #print GroupList[0].gl, GroupList[0].gb
   #print GroupList[0].sgl, GroupList[0].sgb
   #print GroupList[0].R_2t2, GroupList[0].M_v2
   #print "Dist: ", GroupList[0].dist, GroupList[0].ngal
   
   SuperList = SuperGroup(GroupList)
   SuperList = SuperGroup(GroupList, old_SuperList=SuperList)
   
   for qp in range(10):
      mergSGroup(SuperList)   
   
    
   
   #for SGroup in SuperList:
     #print
     #print
     #print SGroup[0].id, SGroup[0].dist, SGroup[0].r1t
     #for Group in SGroup[1]:
       #print "  ",Group.id, Group.dist, Group.r1t
  
   print "\n total super groups #:", len(SuperList)
   
   
   SGroupwrite('south.iter.2.v31.supergroup', SuperList, GroupList)

