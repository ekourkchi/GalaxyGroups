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
from scipy  import interpolate

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
LtoM_func = None


# Returns mass in solar unit, given the absolute K-band luminosity
def Mass(L_k):
  
  if L_k==0:
    return 0
  
  L = L_k / 1.E10
  
  if L <= 0.0927:
    MtoL = 32.0*(L**-0.5)
  elif L > 0.0927 and L < 4.423:
    #MtoL = 58.0*(L**-0.25)
    
    return LtoM_func(L_k)
  elif L >= 4.423:
    MtoL = 32*(L**0.15)
  
  Mass_out = L_k * MtoL
  
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
  
   if d1==0 or d2==0:
       
       dd1 = 1.
       dd2 = 1.
   else: 
       dd1 = d1
       dd2 = d2
 
   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)
   
   x1 = dd1 * cl1 * cb1
   y1 = dd1 * sl1 * cb1
   z1 = dd1 * sb1
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = dd2 * cl2 * cb2
   y2 = dd2 * sl2 * cb2
   z2 = dd2 * sb2   

   if L1==0 and L2==0:
     x = 0.5*(x1+x2)
     y = 0.5*(y1+y2)
     z = 0.5*(z1+z2)
   else: 
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

def galJoint(galNode1, galNode2, ID):
   
   if galNode1==None and galNode2!=None:
     return galNode2
   if galNode2==None and galNode1!=None:
     return galNode1
   if galNode2==None and galNode1==None:
     return None
   

   
   L1 = 10**galNode1.logK
   L2 = 10**galNode2.logK
   
   if L1 == 1: L1 = 0
   if L2 == 1: L2 = 0
   
   L_tot = L1 + L2
   
   try:
     logK_tot = log10(L_tot)
   except:
     logK_tot = 0
     
     
   
   
   d1 = galNode1.dcf2
   d2 = galNode2.dcf2
   
   sgl1 = galNode1.sgl
   sgl2 = galNode2.sgl
   sgb1 = galNode1.sgb
   sgb2 = galNode2.sgb
   
   gl1 = galNode1.gl
   gl2 = galNode2.gl
   gb1 = galNode1.gb
   gb2 = galNode2.gb 
   
   ra1 = galNode1.ra
   ra2 = galNode2.ra
   dec1 = galNode1.dec
   dec2 = galNode2.dec   
      
   
   
   
   v1 = galNode1.Vls
   v2 = galNode2.Vls
   
   if d1==0 and d2==0:
     d1 = 0.
     d2 = 0
     d_final=0.
   elif d1!=0 and d2==0:
     d2 = galNode2.Vls/H0
     d_final=0
   elif d2!=0 and d1==0:
     d1 = galNode1.Vls/H0
     d_final=0
   else:
     d_final=1
     
     
   sgl, sgb, d   = barycenter(L1,  sgl1,  sgb1, d1, L2,  sgl2,  sgb2, d2)
   gl, gb, d = barycenter(L1, gl1, gb1, d1, L2, gl2, gb2, d2)
   ra, dec, d = barycenter(L1, ra1, dec1, d1, L2, ra2, dec2, d2)
   
   if d_final==0:
      d = 0
      dm = 0
   else:
      dm = 5*log(d)+25


   n1 = galNode1.subGalaxies
   n2 = galNode2.subGalaxies   
   Vls = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)
   if galNode1.v_av == 0 and galNode2.v_av != 0 :
     Vls = galNode2.v_av
     if L1 == 0:
       galNode1.logK = m_logK(galNode1.Ks, Vls, d=galNode1.dcf2)
       
   if galNode1.v_av != 0 and galNode2.v_av == 0 :
     Vls = galNode1.v_av   
     if L2 == 0:
       galNode2.logK = m_logK(galNode2.Ks, Vls, d=galNode2.dcf2) 

   newNode = GalxyNode(ID, gl, gb, sgl, sgb, Vls, 0, 0, d, 0, dm, 0)
   newNode.logK = logK_tot
   newNode.ra = ra
   newNode.dec = dec
   newNode.coordinate_src = 'GRPhead'
   newNode.Ty_src  =  'GRPhead'
   newNode.Ks_src  =  'GRPhead'
   newNode.Vls_src =  'GRPhead'
   newNode.objname =  'GRPhead'
   
   newNode.subGalaxies = n1 + n2
   newNode.level = max([galNode1.level, galNode2.level]) + 1
   
   galNode1.topLevel = newNode.level
   galNode2.topLevel = newNode.level
   
   newNode.sumDist = galNode1.sumDist + galNode2.sumDist
   newNode.sumError = galNode1.sumError + galNode2.sumError
   
   if newNode.sumDist !=0 and newNode.sumError!=0:
        meanDist = newNode.sumDist/newNode.sumError
        meanDistErr = sqrt(1/newNode.sumError)/newNode.sumDist
        mDM  = (log10(meanDist)+5.)*5.
        meDM = meanDistErr / (0.2*log(10.))
   else:
        mDM = 0.
        meDM = 0.
        meanDist = 0.
        meanDistErr = 0.
    
   
   
   if galNode1.mDist != 0 and galNode2.mDist != 0:
     if L1 > L2: 
       newNode.mDist = galNode1.mDist
       newNode.mDistErr = galNode1.mDistErr
       newNode.mDM = galNode1.mDM
       newNode.meDM = galNode1.meDM
     else:
       newNode.mDist = galNode2.mDist
       newNode.mDistErr = galNode2.mDistErr
       newNode.mDM = galNode2.mDM
       newNode.meDM = galNode2.meDM
       
   elif galNode1.mDist != 0:
     newNode.mDist = galNode1.mDist
     newNode.mDistErr = galNode1.mDistErr
     newNode.mDM = galNode1.mDM
     newNode.meDM = galNode1.meDM
     
   elif galNode2.mDist != 0:
     newNode.mDist = galNode2.mDist
     newNode.mDistErr = galNode2.mDistErr
     newNode.mDM = galNode2.mDM
     newNode.meDM = galNode2.meDM
   else:
     newNode.mDist = meanDist
     newNode.mDistErr = meanDistErr
     newNode.mDM = mDM
     newNode.meDM = meDM

   
   # Brighter galaxy is the left child
   if L1 >= L2:
      newNode.left = galNode1
      newNode.right = galNode2
      newNode.nest = galNode1.nest
   else:
      newNode.left = galNode2
      newNode.right = galNode1
      newNode.nest = galNode2.nest
   

   newNode.v_av = Vls
   newNode.v2_av = (n1*galNode1.v2_av + n2*galNode2.v2_av) / (n1+n2)
   if galNode1.v2_av == 0 and galNode2.v2_av != 0 :
     newNode.v2_av = galNode2.v2_av
   if galNode1.v2_av != 0 and galNode2.v2_av == 0 :
     newNode.v2_av = galNode1.v2_av   
   
   
   
   if (newNode.v2_av - newNode.v_av**2) > 0:
      newNode.sigma =  sqrt(newNode.v2_av - newNode.v_av**2) 
   
   
   newNode.R_theta = Theta_max(newNode)
   

   if newNode.sigma == 0:
     sig = 1
   else:
     sig = newNode.sigma
   

   mass = Mass(L_tot)
   newNode.M_v2 = mass
   if np.isinf(newNode.M_v2): 
	    print 'Grouping ...'
	    sys.exit()
   newNode.R_2t2 = 0.215*((mass/1.E12)**(1./3))  # Mpc

   return newNode
# **************************************
# The calls definition of a galaxy
# each node contains all esseential property of a galaxy
# when gaalxies get connected along a tree, the new nodes
# are defnided as new entities, which are new galaxies in
# the context of this code
class GalxyNode:
  
  # field variables
  id = 0.
  l = 0.
  b = 0.
  gl = 0.
  gb = 0.
  Vls = 0.
  logK = 0.
  subGalaxies = 1
  level = 0
  R_theta = 0.
  nest = 0
  Ty = 0.
  v_av = 0.
  v2_av = 0.
  sigma = 0.  # velocity dispersion (if a node contains several gaalxies at the bottom of a tree)
  dcf2 = 0.
  ed = 0.
  dcf2Copy = 0.
  edCopy = 0.
  mDist = 0.
  mDistErr = 0.
  Ks = -100000.0
  inGroup = 0.
  
  M_v2 = 0.
  R_2t2 = 0.
  
  sumDist = 0.
  sumError = 0.
  
  
  ra  =  0.
  dec = 0.
  coordinate_src = 'NotSet'
  Ty_src  = 'NotSet'
  B_mag   = -100000.0
  Ks_src  = 'NotSet'
  Vls_src = 'NotSet'
  objname = 'NotSet'
  
  flag = 0
  
  left = None
  right = None
  
  # Class constructor
  def __init__(self, id, gl, gb, sgl, sgb, Vls, Ks, Ty, dcf2, ed, dm, edm):
    
    
    #dcf2 = 0
    self.id = id
    self.gl = gl
    self.gb = gb    
    self.sgl = sgl
    self.sgb = sgb
    self.Vls = Vls
    self.Vhelio = 0.
    self.Ks = Ks
    self.Ty = Ty
    
    
    if id in [42666]:
      dcf2 = 16.
      ed = 0.20
      dm = (log10(dcf2)+5.)*5.
      edm =  ed / (0.2*log(10.))
    
    
    self.dcf2 = dcf2
    self.ed = ed
    self.dm = dm
    self.edm = edm
    
    
    
    
    self.dcf2Copy = dcf2
    self.edCopy = ed    
    self.dmCopy = dm
    self.edmCopy = edm        
    
    self.mDist = dcf2
    self.mDistErr = ed
    self.mDM = dm
    self.meDM = edm    
    
    self.left = None
    self.right = None
    
    self.subGalaxies = 1
    self.level = 0
    self.Rgroup = 0.  
    self.R_theta = 0.
    self.nest = id
    self.sigma = 0. 
    self.v_av = Vls
    self.v2_av = Vls**2
    
    # 0: if the galaxy is NOT in a group
    # 1: if the galaxy falls in a group
    self.inGroup = 0.
    
    self.logK = m_logK(Ks, Vls, d=dcf2)
 
    self.M_v2 = Mass(10**self.logK)
    self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc


    
    if dcf2!=0 and ed!=0:
        
        err = dcf2*ed
        self.sumDist = 1. * dcf2 / (err**2)
        self.sumError = 1. / (err**2)
    else:
        self.sumDist = 0.
        self.sumError = 0.
        
    
    
  def setMeanDist(self, Dist, errDist, GRP_vel = 0, dist=0):
    
      self.mDist = Dist
      self.mDistErr = errDist
      if GRP_vel == 0 : 
	vel = self.Vls
      else: vel = GRP_vel
      self.logK = m_logK(self.Ks, vel, d=dist)
      self.M_v2 = Mass(10**self.logK)
      self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
      
    


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
    
    if Vls == 0 and d==0:
      logK = 0
    
    
    return logK
# **************************************
# It returns all leaf-nodes at the bottom of a gaalxy tree.
# root: is the root of the tree
# It works recursively, so for large trees, it might 
# exhaust the recursion stack, and the codes crashe
# so some exception handling is required in this case
def NodeLeaf(root):
  list = [root]
  NodeLeafCore(root, list)
  return list

def NodeLeafCore(root, list):  
  if root.left == None:

    list.append(root)
  else:
    NodeLeafCore(root.left, list)
    NodeLeafCore(root.right, list)
    
    
   


################################################################
def readgalList(table, skyPatch='none', sgl_min=0, sgl_max=0, sgb_min=0, sgb_max=0, Vls_max=4000, ignore_list=[]):

  id   = table['pgc']
  sgl  = table['sgl']
  sgb  = table['sgb']
  Vls  = table['Vls']
  Vhelio = table['Vhelio']
  Ks   = table['Ks'] 
  dcf2 = table['dcf2'] 
  ed   = table['ed'] 
  dm = table['DM'] 
  edm   = table['eDM'] 
  
  Ty   = table['Ty'] 
  
  gl  = table['gl']
  gb  = table['gb']  
  
  
  ra  =  table['ra']
  dec = table['dec']
  coordinate_src = table['coordinate_src']
  Ty_src  = table['Ty_src']
  B_mag   = table['B_mag']
  Ks_src  = table['Ks_src']
  Vls_src = table['Vls_src']
  objname = table['objname']
   
  N_galaxies = len(id)


  print "Reading the data file .... "
  galList = []
  

    

  i = 0
  
  for i in range(N_galaxies):
    
    
    
    
    ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)
    
    
    Bool = False
    
    
    inHydra = False
    if (180./pi)*angle(sgl[i], sgb[i], 139.3599, -37.7051) < 3 and Vhelio[i] < 6000.: inHydra=True

    
    
    
    inCentaurus = False
    if (180./pi)*angle(sgl[i], sgb[i], 156.2981, -11.6740) < 3 and Vhelio[i] < 6000.: inCentaurus=True

    if  skyPatch=='manual':
      if (Vls[i] < Vls_max or inHydra or inCentaurus) and sgl[i]<sgl_max and sgl[i]>sgl_min and sgb[i]>sgb_min and sgb[i]<sgb_max:
        Bool = True
    elif skyPatch=='local':
      if (Vls[i] < Vls_max or inHydra or inCentaurus):
	Bool = True
    else:
      if (Vls[i] < Vls_max or inHydra or inCentaurus): 
	  Bool = True
	  
    if id[i] == 31599:
      print id[i], (180./pi)*angle(sgl[i], sgb[i], 139.3599, -37.7051) , Vls[i], inHydra, Bool
      
    if Bool:         

	   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[i], Vls[i], Ks[i], Ty[i], dcf2[i], ed[i], dm[i], edm[i])
	   node.ra = ra[i]
	   node.dec = dec[i]
	   node.coordinate_src = coordinate_src[i]
	   node.Ty_src = Ty_src[i]
	   node.B_mag =  B_mag[i]
	   node.Ks_src = Ks_src[i]
	   node.Vls_src = Vls_src[i]
	   node.objname = objname[i]
	   node.Vhelio = Vhelio[i]
           galList.append(node)
           
           if id[i] in ignore_list or Ks[i] < 0 or Vls[i]==0:
	     #print "dis_include:", id[i]
	     node.inGroup = -1
           
  print "No. of galaxies: ", len(galList)
  print "Data file loaded .... "
  

  
  
  return galList
################################################################  
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

def xymap(l,b,l0):
  
  t = theta(b)
  x = 180.-((2*90.)/180.)*(l0-l)*cos(t)
  y = 90.*sin(t)
  
  return x, y 
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

  plt.plot(X, Y, '-', markersize = 1, linewidth=2., color='#0D1891')
  
################################################################ 

def plot_galaxies(inFile, l0):
  
  table = np.genfromtxt( inFile , delimiter=',', filling_values=0, names=True, dtype='float64')
  
  #print table.dtype
  id   = table['pgc']
  gl  = table['sgl']
  gb  = table['sgb']
  N_galaxies = len(id)
  
  X0 = []
  Y0 = []
  for i in range(0, N_galaxies):
    x, y = xymap(gl[i],gb[i],l0)
    X0.append(x)
    Y0.append(y)
  plt.plot(X0, Y0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  
  
################################################################ 
def groupPlot(northGList, southGList):
  
  l0 = 180
  
  fig = plt.figure(figsize=(12, 12*0.5), dpi=100)
  ax = fig.add_axes([0.13, 0.13, 0.80,  0.80]) 
  plt.ylim(-100,100)
  plt.xlim(380,-20)
  plt.xlabel("SGL (deg)", fontsize=20)
  plt.ylabel("SGB (deg)", fontsize=20)
  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)
  
  plot_border(l0)
  plot_galaxies('AllSky.north.csv', l0)
  plot_galaxies('AllSky.south.csv', l0)
  
  # Virgo Border (6.8 deg)
  u = np.arange(0,1,0.010)
  X = u*0.
  Y = u*0.
  for i in range(0,len(u)):
     x = 6.8*np.cos(u[i]*2*pi) + sglV
     y = 6.8*np.sin(u[i]*2*pi) + sgbV
     X[i], Y[i] = xymap(x,y, l0)
  
  line00, = plt.plot(X,Y, 'r.', markersize = 2)   
  line00.set_dashes([2, 2]) 
  
  groupPlotak(northGList, l0)
  groupPlotak(southGList, l0)
##########

def groupPlotak(G_list, l0):
  
 
  NoGroups = len(G_list)
  print "Number of groups: ", NoGroups
  
  col = 0

  for i in range(0, NoGroups):  # for all groups

    
    Key = True
    if Key:
	random.seed(G_list[i][0].nest )
	Red, Blue, Green = random.random(), random.random(), random.random()
	
	
	r = G_list[i][0].R_theta
	
	d_theta = 0.001
	u = np.arange(0,2*pi,d_theta)
	X = u*0.
        Y = u*0.
        for q in range(0,len(u)):
           x = r*np.cos(u[q]) + G_list[i][0].sgl
           y = r*np.sin(u[q]) + G_list[i][0].sgb
           X[q], Y[q] = xymap(x,y, l0)

	
	if r <= 6.3:
	  line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Blue, Green)) 
	  line.set_dashes([8, 3]) 
	

	X  = []
	Y = []
	

	for j in range(1, len(G_list[i])):

	    x, y = xymap(G_list[i][j].sgl,G_list[i][j].sgb, l0)
	    X.append(x)
	    Y.append(y)
	    
	plt.plot(X, Y, 'o', markersize = 3, color=(Red, Blue, Green), markeredgecolor = 'none')
        

#################################################################
################################################################# 
def groupWrite(outfile, G_list, galList):
  
  
  # finding the leave nodes
  #galList = NodeLeaf(root)
  #galList.pop(0)
  
  #for grp in G_list:
    #for gal in grp:
      #if gal.id == 31599:
        #print 'I have in my catalog'
      
      
  NoGroups = len(G_list)
  #print "Number of groups: ", NoGroups  

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
  myTable.add_column(Column(data=empty,name='Ty'))
  myTable.add_column(Column(data=empty,name='B_mag', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='logK', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=empty,name='Vhelio', format='%0.0f'))
  myTable.add_column(Column(data=empty,name='dcf2', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='ed', format='%0.2f'))

  myTable.add_column(Column(data=empty,name='mDist', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='mDistErr', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='R_theta', format='%0.5f'))
  myTable.add_column(Column(data=empty,name='sigmaP_dyn', format='%0.1f'))
  myTable.add_column(Column(data=empty,name='sigmaP_lum', format='%0.1f'))
  
  myTable.add_column(Column(data=empty,name='Mv_dyn', format='%1.2e'))
  myTable.add_column(Column(data=empty,name='Mv_lum', format='%1.2e'))
  myTable.add_column(Column(data=empty,name='Rg_angular', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='Rg_dyn', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='R2t_dyn', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='R2t_lum', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='tX_dyn', format='%1.2e'))
  myTable.add_column(Column(data=empty,name='tX_lum', format='%1.2e'))  
  myTable.add_column(Column(data=empty,name='No_Galaxies',dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='nest', dtype=np.dtype(int)))

  myTable.add_column(Column(data=empty,name='coordinate_src', dtype='S15'))
  myTable.add_column(Column(data=empty,name='Ty_src', dtype='S15'))
  myTable.add_column(Column(data=empty,name='Ks_src', dtype='S15'))
  myTable.add_column(Column(data=empty,name='Vls_src', dtype='S15'))
  myTable.add_column(Column(data=empty,name='objname', dtype='S35'))

  myTable.add_column(Column(data=empty,name='DM', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='eDM', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='mDM', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='meDM', format='%0.2f'))
  
  
  #print "# of all groups: ", NoGroups
  
  for i in range(0, NoGroups):  # for all groups
    
   
    if G_list[i][0].Vls <= 3500:
       dist = G_list[i][0].mDist
       if dist == 0 : dist = G_list[i][0].Vls / H0
       if dist<1: dist = 1
       Rg =  dist * Rg_radian(G_list[i][1:])
       Rg_angular = dist * G_list[i][0].R_theta
       for j in range(0, len(G_list[i])):  # for all galaxies
        
        
        galaxy = G_list[i][j]
        flag = galaxy.flag
        pgc = galaxy.id  
        
        
        
        
        gl = galaxy.gl  
        gb = galaxy.gb
        sgl = galaxy.sgl  
        sgb = galaxy.sgb  
        Vls = galaxy.Vls  
        Vhelio = galaxy.Vhelio
        logK = galaxy.logK  
        Ks = galaxy.Ks
        Ty = galaxy.Ty  
        dcf2 = galaxy.dcf2Copy
        ed = galaxy.edCopy
        dm = galaxy.dmCopy
        edm = galaxy.edmCopy
        mdm = galaxy.mDM
        medm = galaxy.meDM
        
        
        ra = galaxy.ra
        dec = galaxy.dec
        coordinate_src = galaxy.coordinate_src
        Ty_src = galaxy.Ty_src
        B_mag = galaxy.B_mag
        Ks_src = galaxy.Ks_src
        Vls_src = galaxy.Vls_src
        objname = galaxy.objname
        
        mDist = galaxy.mDist
        mDistErr = galaxy.mDistErr
        

        
        subGalaxies = G_list[i][0].subGalaxies  
        R_theta = G_list[i][0].R_theta  
        sigmaP_dyn = G_list[i][0].sigma  
        nest = G_list[i][0].nest  
        Mv_lum = galaxy.M_v2
        R2t_lum = galaxy.R_2t2
        
        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqrt(2.5)
        if j != 0: 
	  Rg = 0  # dynamical virial radius
	  R_theta = 0
	  Rg_angular = 0
        Mv_dyn = (1.E9 * (2.5*pi/G/2.) * (Mpc_km*Rg) * sigmaP_dyn**2)/M_sun  # solar mass
        
        if sigmaP_dyn == 0: 
	  tX_dyn = 0 
	else:
          tX_dyn = Mpc_km*0.5*pi*Rg/sigmaP_dyn/sqrt(2.5)
        R2t_dyn = (Rg*pi*0.5)/1.05/sqrt(1.5)
       
 
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,Ks,logK,Vls, Vhelio, dcf2, ed, \
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_lum, \
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, tX_dyn, tX_lum, subGalaxies, nest, \
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname, dm, edm, mdm, medm])
	

  
  # writing individual galaxies
  for galaxy in galList:
    


      if galaxy.inGroup <= 0. and galaxy.Vls<=3500:
	flag = galaxy.inGroup
	galaxy.flag = galaxy.inGroup
        pgc = galaxy.id  
        gl = galaxy.gl  
        gb = galaxy.gb
        sgl = galaxy.sgl  
        sgb = galaxy.sgb  
        Vls = galaxy.Vls  
        Vhelio = galaxy.Vhelio
        logK = galaxy.logK  
        Ks = galaxy.Ks
        Ty = galaxy.Ty  
        dcf2 = galaxy.dcf2Copy
        ed = galaxy.edCopy
        dm = galaxy.dmCopy
        edm = galaxy.edmCopy
        mdm = galaxy.mDM
        medm = galaxy.meDM
        
        ra = galaxy.ra
        dec = galaxy.dec
        coordinate_src = galaxy.coordinate_src
        Ty_src = galaxy.Ty_src
        B_mag = galaxy.B_mag
        Ks_src = galaxy.Ks_src
        Vls_src = galaxy.Vls_src
        objname = galaxy.objname
        
        mDist = galaxy.mDist
        mDistErr = galaxy.mDistErr

        
        subGalaxies = galaxy.subGalaxies  
        R_theta = galaxy.R_theta  
        sigmaP_dyn = galaxy.sigma  
        nest = galaxy.nest  
        Mv_lum = galaxy.M_v2
        R2t_lum = galaxy.R_2t2
        
        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqrt(2.5)
        Rg = 0  # dynamical virial radius
        Mv_dyn = 0  # solar mass
        tX_dyn = 0
        R2t_dyn = 0
        Rg_angular = 0
        R_theta = 0 
       
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,Ks,logK,Vls, Vhelio, dcf2, ed, \
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_lum, \
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, tX_dyn, tX_lum, subGalaxies, nest, \
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname, dm, edm, mdm, medm])

	
  
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
  dm = 99.99
  edm = 99.99
  mdm = 99.99
  medm = 99.99

  
  myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,Ks,logK,Vls, Vhelio, dcf2, ed, \
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_lum, \
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, tX_dyn, tX_lum, subGalaxies, nest, \
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname, dm, edm, mdm, medm])
  
  
  myTable.write(outfile, format='ascii.fixed_width',delimiter='|', bookend=False)
  
  # removing the last line, (it sits o adjust the column wodths)
  command =  ["csh", "remove_lastline.csh", outfile]
  subprocess.call(command)
################################################################# 
################################################################# 

def Theta_max(root):
  
  if root.left == None: return 0
  
  galList = NodeLeaf(root)
  
  N = len(galList)
  theta = np.zeros(N-1) 
  
  for i in range(1,N):
      
      theta[i-1] = angle(root.sgl, root.sgb, galList[i].sgl, galList[i].sgb)
  
  return np.max(theta)


################################################################# 
# gets a list of galaxies and returns Rg in terms of radian
def Rg_radian(galList):
  
  
  
  
  
  N = len(galList)
  sum = 0.
  BOL = True
  Rg = 0 
  for i in range(0,N-1):
    for j in range(i+1,N):
      
      distij = angle(galList[i].sgl, galList[i].sgb, galList[j].sgl, galList[j].sgb)
      if distij !=0 : 
         sum += 1. / distij
      else: BOL = False    # If one of the didtnaces is zero, therre is change to have duplicates 
      
  if BOL == True and sum != 0:
    Rg = ((N)*(N)) / sum
  
  return Rg




################################################################# 
################################################################# 

        
#################################################################
def mergGroup(G_list, restrict=False):
    

    
    NoGroups = len(G_list)
 
    i = 0
    while i < NoGroups:
      j = i+1
      while j < NoGroups-1:
	
        n1 = G_list[i][0].subGalaxies
        n2 = G_list[j][0].subGalaxies

	d1 = G_list[i][0].Vls/H0
	d2 = G_list[j][0].Vls/H0
	if d1 < 1 : d1 = 1 
	if d2 < 1 : d2 = 1
	r1 = (180.*atan(G_list[i][0].R_2t2/d1)/pi)
	r2 = (180.*atan(G_list[j][0].R_2t2/d2)/pi)


	ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0].gb, G_list[j][0].gl, G_list[j][0].gb)
	
	
	d1 = G_list[i][0].mDist
	e1 = d1 * G_list[i][0].mDistErr
	d2 = G_list[j][0].mDist
	e2 = d2 * G_list[j][0].mDistErr
	delt = abs(d1-d2)
	
	v1 = G_list[i][0].Vls
	v2 = G_list[j][0].Vls

	
	
	# using dynamical velocity dispersions calculated based on the luminosities
        sig1 = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
	sig2 = (G_list[j][0].M_v2 / (2.0E6))**(1./3)
	

	n1 = G_list[i][0].subGalaxies
	n2 = G_list[j][0].subGalaxies
	Bol = False
	

#######################################################################
        Sigquad = sqrt(sig1**2+sig2**2)
#######################################################################

	idd1 = G_list[i][0].nest
	idd2 = G_list[j][0].nest
        
        if restrict:

	  if (idd1==13418 and  idd2==14077):
	      Bol = True
	  if (idd1==36188 and  idd2==36136) or (idd1==36136 and  idd2==36188):
	      Bol = True
	  if (idd1==13324 and  idd2==13108) or (idd1==13108 and  idd2==13324):
	      Bol = True    
	  if (idd1==36699 and  idd2==36875) or (idd1==36875 and  idd2==36699):
	      Bol = True  
	      

	      
#######################################################################
      
	#if ang12 <= 1.1*(r1+r2) and max(r1,r2)<6:
	  #if abs(v1-v2) <= max(sig1,sig2) and (min(r1,r2))**3 < 0.2*(max(r1,r2))**3:
	       #Bol = True    
	       
	if ang12 <= 1.0*(r1+r2) and max(r1,r2)<6 and min(n1,n2)<5:
	  if abs(v1-v2) <= 2.0*max(sig1,sig2):
	      Bol = True
	      
	if ang12 <= 0.6*(r1+r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= 2*Sigquad:
	      Bol = True

	      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True
	
	      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt < min(e1, e2):
	      Bol = True
	      
	
	  
	# one group completely projected on another one (e.g. smal in big)
	if ang12 <= 1.0*(max(r1,r2)-min(r1,r2)) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True	      
#######################################################################

        if (idd1==39241 and  idd2==39600) or (idd1==39600 and  idd2==39241):
	  Bol = False
	if (idd1==13324 and  idd2==13108) or (idd1==13108 and  idd2==13324):
	  Bol = True    
	  

	if (idd1==12651 and  idd2==12651):
	  Bol = False
	if (idd1==13418 and  idd2==12651) or (idd1==12651 and  idd2==13418):
	  Bol = False  

	  
        if Bol:
	  
	  newGroup = []
	  
	  for p in range(1, len(G_list[i])):
	    newGroup.append(G_list[i][p])
	  
	  for q in range(1, len(G_list[j])):
	    newGroup.append(G_list[j][q])
	  
 
	  root = LoosgalListJoint(newGroup, grID = 200000000)
	  G_list[i] = NodeLeaf(root)

	  
	
	  G_list.pop(j) 
	  NoGroups-=1
	  i=0
	  j=0
        j+=1
      i+=1
    
    return G_list

 #################################################################
 # Function name: addGalGroup
 #
 # This function tries to find those galaxies which do not fall 
 # into any Group. If the galaxy is close enough to the center of any group,
 # according to its radial velocity (i.e. within 2*sigma of the group radial velocity) 
 # it might get absorbed by that group...
 #################################################################

def addGalGroup(G_list, galList):
  

    singles = []
    for galaxy in galList:
       if galaxy.inGroup == 0:
	  singles.append([galaxy, -1, -10000])  # [galaxy, new.group.id, angular.separation]
    
    N = len(G_list)
    #BB = []
    #LL = []
    #for i in range(N):
      #LL.append(G_list[i][0].gl)
      #BB.append(G_list[i][0].gb)
    
    for entity in singles:
        if  entity[0].id != 40001:
            gl = entity[0].gl
            gb = entity[0].gb
    
	    #q1 = np.zeros((N,), dtype=np.int)
            #q2 = np.zeros((N,), dtype=np.int)
            #q3 = np.zeros((N,), dtype=np.int)
            #q4 = np.zeros((N,), dtype=np.int)
	    #q1[np.where(LL<=gl+10)] = 1
	    #q2[np.where(LL>=gl-10)] = 1
	    #q3[np.where(BB<=gb+10)] = 1
	    #q4[np.where(BB>=gb-10)] = 1
	    #qq = q1+q2+q3+q4
	    #group_indices = np.where(qq==4)
	       
	    for i in range(len(G_list)): # group_indices[0]:
	      if  G_list[i][0].nest != 41220:  # Virgo
		    		    
		    ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0].gb, gl, gb)
		    d = G_list[i][0].Vls / H0
		    if d < 1: d = 1 
		    r = (180.*atan(G_list[i][0].R_2t2/(d))/pi)	      
		    
		    # if the galaxy is close enough to the center of the group
		    # and it's not already in any other group, then it would be added to the current group
		    sig_p = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
		    
		    
		    #if G_list[i][0].nest == 40692: 
		      #print "Esn TEst: ", G_list[i][0].Vls, sig_p, " <> ", entity[0].id, entity[0].Vls
		    
		    if ang12 <= 2*r:


		      if ang12 <= 1.01*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d) and abs(entity[0].Vls - G_list[i][0].Vls) < 2.5*sig_p:
			  d1 = G_list[i][0].mDist
			  e1 = d1 * G_list[i][0].mDistErr
			  d2 = entity[0].dcf2
			  e2 = d2 * entity[0].ed
			  delt = abs(d1-d2)
			  join = False
			  if ang12 > 1.01*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d) and d1!=0 and d2!=0 and delt<max(e1,e2):
			    join = True
			  if d1==0 or d2==0:
			    join = True
			  if ang12 <= 1.01*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d):
			    join = True
			    

			  if join==True and  sig_p > entity[2] and isAllowed(G_list[i][0].nest, entity[0].id):
			    entity[1] = i
			    entity[2] = sig_p
		  
            if entity[1] > -1: # if the galaxy is absorbed (check the id of the target group)
              entity[0].inGroup = 1
              G_list[entity[1]].append(entity[0])
    for i in range(N):
       G_list[i].pop(0)
       root  = LoosgalListJoint(G_list[i], grID = 300000000)
       G_list[i] = NodeLeaf(root) 
    return G_list


#################################################################
# Sometimes when groups grow the order of galaxy absorbation must be different:
def IndividualGal_modifier(G_list, galList):
  

  

  tmp_galList = []
  N = len(galList)
  for galaxy in galList:
    if galaxy.inGroup > 0:
      galaxy.inGroup = 0
      tmp_galList.append(galaxy)
    

  new_G_list=[]
  for Group in G_list:
    Group[1].inGroup = 1  # first galaxy of group, dominant one
    new_G_list.append([Group[0], Group[1]])
  

  
  
  Lum_mass = []
  for i in range(0, len(new_G_list)): Lum_mass.append(new_G_list[i][0].M_v2) # look at the mass of each group
  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)
  
  NEWgpLIST = []
  for i in indices[::-1]:
            
            
            grp = new_G_list[i][0]
            gl = grp.gl
            gb = grp.gb
    

	    j = 0 
	    N = len(tmp_galList)
	    while j < N: 
                     
                  gal = tmp_galList[j]

	          ang12 = (180./pi)*angle(gl, gb, gal.gl, gal.gb)
                  d = grp.Vls / H0
                  if d < 1: d = 1 
                  r = (180.*atan(0.95*grp.R_2t2/(d))/pi)	   
                  
                  if ang12 <= r and gal.inGroup == 0:
 
	              sig_p = (grp.M_v2 / (2.0E6))**(1./3)
	              if (abs(grp.Vls-gal.Vls) <= 2.0*sig_p):	
			new_G_list[i].append(gal)
			gal.inGroup = 1
		   	tmp_galList.pop(j)
			j -= 1
			N -= 1
	          j += 1

	    
    
	    new_G_list[i].pop(0)
	    if len(new_G_list[i]) > 1: 
	        grp =  LoosgalListJoint(new_G_list[i], grID = 400000000)
	        NEWgpLIST.append(NodeLeaf(grp))
	    else:
	        new_G_list[i][0].inGroup = 0
	        
  
  virg_w = None
  for group in NEWgpLIST:
    for galaxy in group[1:]:
      if galaxy.id == 39659:
	virg_w = group
        break
  
  if virg_w != None:
    newG = []
    #print "Here I am"
    for galaxy in virg_w[1:]:
      if galaxy.Vls < 1450 and galaxy.sgb > -7:
	galaxy.inGroup = 0
      else:
        newG.append(galaxy)
   
    grp = LoosgalListJoint(newG, grID = 400000000)
    virg_w = NodeLeaf(grp)
  
  
  
  # This is the output	    
  #G_list = NEWgpLIST
  return NEWgpLIST


#################################################################
## Trying to find linkages when both galaxies have distances

def find_pair_dist(galList_org):
  
  galList = []
  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 
	galList.append(galaxy)  
  
  individuals = []
  Lum_mass = []
  for i in range(0, len(galList)): Lum_mass.append(galList[i].M_v2)
  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)
  for i in indices[::-1]: individuals.append(galList[i])  
  
  NEWgpLIST =[]
  p = 0 
  N = len(individuals)
  while p < N-1:
    #print p, N
    galaxies = []
    galaxies.append(individuals[p])

    q = (p+1)%N
    pivot = p
    grp = galaxies[0]
    Bol = False
    
    if grp.dcf2==0:
      p+=1
      continue
    
    while q!=pivot:
      
      if individuals[q].dcf2==0:
	q = (q+1)%N
	continue
      
      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].gl, individuals[q].gb)
      
      d = grp.dcf2
      if d == 0:
	q = (q+1)%N
	continue
      
      
      mass = Mass(10**m_logK(grp.Ks, grp.Vls, d=d))
      R_2t2 = 0.215*((mass/1.E12)**(1./3))

      thet = (180.*atan(1.0*R_2t2/d)/pi)

      coeff = 1.3

      if abs(individuals[q].dcf2 - d)/d < 0.1 and ang12 < coeff*thet :
            
            mass = Mass(10**m_logK(individuals[q].Ks, individuals[q].Vls, d=individuals[q].dcf2))
            R_2t2_2 = 0.215*((mass/1.E12)**(1./3))
            
            Bol = True

	    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 400000000)
	    individuals.pop(q)

	    N-=1
	    
	    if p>=q:  p-=1
	    
	    
	    if len(galaxies) == 2:
	      individuals.pop(p)
	      N-=1	
	      if q>p: q = (q-1)%N
	      break
            
            
	    q = (q-1)%N
	    pivot = q 
	    
      q = (q+1)%N
    if Bol: 
      NEWgpLIST.append(NodeLeaf(grp))
    p+=1
    
    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1
  
  
  return NEWgpLIST  


#################################################################

def IndividualGal(galList_org, pairs=False):
  
  
  galList = []
  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 
	galList.append(galaxy)   
  

  individuals = []
  Lum_mass = []
  for i in range(0, len(galList)): Lum_mass.append(galList[i].M_v2)
  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)
  for i in indices[::-1]: individuals.append(galList[i])  
  
  NEWgpLIST =[]
  p = 0 
  N = len(individuals)
  while p < N-1:
    #print p, N
    galaxies = []
    galaxies.append(individuals[p])
    #individuals.pop(p)
    #N-=1
    q = (p+1)%N
    pivot = p
    grp = galaxies[0]
    Bol = False
    while q!=pivot:
      
      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].gl, individuals[q].gb)
      
      d = grp.Vls/H0
      if d <=1: d =1

      thet = (180.*atan(1.0*grp.R_2t2/d)/pi)
      #if grp.Vls<75:
	#thet = (180.*atan(1.0*grp.R_2t2/(1.))/pi)
      
      coeff = 1.0 # 1.05
      delta_v = abs(grp.Vls-individuals[q].Vls)
      
      
      test = False
      if pairs:
	sig_sig = 100000
	if len(galaxies) == 1:
	  d = individuals[q].Vls/H0
	  if d <=1: d =1
	  thet +=  (180.*atan(individuals[q].R_2t2/d)/pi)
	  L1 = 10**grp.logK
	  L2 = 10**individuals[q].logK
	  L_tot = L1 + L2
	  mass = Mass(L_tot)
	  sig_sig = (mass / (2.0E6))**(1./3)
	  if delta_v < 2.0*sig_sig: 
	    test = True
	
      
      
      if ang12 <= coeff*thet:
	sig_p = (grp.M_v2 / (2.0E6))**(1./3)

	if test == True or delta_v <= 2.0*sig_p:

	    Bol = True
	    

	    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 400000000)
	    individuals.pop(q)

	    N-=1
	    
	    if p>=q:  p-=1
	    
	    
	    if len(galaxies) == 2:
	      individuals.pop(p)
	      N-=1	
	      if q>p: q = (q-1)%N
            
            
	    q = (q-1)%N
	    pivot = q 
	    
      #if grp.id == 40001: Bol = False
      q = (q+1)%N
    if Bol: #len(galaxies) > 1:
      NEWgpLIST.append(NodeLeaf(grp))
    p+=1
    
    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1
  
  
  ##3 Dealingwith Virgo-W
  virg_w = None
  for group in NEWgpLIST:
    for galaxy in group[1:]:
      if galaxy.id == 39659:
	virg_w = group
        break
  
  if virg_w != None:
    newG = []
    #print "Here I am"
    for galaxy in virg_w[1:]:
      if galaxy.Vls < 1450 and galaxy.sgb > -7:
	galaxy.inGroup = 0
      else:
        newG.append(galaxy)
   
    grp = LoosgalListJoint(newG, grID = 400000000)
    virg_w = NodeLeaf(grp)
  
  
  return NEWgpLIST

#################################################################

def LoosgalListJoint(galList, grID = 500000000, noSort=False):
   
   if len(galList)==0: return None
   
   Lum_mass = []
   for i in range(0, len(galList)):
        Lum_mass.append(galList[i].Ks)
   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass)
   
   if noSort:
      Lum_mass = []
      for i in range(0, len(galList)):
	    Lum_mass.append(galList[i].logK)
      Lum_mass = np.asarray(Lum_mass)
      indices = np.argsort(Lum_mass)
      indices = indices[::-1]
   
   root = None
   for i in indices:
       root = galJoint(root, galList[i], grID + galList[i].id)
   
   
   sumDist = 0
   sumError = 0
   V_mean = 0 
   for galaxy in galList:
    sumDist += galaxy.sumDist
    sumError += galaxy.sumError
    V_mean += galaxy.Vls
   
   V_mean /= len(galList)
   mDM = 0.
   meDM = 0.
   mDist = 0.
   mDistErr = 0.
   if sumDist !=0 and sumError!=0:
      mDM = sumDist/sumError
      meDM = sqrt(1/sumError)
      mDist = 10**(mDM/5.-5.)
      mDistErr = (0.2*log(10.))*meDM

   
   
   for gal in galList:
      gal.logK = m_logK(gal.Ks, V_mean)
      gal.mDist = mDist
      gal.mDistErr = mDistErr
      gal.mDM = mDM
      gal.meDM = meDM
     
   root = None
   for i in indices:
      root = galJoint(root, galList[i], grID + galList[i].id)   
   
   
   root.dcf2 = 0
   root.ed = 0
   
   return root

#################################################################

def Theta_gr(head, galList):
  
  
  
  N = len(galList)
  theta = np.zeros(N-1) 
  
  for i in range(1,N):
      
      theta[i-1] = angle(head.l, head.b, galList[i].l, galList[i].b)
  
  return np.max(theta)



#################################################################

def galList_moderator(galList):
  
  for galaxy in galList:
    if galaxy.inGroup <= 0:
      galaxy.mDist = galaxy.dcf2Copy
      galaxy.mDistErr = galaxy.edCopy
      galaxy.mDM = galaxy.dmCopy
      galaxy.meDM = galaxy.edmCopy
      
      if galaxy.id in [42666]:
	galaxy.mDist = 16.
	galaxy.mDistErr = 0.2
	galaxy.mDM  = (log10(16.)+5.)*5.
	galaxy.meDM = 0.2 / (0.2*log(10.))
	meanDist = 16.
	meanDistErr = 0.2
	galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = galaxy.Vls, dist=galaxy.dcf2)
	
      
      
     
      #if galaxy.sgl>110 and galaxy.sgl<125 and galaxy.sgb>-8 and galaxy.sgb<2:
	#if  galaxy.Vls<=1750: 
           #galaxy.setMeanDist(galaxy.mDist, galaxy.mDistErr, GRP_vel = galaxy.Vls, dist=16.)
	#else:
	   #galaxy.setMeanDist(galaxy.mDist, galaxy.mDistErr, GRP_vel = galaxy.Vls, dist=32)

#################################################################
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
#################################################################

def group_moderator(G_list, use_dist=False):

  for group in G_list:
	  groupHeader = group[0]
	  groupHeader.flag = 2
	  groupHeader.dcf2 = 0
	  groupHeader.ed = 0
	  groupHeader.dcf2Copy = 0
          groupHeader.edCopy = 0
          groupHeader.dm = 0
          groupHeader.edm = 0
          groupHeader.dmCopy = 0
          groupHeader.edmCopy = 0
          
          

          sum_v = 0.
          sum_v2 = 0.
          helio_sum_v = 0.
          helio_sum_v2 = 0.          
          n_v = 0
          for galaxy in group[1:]:

	    if galaxy.Vhelio == 0:
	      galaxy.Vls = 0
	    if galaxy.Vls !=0:
	      sum_v += galaxy.Vls
	      sum_v2 += galaxy.Vls**2
	      helio_sum_v += galaxy.Vhelio
	      helio_sum_v2 += galaxy.Vhelio**2
	      n_v += 1
	  
	  if n_v != 0:
	    v_av = sum_v / n_v
	    helio_v_av = helio_sum_v / n_v
	    
	    groupHeader.Vls = v_av
	    groupHeader.sigma = sqrt(sum_v2/n_v - v_av**2) 
	    
	    groupHeader.Vhelio = helio_v_av

	  else:
	    groupHeader.Vls = 0
	    groupHeader.sigma = 0  # in Vls
	    groupHeader.Vhelio = 0
	  
	  
	  mDM = 0.
	  meDM = 0.
	  meanDist = 0.
	  meanDistErr = 0.
	  sumDist = 0
          sumError = 0
	  for galaxy in group[1:]:
	    if galaxy.dcf2 != 0 and galaxy.ed != 0:
	      if not galaxy.id in [43775, 42447, 13368]:
		    err =  galaxy.dcf2*galaxy.ed
		    sumDist += galaxy.dcf2/(err**2)
		    sumError += 1./(err**2)
	  
	  
	  if sumDist !=0 and sumError != 0:
            meanDist = sumDist/sumError
            meanDistErr = sqrt(1./sumError)/meanDist
            
           
              
          if meanDist>0:
            groupHeader.mDM  = (log10(meanDist)+5.)*5.
            groupHeader.meDM = meanDistErr / (0.2*log(10.))
          else:
	    groupHeader.mDM = 0
	    groupHeader.meDM = 0
	  
	  groupHeader.dm  = groupHeader.mDM
          groupHeader.edm = groupHeader.meDM
          groupHeader.mDist = meanDist
          groupHeader.mDistErr = meanDistErr
          
          if groupHeader.nest in [40229, 41101, 41893, 43798, 44797]:
	    groupHeader.dm = 16.
	    groupHeader.edm = 0.2
	    groupHeader.mDist = 16.
	    groupHeader.mDistErr = 0.2
	    groupHeader.mDM  = (log10(16.)+5.)*5.
	    groupHeader.meDM = 0.2 / (0.2*log(10.))
	    meanDist = 16.
	    meanDistErr = 0.2

	    
 
	  L_tot = 0
	  ID = groupHeader.id 
	  groupHeader.id = ID - ID % 100000000 + groupHeader.nest
	  
	  V_hel_sum = 0.
	  
	  for galaxy in group[1:]:
	    galaxy.flag = 1
	    # and also modify the absolute luminosities
	    if use_dist:
	       if galaxy.Vls<200: #galaxy.nest == 5064336 or galaxy.nest == 2557:  # MilkyWay or Andromeda
	         galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls, dist=galaxy.dcf2)
	         A= True
	         
	       #elif groupHeader.sgl>110 and groupHeader.sgl<125 and groupHeader.sgb>-8 and groupHeader.sgb<2:
		   #if groupHeader.mDist>0:
		     #galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls, dist=galaxy.mDist)
		   #else:
		     #if  groupHeader.Vls<=1750: 
		       #galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls, dist=16.)
		     #else:
		       #galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls, dist=32)
	       else: 
		 galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls)
	    else:
	       galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls)
	   
	    galaxy.mDM = mDM
	    galaxy.meDM = meDM
	    L_tot += 10**galaxy.logK
	    galaxy.mDM  = groupHeader.mDM
	    galaxy.meDM = groupHeader.meDM
	    
	  
	  
          groupHeader.logK = log10(L_tot)
          Dist_v =  groupHeader.Vls / H0
          if Dist_v<1. : Dist_v=1
          Mk_sun = 3.28   # K-band
          M = Mk_sun - (groupHeader.logK / 0.4)
          groupHeader.Ks = M + 5*log10(Dist_v) + 30 - 5
          
          groupHeader.M_v2 = Mass(10**groupHeader.logK)
          if np.isinf(groupHeader.M_v2): 
	    print 'Modification ...'
	    print Mass(10**10.6453)
	    print 10**groupHeader.logK
	    sys.exit()
          groupHeader.R_2t2 = 0.215*((groupHeader.M_v2/1.E12)**(1./3))  # Mpc
  
  return G_list

#################################################################
#################################################################
#################################################################
# removing a galaxy from all group
def removeFrom(G_list, id):
  
  for i in range(len(G_list)):
    group = G_list[i]
    for galaxy in group[1:]:
      if galaxy.id == id:
        newG = []
        for galaxy in group[1:]:
           if galaxy.id == id:
        	galaxy.inGroup = 0
           else:
                newG.append(galaxy)
        
	if len(newG) > 1 : 
          grp = LoosgalListJoint(newG, grID = 500000000)
          group = NodeLeaf(grp)
          G_list.pop(i)
          G_list.append(group)
        elif  len(newG) == 1 :
	  newG[0].inGroup = 0
	  G_list.pop(i)
	return
#################################################################
def destroy_group(G_list, gal_id):
  
  j = -1
  N_grp = len(G_list)
  i = 0
  while i < N_grp:
    group = G_list[i]
    for galaxy in group[1:]:
      if galaxy.id == gal_id:
	j = i
	i = N_grp+1
	break
    i+=1
  
  if j>=0:
    for galaxy in G_list[j][1:]:
      galaxy.inGroup = 0
    
    G_list.pop(j)
  
  return
#################################################################
# removing a galaxy from a group
def removeFrom2(G_list, gr, gal_id):

  id = None
  for i in range(len(G_list)):
    group = G_list[i]
    if group[0].nest == gr:
      id = i
      break
  
  if id!= None:
    group = G_list[id]
    for galaxy in group[1:]:
      if galaxy.id == gal_id:
	  newG = []
	  for galaxy in group[1:]:
	    if galaxy.id == gal_id:
		  galaxy.inGroup = 0
	    else:
		  newG.append(galaxy)
	  
	  if len(newG) > 1 : 
	    grp = LoosgalListJoint(newG, grID = 500000000)
	    group = NodeLeaf(grp)
	    G_list.pop(i)
	    G_list.append(group)
	  elif  len(newG) == 1 :
	    newG[0].inGroup = 0
	    G_list.pop(i)
	    
	  return
#################################################################


##############################################################
# removing a galaxy from a group
def trimGroup_vel(G_list, gr, Vmin=-100000, Vmax=100000):
  
  id = None
  for i in range(len(G_list)):
     group = G_list[i]
     if group[0].nest == gr:
       id = i
       break
  
  trimed_gal = []
  if id != None:
     newGroup = []
     for p in range(1, len(G_list[id])):
        galaxy = G_list[id][p]
        if galaxy.Vls < Vmax and galaxy.Vls > Vmin:
	  newGroup.append(galaxy)
	else:
	  galaxy.inGroup = 0
	  trimed_gal.append(galaxy)
 
     root = LoosgalListJoint(newGroup, grID = 700000000) 
     group = NodeLeaf(root)   
     G_list.pop(id)
     G_list.append(group)
     return trimed_gal
  else:
    return None
#################################################################
# removing a galaxy from a group
# rad (radius) --> in degrees
def trimGroup_rad(G_list, gr, factor=1, rad=None):
  
  id = None
  for i in range(len(G_list)):
     group = G_list[i]
     if group[0].nest == gr:
       id = i
       break
     
  if id != None:
     newGroup = []
     for p in range(1, len(G_list[id])):
        galaxy = G_list[id][p]
        ang12 = (180./pi)*angle(G_list[id][0].gl, G_list[id][0].gb, galaxy.gl, galaxy.gb)
        d = G_list[i][0].Vls / H0
        if d < 1: d = 1 
        r = (180.*atan(G_list[i][0].R_2t2/(d))/pi)
        
        if rad != None:
	  radius = rad*pi/180.  # --> radius in radian
	else: 
	  radius = factor*G_list[id][0].R_2t2/d
        
        if ang12 <= 1.05*(180./pi)*atan(radius):
	  newGroup.append(galaxy)
	else:
	  galaxy.inGroup = 0
 
     root = LoosgalListJoint(newGroup, grID = 700000000) 
     group = NodeLeaf(root)   
     G_list.pop(id)
     G_list.append(group)   
        

#################################################################
# removing a galaxy from a group
def manualBlend(G_list, gr1, gr2):
  
   id1 = None; id2 = None
   
   for i in range(len(G_list)):
     group = G_list[i]
     if group[0].nest == gr1:
       id1 = i
       break
   for j in range(len(G_list)):
     group = G_list[j]
     if group[0].nest == gr2:
       id2 = j
       break   
   
   if id1 != None and id2 != None:
     newGroup = []
     for p in range(1, len(G_list[id1])):
	newGroup.append(G_list[id1][p])
     for q in range(1, len(G_list[id2])):
	newGroup.append(G_list[id2][q])   
     
     root = LoosgalListJoint(newGroup, grID = 600000000) 
     group = NodeLeaf(root)   
     i = max(id1, id2); j = min(id1, id2)
     G_list.pop(i); G_list.pop(j)
     G_list.append(group)
   
    
#################################################################
# Note that dcf2 and ed would remain the same
# Bu their effect is discarded
#def distanceDiscard(G_list, id):
  
  #for i in range(len(G_list)):
    #group = G_list[i]
    #for galaxy in group[1:]:
      #if galaxy.id == id:
        #newG = []
        #for galaxy in group[1:]:
           #if galaxy.id == id:
        	#galaxy.sumDist = 0
        	#galaxy.sumError = 0 
           #newG.append(galaxy)
        
        #grp = LoosgalListJoint(newG, grID = 500000000)
        #group = NodeLeaf(grp)
        #G_list.pop(i)
        #G_list.append(group)
	#return  

#################################################################
def addTo(G_list, galList, gr_id, gal_id):
  
  removeFrom(G_list, gal_id)
  
  id = None
  B = False
  for i in range(len(G_list)):
     group = G_list[i]
     if group[0].nest == gr_id:
       id = i
       break
     else:
       for gal in group[1:]:
	 if gal.id ==  gr_id:
	   id = i
	   B = True
	   break
       if B: break   
     
  idd = None
  for p in range(len(galList)):
     galaxy = galList[p]
     if galaxy.id == gal_id:
       idd = p
       break
     
  
  if id != None and idd != None:
     newGroup = G_list[id][1:]
     newGroup.append(galList[p])
     galList[p].inGroup = 1

 
     root = LoosgalListJoint(newGroup, grID = 700000000, noSort=True) 
     group = NodeLeaf(root)   
     G_list.pop(id)
     G_list.append(group)       
#################################################################
def formGroup(G_list, galList, galaxies):
  
  newGroup=[]
  for gal_id in galaxies:
    brake = False
    for group in G_list:
      for p in range(1, len(group)):
	if group[p].id == gal_id:
	  brake = True
	  break
      if brake: break
    if brake == False:
      for galaxy in galList:
        if galaxy.id == gal_id:
	    if galaxy.inGroup == 1: 
	      break	  
	    else: 
	      newGroup.append(galaxy)
	      galaxy.inGroup = 1
	      break

  if len(newGroup) > 1:
    root = LoosgalListJoint(newGroup, grID = 700000000) 
    group = NodeLeaf(root)   
    G_list.append(group)  
  

#################################################################

def brent_modifications(G_list, galList):
  
  ### South
  removeFrom(G_list, 12662)  # gal_id
  removeFrom(G_list, 2578)   # gal_id
  manualBlend(G_list, 62836, 62453) # group_id1, group_id2
  trimGroup_vel(G_list, 13255,  Vmax=1600) # gr_id1, remove those with Vls>1600
  #distanceDiscard(G_list, 13368) # gal_id
  
  removeFrom2(G_list, 2557, 4126)
  removeFrom2(G_list, 2789, 2578)
  addTo(G_list, galList, 13505, 13470)
  
  removeFrom2(G_list, 70090, 70069)
  removeFrom2(G_list, 70090, 70184)
  From_To(G_list, galList, 70080, 70184, 70069)
  
  addTo(G_list, galList, 70069, 70306)
  addTo(G_list, galList, 70069, 135466)
    
  addTo(G_list, galList, 13434, 13154)
  
  addTo(G_list, galList, 13419, 13342)
  removeFrom(G_list, 13368)
  addTo(G_list, galList, 13419, 13304)
  
  
  ### North
  trimGroup_vel(G_list, 34695,  Vmin=400, Vmax=1000)
  trimGroup_rad(G_list, 34695,  factor = 1.05)
  
  addTo(G_list, galList, 34426, 86673)
    
  trimGroup_vel(G_list, 32256,  Vmax=800) # 
  
  trimGroup_vel(G_list, 39600,  Vmax=700) # 
  addTo(G_list, galList, 39241, 39344) # gr_id, gal_id
  addTo(G_list, galList, 39241, 40228)
  addTo(G_list, galList, 39600, 40537)
  addTo(G_list, galList, 39600, 39191)
  
  addTo(G_list, galList, 39600, 166131)
  addTo(G_list, galList, 39600, 2832112)
  addTo(G_list, galList, 39600, 5057019)
  addTo(G_list, galList, 39600, 5057020)
  addTo(G_list, galList, 39600, 5061328)
  
  addTo(G_list, galList, 39600, 2832111)
  removeFrom2(G_list, 39600, 40665)
  
  addTo(G_list, galList, 39241, 39237)
  addTo(G_list, galList, 39241, 39864)
  addTo(G_list, galList, 39241, 166129)

  
  trimGroup_vel(G_list, 38440,  Vmin=700) # 
  addTo(G_list, galList, 38440, 38068)
  
  manualBlend(G_list, 28630, 28655)
  
  removeFrom(G_list, 24213) 
  addTo(G_list, galList, 23324, 24213)
  
  manualBlend(G_list, 21396, 21102)
  
  manualBlend(G_list, 46957, 45279)
  removeFrom2(G_list, 46957, 47847)
  removeFrom2(G_list, 46957, 166152)
  removeFrom2(G_list, 46957, 592761)
  addTo(G_list, galList, 46957, 46680)
  addTo(G_list, galList, 46957, 166158)
  addTo(G_list, galList, 46957, 166167)
  addTo(G_list, galList, 46957, 166172)
  addTo(G_list, galList, 46957, 166175)
  addTo(G_list, galList, 46957, 2815820)
  addTo(G_list, galList, 46957, 2815822)
  addTo(G_list, galList, 46957, 2815823)
  addTo(G_list, galList, 46957, 4689187)
  
  removeFrom2(G_list, 48082, 48334)
  
  
  
  trimGroup_vel(G_list, 43495,  Vmax=500)
  trimGroup_rad(G_list, 43495,  factor = 1.1)
  removeFrom2(G_list, 40973, 40904)
  addTo(G_list, galList, 43495, 40973)
  addTo(G_list, galList, 43495, 45314)
  addTo(G_list, galList, 43495, 166146)
  
  addTo(G_list, galList, 42575, 166140)
  addTo(G_list, galList, 42575, 42045)
  removeFrom(G_list, 41902)
  addTo(G_list, galList, 42575, 41902)
  
  addTo(G_list, galList, 39225, 38881)
  
  manualBlend(G_list, 39225, 39023)
  
  
  trimGroup_vel(G_list, 47404,  Vmin=400)
  
  formGroup(G_list, galList, [46039, 45939, 45506, 46127, 5057032])
  
  removeFrom2(G_list, 29265, 29033)
  
  manualBlend(G_list, 13826, 15345)
  
  manualBlend(G_list, 25950, 26259)
  
  addTo(G_list, galList, 24930, 166096)
  

  
  addTo(G_list, galList, 39724, 40692)
  addTo(G_list, galList, 40692, 38742)
  addTo(G_list, galList, 40692, 39690)
  addTo(G_list, galList, 40692, 38598)
  addTo(G_list, galList, 40692, 40495)
  addTo(G_list, galList, 40692, 4326021)
  #addTo(G_list, galList, 40692, )
  
  From_To(G_list, galList, 30083, 30445, 30493)

  From_To(G_list, galList, 30744, 31166, 31029)
  
  removeFrom2(G_list, 39225, 38685)
  
  removeFrom2(G_list, 41789, 1181814)
  #distanceDiscard(G_list, 43775)
  #distanceDiscard(G_list, 42447)
  
  addTo(G_list, galList, 38440, 38356)
  
  ### NGC784 Group (7671)
  #addTo(G_list, galList, 7671, 166064)
  #addTo(G_list, galList, 7671, 6699)
  #addTo(G_list, galList, 7671, 8484)
  
  #### Cen A Group (46957)
  addTo(G_list, galList, 46957, 44110)
  
  ### NGC1313 Group (12286)
  addTo(G_list, galList, 12286, 166073)
  
  ### M83 Group (48082)
  addTo(G_list, galList, 48082, 166170)
  addTo(G_list, galList, 48082, 166176)
  
  ### M81 Group (28630)
  addTo(G_list, galList, 28630, 29231)
  addTo(G_list, galList, 28630, 31286)
  addTo(G_list, galList, 28630, 166101)
  addTo(G_list, galList, 28630, 5056932)
  addTo(G_list, galList, 28630, 5056942)
  addTo(G_list, galList, 28630, 5056933)
  addTo(G_list, galList, 28630, 5056947)
  addTo(G_list, galList, 28630, 5056938)
  addTo(G_list, galList, 28630, 5056943)
  addTo(G_list, galList, 28630, 5056937)
  addTo(G_list, galList, 28630, 5056934)
  addTo(G_list, galList, 28630, 5056931)
  addTo(G_list, galList, 28630, 5056944)
  addTo(G_list, galList, 28630, 2807133)
  addTo(G_list, galList, 28630, 3097828)
  addTo(G_list, galList, 28630, 5056941)  
  
  
  addTo(G_list, galList, 32256, 4689201)  
  addTo(G_list, galList, 32256, 4689211) 
  addTo(G_list, galList, 32256, 4689206)  
  addTo(G_list, galList, 32256, 4689203)   
  addTo(G_list, galList, 32256, 4689205)  
  addTo(G_list, galList, 32256, 4689202) 
  addTo(G_list, galList, 32256, 4689197)  
  addTo(G_list, galList, 32256, 83355)  
  addTo(G_list, galList, 32256, 4689213)  
  addTo(G_list, galList, 32256, 4689217) 
  addTo(G_list, galList, 32256, 83321)  
  addTo(G_list, galList, 32256, 2806988)   
  addTo(G_list, galList, 32256, 166107)  


  addTo(G_list, galList, 46957, 2815821)  
  addTo(G_list, galList, 46957, 2815819)  
  addTo(G_list, galList, 46957, 2815824)  
  addTo(G_list, galList, 46957, 166164)   
  addTo(G_list, galList, 46957, 3097729)  
  
  addTo(G_list, galList, 48082, 3097728)   
  
  addTo(G_list, galList, 40692, 166133)  
  
  addTo(G_list, galList, 46153, 166161)    
  
  addTo(G_list, galList, 33550, 3097701) 
  
  addTo(G_list, galList, 62836, 2815833) 
  addTo(G_list, galList, 62836, 2815832) 
  addTo(G_list, galList, 62836, 62815) 
  
  
  addTo(G_list, galList, 46039, 166162) 
  
  addTo(G_list, galList, 42637, 166145)   
  addTo(G_list, galList, 42637, 5056992)     
  
  addTo(G_list, galList, 39600, 166127)       
  
  addTo(G_list, galList, 42407, 3097710)   
  addTo(G_list, galList, 42407, 3097711)     
  addTo(G_list, galList, 42407, 3097712)     
  addTo(G_list, galList, 42407, 3097709)   
  
  addTo(G_list, galList, 50063, 49404)    
  
  addTo(G_list, galList, 13826, 2807114)      
  
  addTo(G_list, galList, 9332, 83333)    
  
  addTo(G_list, galList, 51233, 2801016)    
  
  addTo(G_list, galList, 32256, 83338)  
  
  addTo(G_list, galList, 29265, 3097699)    
  
  addTo(G_list, galList, 25950, 166098)  
  addTo(G_list, galList, 25950, 166097)    
  addTo(G_list, galList, 25950, 166099)    
  
  addTo(G_list, galList, 33550, 5057225)     
  addTo(G_list, galList, 32256, 87258)    
  addTo(G_list, galList, 47404, 47658)     
  addTo(G_list, galList, 25950, 4689198)  
  
  formGroup(G_list, galList, [30087, 5056988])     
  formGroup(G_list, galList, [40596, 2807147])    
  formGroup(G_list, galList, [50961, 2801017])     

  

  
  
  ### Hydra group - V_ls < 6000
  addTo(G_list, galList, 31478, 31599)
  addTo(G_list, galList, 31478, 31951)

  ### Centaurus group - V_ls < 6000  
  addTo(G_list, galList, 43296, 43466)
  
  manual_group(G_list, galList, [39225, 39023, 38881, 39145])
  
  destroy_group(G_list, 63616)
  destroy_group(G_list, 73049)
  
  
  # Based on the file: new_nearby_galaxies.20170220
  addTo(G_list, galList, 29128, 5098252)
  addTo(G_list, galList, 21396, 6656984)
  addTo(G_list, galList, 41220, 5058305)  #  (Virgo)
  addTo(G_list, galList, 50063, 5067385)
  addTo(G_list, galList, 50063, 5067386)
  addTo(G_list, galList, 50063, 5067387)

  # isolated ones:
  removeFrom(G_list, 5072710)
  removeFrom(G_list, 5072714)
  removeFrom(G_list, 5072715)
  
  
#################################################################

def From_To(G_list, galList, gr1, gr2, gal_id):
  
  removeFrom2(G_list, gr1, gal_id)
  addTo(G_list, galList, gr2, gal_id)
  

def isAllowed(group_id, galaxy_id):
  
  allowed = True
  
  ## not allowed list
  if group_id == 34695 and galaxy_id in [3471336, 4098734, 2807141]:
    allowed = False
  
  return allowed
  

#################################################################
def manual_group(G_list, galList, id_list, verbose=False):  # my_group: id list of members
  
  
    my_group = []
    for galaxy in galList: 
      id = galaxy.id
      if id in id_list: 
	removeFrom(G_list, id)
	my_group.append(galaxy)
	galaxy.inGroup = 1
	galaxy.logK = m_logK(galaxy.Ks, galaxy.Vls, d=galaxy.dcf2)
	
    if len(my_group)>0: 
      root = LoosgalListJoint(my_group, grID = 900000000)
      new_Group = NodeLeaf(root)
      G_list.append(new_Group)  
    else: return None
    
    if verbose:
      for all in new_Group:
	print all.id
      

    return my_group
  

#################################################################

if __name__ == '__main__':
  
  
  # Originally made by M_L_ratio_curve_v2.py
  table = np.genfromtxt('M_L_curve_v2.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  Lumin_in   = 10**table['log10_L']
  Mass_out   = 10**table['log10_M']
  
  LtoM_func = interpolate.interp1d(Lumin_in, Mass_out)
  

  R_max=50. # R_max=50.
  cluster = ''

  
  if len(sys.argv) < 2:
    print "\nEnter the cluster name as the input ..." 
    print "\nexample: \n > python " + str(sys.argv[0])+ " manual"
    print "\nPossible options: all, manual, local" 
    print "Use north/south for whole sky.\n"
    sys.exit(1)
  
  cluster = str(sys.argv[1])
  
  
  inFile = 'AllSky.v24.cf3.csv'
  
  
  table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
  
  


  ignore_list = [40001, 42447, 40240, 40851, 3085, 70090]   # April 18, 2016 (added pgc3085 SMC)
  Vh_zero = [29231,  31286,  40747,  46680, 166101, 166146, 166158, 166172, 166175, 166176]
  
 
  ignore_list = np.asarray(ignore_list)
  Vh_zero = np.asarray(Vh_zero)
  ignore_list = np.concatenate((ignore_list, Vh_zero))
  
  # added May 3, 2016
  Vh_zero = [166167, 2815820, 2815822, 2815823, 4689187]
  Vh_zero = np.asarray(Vh_zero)
  ignore_list = np.concatenate((ignore_list, Vh_zero))
  


  # <> ---------- Local groups
  MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  id_MW = MW_g['PGC']
  ignore_list = np.concatenate((ignore_list, id_MW))
  
  M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=',', filling_values=0, names=True, dtype=None)    
  id_M31 = M31_g['PGC']
  ignore_list = np.concatenate((ignore_list, id_M31))
  # </> ---------- Local groups
  

  
  if cluster == 'all' :
     galList = readgalList(table, ignore_list=ignore_list, Vls_max=4000)
     print 'No of galaxies with Vls < 4000: ', len(galList)
     sys.exit()
  elif cluster == 'manual' :
     alpha =  102.8806
     delta = -2.3479
     step = 15
     galList = readgalList(table, skyPatch='manual', sgl_min=alpha-step, sgl_max=alpha+step, sgb_min=delta-step, sgb_max=delta+step, ignore_list=ignore_list, Vls_max=4000)    
  elif cluster == 'local':
     galList = readgalList(table, skyPatch='local', Vls_max=200)
  else:
    sys.exit(1)

  
  ### Virgo: Ignoring galaxies inside 6.8 deg
  Virgo_ignore = []
  for galaxy in galList:
    ang12 = (180./pi)*angle(galaxy.sgl, galaxy.sgb, sglV, sgbV)
    if ang12<= 6.8:
      Virgo_ignore.append(galaxy.id)
  
  if len(Virgo_ignore)>0:
    Virgo_ignore = np.asarray(Virgo_ignore)
    ignore_list = np.concatenate((ignore_list, Virgo_ignore))
  
  Virgo_W1 = [40375, 40122, 40886, 40516, 40033, 40119,41050]
  Virgo_W1 = np.asarray(Virgo_W1)
  
  Virgo_W2 = [42741, 42619, 43001]
  Virgo_W2 = np.asarray(Virgo_W2)
  
  Virgo_M   = [38890, 38916, 39002, 39025, 39040, 39152, 39390, 38749, 38792, 38885, 39124]
  Virgo_M   = np.asarray(Virgo_M)
  

  Virgo_W_g = np.genfromtxt('brent_virgoW.csv' , delimiter=',', filling_values=0, names=True, dtype=None)   
  id_Virgo_w = Virgo_W_g['PGC'][0:41]
  ignore_list = np.concatenate((ignore_list, id_Virgo_w))


  Virgo_foreback = [40045, 42081, 43072, 44491, 43601, 43254, 42833]
  Virgo_foreback = np.asarray(Virgo_foreback)
  Virgo_foreback = np.concatenate((Virgo_foreback, Virgo_W_g['PGC'][41:50]))
  #ignore_list = np.concatenate((ignore_list, Virgo_foreback))
  
  virgo_g1 = [39620,40004,40087,40240,40339,40469,40494,40801,1279440,1293380]
  virgo_g1 = np.asarray(virgo_g1)
  ignore_list = np.concatenate((ignore_list, virgo_g1))
  
  virgo_g2 = [39483,39537,39646,39809,39878,40109,40229,40321,3394197]
  virgo_g2 = np.asarray(virgo_g2)
  ignore_list = np.concatenate((ignore_list, virgo_g2))
  
# ===================================================================== 
  #id_all   = table['pgc']
  #Ks_all   = table['Ks'] 
  
  
  #MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  #id = MW_g['PGC']
  #Vmag = MW_g['Vmag']
  #Ks0=[]
  #Vm0=[]
  #for i in range(len(id)):
    #for j in range(len(id_all)):
      #if id_all[j] == id[i] and Vmag[i]!=Ks_all[j] and Vmag[i]>10:
	#Ks0.append(Vmag[i]-Ks_all[j])
	#Vm0.append(Vmag[i])
	#break
    
  #plt.plot(Vm0, Ks0, '+')
  
  
  
  #M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  #id = M31_g['PGC']
  #Vmag = M31_g['Vmag']
  #Ks=[]
  #Vm=[]
  #for i in range(len(id)):
    #for j in range(len(id_all)):
      #if id_all[j] == id[i] and Vmag[i]!=Ks_all[j] and Vmag[i]>10:
	#Ks.append(Vmag[i]-Ks_all[j])
	#Vm.append(Vmag[i])
	#break
    
  #plt.plot(Vm, Ks, '*', color='red')  

  #Ks = np.asarray(Ks)
  #Vm = np.asarray(Vm)
  #Ks0 = np.asarray(Ks0)
  #Vm0 = np.asarray(Vm0)  
  #Ks = np.concatenate((Ks0, Ks))
  #Vm = np.concatenate((Vm0, Vm))
  #sig = np.std(Ks)
  #med = np.median(Ks)
  #print med, sig
  #flag = np.zeros(len(Ks))
  #flag[np.where(Ks>med+1.5*sig)] = 1.
  #flag[np.where(Ks<med-1.5*sig)] = 1.
  #ind = np.where(flag == 0)
  #Ks = Ks[ind]
  #Vm = Vm[ind]
  #med = np.median(Ks)
  #plt.plot([10, 20],[med, med], '--')
  #print med
  
  
  #plt.xlabel('Vmag')
  #plt.ylabel('Vmag-Ks')
  #plt.show()
  #sys.exit(1)
# ===================================================================== 
  

  for iter in range(4):
   
        for gal in galList:
	  if gal.id in ignore_list:
	    gal.inGroup = -1
	    #print gal.id

        # these are then added manually
        # Their Vh=0 but the belong to gr39600
        for gal in galList:
	  if gal.id in [166131, 2832112, 5057019, 5057020, 5061328]:
	    gal.inGroup = -1
    
        print "working on grouping ... iter:", iter
	##  Finding groups based on the criteria
        G_list = []
        

        
        
        for qp in range(0,3):
	   G_list += IndividualGal(galList)
	   G_list += find_pair_dist(galList)
	
	## fixed
        G_list = addGalGroup(G_list, galList)
        
        group_moderator(G_list)
        
        inject_list=[40001, 40240, 40851, 42447, 70090]   
        for gal in galList:
	  if gal.id in inject_list:
	    gal.inGroup = 0
        
        
        if True:
          for pq in range(0,10): 
	    G_list = IndividualGal_modifier(G_list, galList)
	    G_list += find_pair_dist(galList)
	     
	  for pq in range(0,10): 
	    G_list = addGalGroup(G_list, galList)
	     
          for qp in range(0,10):
            mergGroup(G_list)
          
          group_moderator(G_list)
          mergGroup(G_list, restrict=True)
           
	  for pq in range(0,10): 
	    G_list = addGalGroup(G_list, galList)  
	     
          for qp in range(0,2):
            mergGroup(G_list)
             
          mergGroup(G_list, restrict=True)    
           
          group_moderator(G_list)
          


####################################################################

	##Error: galList is modified ?
	for galaxy in galList:
	  if galaxy.inGroup == 0:
	    galaxy.mDist = 0
	    galaxy.setMeanDist(galaxy.dcf2, galaxy.ed, 0)
	    

	brent_modifications(G_list, galList)  
	group_moderator(G_list)
	
	G_list += IndividualGal(galList, pairs=True)
	G_list = addGalGroup(G_list, galList)
	brent_modifications(G_list, galList) 
	group_moderator(G_list)
	G_list = addGalGroup(G_list, galList)
	
	for qp in range(0,2):
	  mergGroup(G_list, restrict=True)    
        group_moderator(G_list)
        
        brent_modifications(G_list, galList) 





        if True: # cluster=='manual':

	    destroy_group(G_list, 40001)
	    
	    for galaxy in galList: 
	      if galaxy.inGroup == 1 and galaxy.sgl>106 and galaxy.sgl<110 and galaxy.sgb<-5.5 and galaxy.sgb>-8:
		destroy_group(G_list, galaxy.id)
	    
	        
	    
	    ##manual_group(G_list, galList, virgo_g1)
	    g_lst = []
	    for galaxy in galList:
	      if galaxy.id in virgo_g1:
		galaxy.inGroup = 0
		g_lst.append(galaxy)
		galaxy.setMeanDist(16, 0.1, GRP_vel = galaxy.Vls, dist=16.)
	    
	    new_grp = IndividualGal(g_lst)
	    G_list += new_grp
	    meanDist    = new_grp[0][0].mDist
	    meanDistErr = new_grp[0][0].mDistErr
	    for galaxy in g_lst:
	      if galaxy.inGroup == 0: 
		#galaxy.inGroup=-1
		galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = galaxy.Vls, dist=meanDist)
	      
	    g_lst = []
	    for galaxy in galList:
	      if galaxy.id in virgo_g2:
		galaxy.inGroup = 0
		g_lst.append(galaxy)
		galaxy.setMeanDist(16, 0.1, GRP_vel = galaxy.Vls, dist=15.7)
	    
	    new_grp = IndividualGal(g_lst)
	    G_list += new_grp
	    meanDist    = new_grp[0][0].mDist
	    meanDistErr = new_grp[0][0].mDistErr
	    for galaxy in g_lst:
	      if galaxy.inGroup == 0: 
		#galaxy.inGroup=-1
		galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = galaxy.Vls, dist=meanDist)
	    
	    
	    addTo(G_list, galList, 39483, 39537)
	    addTo(G_list, galList, 39483, 39620)
	    addTo(G_list, galList, 39483, 39809)
	    #addTo(G_list, galList, 39483, 3394197)
	    
	    G_list = addGalGroup(G_list, galList)
	    
	    
	    manual_group(G_list, galList, Virgo_W1)
	    manual_group(G_list, galList, Virgo_W2)
	    manual_group(G_list, galList, Virgo_M)
	    
	
		
	    #print id_Virgo_w
	    VirW_grp = manual_group(G_list, galList, id_Virgo_w, verbose=False)
	    #print "Virgo W length: " ,len(VirW_grp)
	    

	    
	    
            ### VVVVVIIIIIRRRRGGGGOOOOO
            Virgo_grp = []
            for galaxy in galList: 
	      if galaxy.id in Virgo_ignore and galaxy.inGroup == -1:
		Virgo_grp.append(galaxy.id)
            if len(Virgo_grp)>0:
	      print "Working on generating the Virgo cluster ..."
	      manual_group(G_list, galList, Virgo_grp)
	    
	    
            ### Assumption: any galaxy without a known distance that is within 1 deg of 39659
            ### and has V_ls > 1200 km/s is in 39659 (Virgo W)
	    for galaxy in galList:
	      if (180./pi)*angle(galaxy.sgl, galaxy.sgb, 108.3771, -6.9394)<1.0 and  galaxy.Vls>1200 and galaxy.dcf2==0:
		addTo(G_list, galList, 39659, galaxy.id)
              # Virgo negative velocities
	      if (180./pi)*angle(galaxy.sgl, galaxy.sgb, sglV, sgbV)<6.8 and  galaxy.Vls<0:
		addTo(G_list, galList, 41220, galaxy.id)
              

	    G_list = addGalGroup(G_list, galList)
	    G_list += IndividualGal(galList)
	    removeFrom2(G_list, 39659, 40001)
	    trimed_virgoW = trimGroup_vel(G_list, 39659,  Vmin=1200)
	    
	  
	    manualBlend(G_list, 40001, 40240)
	    addTo(G_list, galList, 40001, 40851)
	    
	    
	    Grp40001 = []
	    for i in range(len(G_list)):
	      if G_list[i][0].nest == 40001:
		print "group 40001 exists ..."
		Grp40001.append(G_list[i])
		G_list.pop(i)
		break
	      
      
	    if len(Grp40001) > 0 and trimed_virgoW!=None: 
	       print "adding to 40001 group ..."
	       Grp40001 = addGalGroup(Grp40001, trimed_virgoW)
	       G_list.append(Grp40001[0])
	    
	    
	    brent_modifications(G_list, galList)
    
   
	    for id in Virgo_foreback:
	      removeFrom(G_list, id)

        ## Very last stage ... No group formation after this stage
        ### Hard coding groups
        if cluster!='manual':
	    manual_group(G_list, galList, id_MW)
	    manual_group(G_list, galList, id_M31) 
	    addTo(G_list, galList, 2557, 3097691)   # andromeda 
	    addTo(G_list, galList, 2557, 5060429) 
	    addTo(G_list, galList, 2557, 4608690) 
	    addTo(G_list, galList, 2557, 5056919) 
	    addTo(G_list, galList, 2557, 5057226) 
	    addTo(G_list, galList, 2557, 5057227) 
	    addTo(G_list, galList, 2557, 5057228) 
	    addTo(G_list, galList, 2557, 5057229) 
	    addTo(G_list, galList, 2557, 5057230)   
   
	
	
	
        group_moderator(G_list, use_dist=True)
        galList_moderator(galList)


        groupWrite(cluster+'.iter.'+str(int(iter))+'.v43.group', G_list, galList)
        #sys.exit(1)

	for galaxy in galList: 
	  if galaxy.inGroup > 0 : 
	     galaxy.inGroup = 0
        

  
