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
   L_tot = L1 + L2
   
   logK_tot = log10(L_tot)
   
   
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


   n1 = galNode1.subGalaxies
   n2 = galNode2.subGalaxies   
   Vls = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)

   newNode = GalxyNode(ID, gl, gb, sgl, sgb, Vls, 0, 0, d, 0)
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
        meanDistErr = sqrt(1/newNode.sumError)/meanDist
   else:
        meanDist = 0.
        meanDistErr = 0.
    
   
   
   if galNode1.mDist != 0 and galNode2.mDist != 0:
     if L1 > L2: 
       newNode.mDist = galNode1.mDist
       newNode.mDistErr = galNode1.mDistErr
     else:
       newNode.mDist = galNode2.mDist
       newNode.mDistErr = galNode2.mDistErr
   elif galNode1.mDist != 0:
     newNode.mDist = galNode1.mDist
     newNode.mDistErr = galNode1.mDistErr
   elif galNode2.mDist != 0:
     newNode.mDist = galNode2.mDist
     newNode.mDistErr = galNode2.mDistErr
   else:
     newNode.mDist = meanDist
     newNode.mDistErr = meanDistErr

   
   # Brighter galaxy is the left child
   if L1 >= L2:
      newNode.left = galNode1
      newNode.right = galNode2
      newNode.nest = galNode1.nest
   else:
      newNode.left = galNode2
      newNode.right = galNode1
      newNode.nest = galNode2.nest
   

   newNode.v_av = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)
   newNode.v2_av = (n1*galNode1.v2_av + n2*galNode2.v2_av) / (n1+n2)
   
   if (newNode.v2_av - newNode.v_av**2) > 0:
      newNode.sigma =  sqrt(newNode.v2_av - newNode.v_av**2) 
   
   
   newNode.R_theta = Theta_max(newNode)
   

   if newNode.sigma == 0:
     sig = 1
   else:
     sig = newNode.sigma
   

   mass = Mass(L_tot)
   newNode.M_v2 = mass
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
  def __init__(self, id, gl, gb, sgl, sgb, Vls, Ks, Ty, dcf2, ed):
    
    
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
    
    self.dcf2 = dcf2
    self.ed = ed
    self.dcf2Copy = dcf2
    self.edCopy = ed    
    
    self.mDist = dcf2
    self.mDistErr = ed
    
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
    
    self.logK = m_logK(Ks, Vls)
 
    self.M_v2 = Mass(10**self.logK)
    self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc


    
    if dcf2!=0 and ed!=0:
        errDist = dcf2 * ed
        self.sumDist = 1. * dcf2 / (errDist**2)
        self.sumError = 1. / (errDist**2)
    else:
        self.sumDist = 0.
        self.sumError = 0.
        
    
    
  def setMeanDist(self, Dist, errDist, GRP_vel = 0):
    
      self.mDist = Dist
      self.mDistErr = errDist
      if GRP_vel == 0 : 
	vel = self.Vls
      else: vel = GRP_vel
      self.logK = m_logK(self.Ks, vel)
      self.M_v2 = Mass(10**self.logK)
      self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
      
    


# **************************************
# Given the apparent magnitude of a galaxy, m
# it returns the absolute luminosity 
# Vls is the radial velocity of the galaxy
# So, the distance to a galaxy is estimated 
# based on its velocity
def m_logK(m, Vls):
    
    Mk_sun = 3.28   # K-band
    distance = Vls / H0  # Mpc
    
    if distance < 0:
      distance = 1 # Mpc
    
    M = m - 5*log10(distance) - 30 + 5
    logK = -0.4 * (M - Mk_sun)
  
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
def readgalList(table, skyPatch='none', R_min=6.0, R_max=100000000., sgl_min=0, sgl_max=0, sgb_min=0, sgb_max=0):

  id   = table['pgc']
  sgl  = table['sgl']
  sgb  = table['sgb']
  Vls  = table['Vls']
  Vhelio = table['Vhelio']
  Ks   = table['Ks'] 
  dcf2 = table['dcf2'] 
  ed   = table['ed'] 
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
  
  fornax = False
  virgo_w = False
  manual = False
  if skyPatch == 'fornax':
    fornax = True
  elif skyPatch == 'virgo_w':
    virgo_w = True
  elif skyPatch=='manual':
    manual = True
     

  i = 0
  
  for i in range(N_galaxies):
    
    Bool = False
    if fornax:
      if Vls[i] < 4000 and Ks[i] > 0 and sgl[i]<303 and sgl[i]>240 and sgb[i]>-57 and sgb[i]<-20:
        Bool = True
    
    elif virgo_w:
      ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)
      if Vls[i] < 4000 and Ks[i] > 0 and sgl[i]<117 and sgl[i]>105 and sgb[i]>-10 and sgb[i]<2 and ang12 > R_min and ang12 < R_max:
        Bool = True
    elif  manual:
      ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)
      if Vls[i] < 4000 and Ks[i] > 0 and sgl[i]<sgl_max and sgl[i]>sgl_min and sgb[i]>sgb_min and sgb[i]<sgb_max and ang12 > R_min and ang12 < R_max:
        Bool = True
    else:
      ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)
      if Vls[i] < 4000 and Ks[i] > 0 and ang12 > R_min and ang12 < R_max: 
	#if (180./pi)*angle(sgl[i], sgb[i], 328.7093, -6.0825) < 5: # 255.3245 |   5.6045  --- 334.8504 |  27.0145 -- 328.7093 |  -6.0825 
	  Bool = True
    if Bool:         
	   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[i], Vls[i], Ks[i], Ty[i], dcf2[i], ed[i])
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
  
  
  NoGroups = len(G_list)
  #print "Number of groups: ", NoGroups  

  myTable = Table()
    
    
  empty = []
  myTable.add_column(Column(data=empty,name='pgc', dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='flag', dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='ra', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='dec', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='gl', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='gb', format='%0.4f'))    
  myTable.add_column(Column(data=empty,name='sgl', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='sgb', format='%0.4f'))
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
  


  
  for i in range(0, NoGroups):  # for all groups
    if G_list[i][0]. Vls <= 3500:
       dist = G_list[i][0].mDist
       if dist == 0 : dist = G_list[i][0].Vls / H0
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
        sumDist = galaxy.sumDist
        sumError = galaxy.sumError
        
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
        tX_dyn = Mpc_km*0.5*pi*Rg/sigmaP_dyn/sqrt(2.5)
        R2t_dyn = (Rg*pi*0.5)/1.05/sqrt(1.5)
       
 
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,Ks,logK,Vls, Vhelio, dcf2, ed, \
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_lum, \
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, tX_dyn, tX_lum, subGalaxies, nest, \
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname])
	

  
  # writing individual galaxies
  for galaxy in galList:


      if galaxy.inGroup == 0. and galaxy.Vls<=3500:
	flag = 0
	galaxy.flag = 0
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
        sumDist = galaxy.sumDist
        sumError = galaxy.sumError
        
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
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname])

	

  
  #myTable.write(outfile, format='ascii',delimiter=',')
  
  myTable.write(outfile, format='ascii.fixed_width',delimiter='|', bookend=False)
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

def fornaxPlot(G_list, galList):
  
    

    
    id   = []

    RA0  = []
    DEC0 = []    
    dcf2 = []
    Vls  = []
    
    
    for i in range(len(galList)):
       id.append(galList[i].id)
       RA0.append(galList[i].ra)
       DEC0.append(galList[i].dec)
       dcf2.append(galList[i].dcf2)
       Vls.append(galList[i].Vls)
    
    N_galaxies = len(id)
   
    
    
    NoGroups = len(G_list)
    print "Number of groups: ", NoGroups
  

    
    
    fig = plt.figure(figsize=(7, 7), dpi=100)
    
    
    ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))  
  
  
    plt.minorticks_on()
    plt.tick_params(which='major', length=7, width=1.5)
    plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     
    
    plt.plot(RA0, DEC0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
    #plt.ylim(-45,-15)
    #plt.xlim(69,39)
    
    plt.ylim(-45,-10)
    plt.xlim(70,35)
    
    
    plt.xlabel("RA (deg)", fontsize=20)
    plt.ylabel("DEC (deg)", fontsize=20)
    
    for i in range(0, NoGroups):  # for all groups
  
      Key = True
      if Key:
	  random.seed(G_list[i][0].nest )
	  

	  Red, Green, Blue = random.random(), random.random(), random.random()
	  
	  

	  Dist = G_list[i][0].Vls/H0
	  if Dist<1:
	       Dist = 1.
	  r = 180.*atan(G_list[i][0].R_2t2/Dist)/pi
	  
	  
	  
	  d_theta = 0.01
	  u = np.arange(0,2*pi,d_theta)
	  X = u*0.
	  Y = u*0.

	  RA0 = G_list[i][0].ra
	  DEC0 = G_list[i][0].dec
	  
	  for q in range(0,len(u)):
	    x = r*np.cos(u[q]) + RA0
	    y = r*np.sin(u[q]) + DEC0
	    X[q] = x
	    Y[q] = y

	  
	  if r <= 10000:
	    if G_list[i][0].nest == 13418:   # NGC_1399
	      Blue = 0; Red = 0.2 ; Green = 0.8
	      
	    line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Green, Blue)) 
	    line.set_dashes([8, 3]) 
	    plt.text(RA0, DEC0+r, int(G_list[i][0].nest), fontsize=8, color=(Red, Green, Blue))
	    

	  X  = []
	  Y  = []
	  

	  

	  for j in range(1, len(G_list[i])):

	      RA0 = G_list[i][j].ra
	      DEC0 = G_list[i][j].dec
	      X.append(RA0)
	      Y.append(DEC0)
	     

	  plt.plot(X, Y, 'o', markersize = 3, color=(Red, Green, Blue), markeredgecolor = 'none')
	  if G_list[i][0].nest == 13418:   # NGC_1399
	    plt.plot(G_list[i][0].ra, G_list[i][0].dec, '*', markersize = 10, color= 'black')

################################################################# 
################################################################# 
#################################################################
#################################################################

def virgoPlot(G_list, galList, R_min=6.0, R_max=100000000., sgl_min=0, sgl_max=0, sgb_min=0, sgb_max=0):
  
  id   = []
  sgl  = []
  sgb  = []
  Vls  = []
  dcf2 = []
  
  for i in range(len(galList)):
     id.append(galList[i].id)
     sgl.append(galList[i].sgl)
     sgb.append(galList[i].sgb)
     dcf2.append(galList[i].dcf2)
     Vls.append(galList[i].Vls)
    
  N_galaxies = len(id)
  NoGroups = len(G_list)
  print "Number of groups: ", NoGroups

  
  fig = plt.figure(figsize=(7, 7), dpi=100)
  ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
  ax.xaxis.set_major_locator(MultipleLocator(10))
  ax.yaxis.set_major_locator(MultipleLocator(10))
  ax.xaxis.set_minor_locator(MultipleLocator(1))
  ax.yaxis.set_minor_locator(MultipleLocator(1))  
  
  
  
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)  
  
  if sgl_max == 0: 
      xmin = 5*floor((sglV+35)/5.)+3
      xmax = 5*ceil((sglV-35)/5.)-3
      #ymax = 5*floor((sgbV+35)/5.)+3
      ymin = 5*ceil((sgbV-35)/5.)-3
      ymax = ymin - (xmax-xmin)
  else:
      xmin = sgl_max
      xmax = sgl_min
      ymin = sgb_min
      ymax = sgb_max
  
  
  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)
  
  
  plt.xlabel("SGL (deg)", fontsize=20)
  plt.ylabel("SGB (deg)", fontsize=20)
  plt.yticks(fontsize=14)
  plt.xticks(fontsize=14)
  
  # plotting big slod circle
  theta = np.arange(0,1,0.001)
  Circlx = 30.*np.cos(theta*2*pi) + sglV
  Circly = 30.*np.sin(theta*2*pi) + sgbV
  plt.plot(Circlx, Circly, 'k-')
  
  # plotting small balck dotted circle
  theta = np.arange(0,1,0.010)
  Circlx = 6.8*np.cos(theta*2*pi) + sglV
  Circly = 6.8*np.sin(theta*2*pi) + sgbV
  line00, = plt.plot(Circlx, Circly, 'k.', markersize = 2)   
  line00.set_dashes([2, 2]) 
  
  # plotting small red dotted circle
  Circlx = R_min*np.cos(theta*2*pi) + sglV
  Circly = R_min*np.sin(theta*2*pi) + sgbV
  line0, = plt.plot(Circlx, Circly, 'r.', markersize = 2)   
  line0.set_dashes([2, 3]) 
  
  # plotting small red dotted circle
  Circlx = R_max*np.cos(theta*2*pi) + sglV
  Circly = R_max*np.sin(theta*2*pi) + sgbV
  line0, = plt.plot(Circlx, Circly, 'g.', markersize = 2)   
  line0.set_dashes([2, 3])   
  

  NoGroups = len(G_list)
  print "VirgoPlot >> Number of groups: ", NoGroups
  plt.plot(sgl, sgb, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  
  col = 0
  colors= ['#0000FF','#CC00CC','#999900','#006600','#FF0000','#FF8000', '#6600CC'] 
  for i in range(0, NoGroups):  # for all groups
  
    groupFromVirgo = angle(G_list[i][0].sgl, G_list[i][0].sgb, sglV, sgbV)*180./pi
    
    if groupFromVirgo>=6.8 and groupFromVirgo<=1000:
    #if G_list[i][0].id==200041220:
	random.seed(G_list[i][0].nest )
	Red, Blue, Green = random.random(), random.random(), random.random()
	
	

	Dist = G_list[i][0].Vls/H0
	if Dist <= 1 : Dist = 1
	r = 180.*atan(G_list[i][0].R_2t2/Dist)/pi
	
	d_theta = 0.001
	theta = np.arange(0,2*pi,d_theta)
	Circlx = r*np.cos(theta) + G_list[i][0].sgl
	Circly = r*np.sin(theta) + G_list[i][0].sgb
	
	if r <= 7:
	  line, = plt.plot(Circlx, Circly, '-', markersize = 2, color=(Red, Blue, Green)) 
	  line.set_dashes([8, 3]) 
	

	sRA  = []
	sDEC = []

	for j in range(1, len(G_list[i])):
 
	    sRA.append(G_list[i][j].sgl)
	    sDEC.append(G_list[i][j].sgb)
        plt.plot(sRA, sDEC, 'o', markersize = 4, color=(Red, Blue, Green), markeredgecolor = 'none')
        
        
#################################################################
def mergGroup(G_list, restrict=False):
    

    
    NoGroups = len(G_list)
 
    i = 0
    while i < NoGroups:
      j = i+1
      while j < NoGroups-1:
	


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
        
        
	if ang12 <= 1.1*(r1+r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(sig1,sig2) and (min(r1,r2))**3 < 0.2*(max(r1,r2))**3:
	       Bol = True        
        
	if ang12 <= 0.6*(r1+r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True

	      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True
	
	      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt < min(e1, e2):
	      Bol = True
	      
	      
	# one group completely projected on another one (e.g. smal in big)
	if ang12 <= 1.1*(max(r1,r2)-min(r1,r2)) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(4*Sigquad,4.0*max(sig1,sig2)):
	      Bol = True	      
#######################################################################

        if (idd1==39241 and  idd2==39600) or (idd1==39600 and  idd2==39241):
	  Bol = False
	if (idd1==13324 and  idd2==13108) or (idd1==13108 and  idd2==13324):
	  Bol = True    	  
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

def addGalGroup(G_list, galList, ignore_list=[]):
  

    singles = []
    for galaxy in galList:
       if galaxy.inGroup == 0:
	  singles.append([galaxy, -1, -10000])  # [galaxy, new.group.id, angular.separation]
    
    N = len(G_list)
    BB = []
    LL = []
    for i in range(N):
      LL.append(G_list[i][0].gl)
      BB.append(G_list[i][0].gb)
    
    for entity in singles:
            
            gl = entity[0].gl
            gb = entity[0].gb
    
	    q1 = np.zeros((N,), dtype=np.int)
            q2 = np.zeros((N,), dtype=np.int)
            q3 = np.zeros((N,), dtype=np.int)
            q4 = np.zeros((N,), dtype=np.int)
	    q1[np.where(LL<=gl+10)] = 1
	    q2[np.where(LL>=gl-10)] = 1
	    q3[np.where(BB<=gb+10)] = 1
	    q4[np.where(BB>=gb-10)] = 1
	    qq = q1+q2+q3+q4
	    group_indices = np.where(qq==4)
	       
	    for i in group_indices[0]:

	      ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0].gb, gl, gb)
              d = G_list[i][0].Vls / H0
              if d < 1: d = 1 
              r = (180.*atan(G_list[i][0].R_2t2/(d))/pi)	      
	      
	      # if the galaxy is close enough to the center of the group
	      # and it's not already in any other group, then it would be added to the current group
	      if ang12<=10 and ang12 <= 2*r:
	        
	        sig_p = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
	        
	        #force = G_list[i][0].M_v2 / (ang12**2)

	        if ang12 <= 1.05*(180./pi)*atan(1.3*G_list[i][0].R_2t2/d) and abs(entity[0].Vls - G_list[i][0].Vls) < 3.0*sig_p:
	            d1 = G_list[i][0].mDist
	            e1 = d1 * G_list[i][0].mDistErr
	            d2 = entity[0].dcf2
	            e2 = d2 * entity[0].ed
	            delt = abs(d1-d2)
	            join = False
                    if ang12 > 1.05*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d) and d1!=0 and d2!=0 and delt<max(e1,e2):
		      join = True
		    if d1==0 or d2==0:
		      join = True
		    if ang12 <= 1.05*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d):
		      join = True
		    if join==True and  sig_p > entity[2]:
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
	        
	    
  # This is the output	    
  #G_list = NEWgpLIST
  return NEWgpLIST


#################################################################
def IndividualGal(galList_org, coeff=1.0, ignore_list=[]):
  
  
  galList = []
  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000 and (galaxy.id not in ignore_list): 
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

      thet = (180.*atan(1.0*grp.R_2t2/(grp.Vls/H0))/pi)
      if grp.Vls<75:
	thet = (180.*atan(1.0*grp.R_2t2/(1.))/pi)
      
      coeff = 1.05
      if ang12 <= coeff*thet:
	sig_p = (grp.M_v2 / (2.0E6))**(1./3)

	if (abs(grp.Vls-individuals[q].Vls) <= 2.0*sig_p):

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
	    
      if grp.id == 40001: Bol = False
      q = (q+1)%N
    if Bol: #len(galaxies) > 1:
      NEWgpLIST.append(NodeLeaf(grp))
    p+=1
    
    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1
    
  return NEWgpLIST

#################################################################

def LoosgalListJoint(galList, grID = 500000000):
   
   if len(galList)==0: return None
   
   Lum_mass = []
   for i in range(0, len(galList)):
        Lum_mass.append(galList[i].Ks)
   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass)
   
   
   
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
   meanDist = 0.
   meanDistErr = 0.
   if sumDist !=0 and sumError!=0:
      meanDist = sumDist/sumError
      meanDistErr = sqrt(1/sumError)/meanDist
    
   G_dist = meanDist
   G_dist_err = meanDistErr
   for gal in galList:
      gal.logK = m_logK(gal.Ks, V_mean)
      gal.mDist = G_dist
      gal.mDistErr = G_dist_err
     
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
def group_moderator(G_list):

  for group in G_list:
	  groupHeader = group[0]
	  groupHeader.flag = 2
	  groupHeader.dcf2 = 0
	  groupHeader.ed = 0
	  groupHeader.dcf2Copy = 0
          groupHeader.edCopy = 0
          
          sumDist = 0
          sumError = 0
          for galaxy in group[1:]:
	    sumDist += galaxy.sumDist
	    sumError += galaxy.sumError
	  
	  meanDist = 0.
	  meanDistErr = 0.
	  if sumDist !=0 and sumError!=0:
	      meanDist = sumDist/sumError
	      meanDistErr = sqrt(1/sumError)/meanDist
          

	  L_tot = 0
	  ID = groupHeader.id 
	  groupHeader.id = ID - ID % 100000000 + groupHeader.nest
	  
	  V_hel_sum = 0.
	  
	  for galaxy in group[1:]:
	    galaxy.flag = 1
	    # and also modify the absolute luminosities
	    galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls)
	    L_tot += 10**galaxy.logK
	    V_hel_sum += galaxy.Vhelio
	  
	  groupHeader.Vhelio = V_hel_sum / (len(group)-1)
          groupHeader.logK = log10(L_tot)
          Dist_v =  groupHeader.Vls / H0
          if Dist_v<1. : Dist_v=1
          Mk_sun = 3.28   # K-band
          M = Mk_sun - (groupHeader.logK / 0.4)
          groupHeader.Ks = M + 5*log10(Dist_v) + 30 - 5
  
  return G_list
      
#################################################################

if __name__ == '__main__':
  

  R_max=50. # R_max=50.
  cluster = ''

  
  if len(sys.argv) < 2:
    print "\nEnter the cluster name as the input ..." 
    print "\nexample: \n > python " + str(sys.argv[0])+ " fornax"
    print "\nPossible options: virgo, fornax, north, south, manual" 
    print "Use north/south for whole sky.\n"
    sys.exit(1)
  
  cluster = str(sys.argv[1])
  
  if cluster == 'virgo':
    inFile  = 'AllSky.north.v13.csv'
  elif cluster == 'virgo_w':
    inFile  = 'AllSky.north.v13.csv'
  elif cluster == 'fornax':
    inFile  = 'AllSky.south.v13.csv'
  elif cluster == 'north':
    inFile  = 'AllSky.north.v13.csv'
  elif cluster == 'south':
    inFile  = 'AllSky.south.v13.csv'   
  else:
    print "\nInvalid Cluster Name ..." 
    print "Possible options: virgo, fornax, north and south" 
    print "Use north/south for whole sky.\n"
    
    inFile = 'AllSky.north.v13.csv'
    #sys.exit(1)
  
  table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
  

  
  if   cluster == 'virgo' :
     galList = readgalList(table, R_max=R_max)   
  elif cluster == 'virgo_w' :
     galList = readgalList(table, R_max=R_max, skyPatch='vigo_w')  
  elif cluster == 'fornax': 
     galList =  readgalList(table, skyPatch='fornax')
  elif cluster == 'north' :
     galList = readgalList(table)
  elif cluster == 'south' :
     galList = readgalList(table)
  else:
     alpha = 102.8806
     delta = -2.3479
     step = 10
     
     galList = readgalList(table, skyPatch='manual', sgl_min=alpha-step, sgl_max=alpha+step, sgb_min=delta-step, sgb_max=delta+step)    
  
  

  
  
  for iter in range(4):
  
        print "working on grouping ... iter:", iter
	##  Finding groups based on the criteria
        G_list = []
        ignore_list=[41618, 40001, 42447, 40240, 40851]
        for qp in range(0,3):
	  G_list += IndividualGal(galList, ignore_list = ignore_list)
	  
	# fixed
        G_list = addGalGroup(G_list, galList, ignore_list = ignore_list)
        if True:
          for pq in range(0,10): 
	    G_list = IndividualGal_modifier(G_list, galList)
	  for pq in range(0,10): 
	    G_list = addGalGroup(G_list, galList, ignore_list = ignore_list)
          for qp in range(0,10):
            mergGroup(G_list)
          mergGroup(G_list, restrict=True)
	  for pq in range(0,10): 
	    G_list = addGalGroup(G_list, galList, ignore_list = ignore_list)        
          for qp in range(0,2):
            mergGroup(G_list)
          mergGroup(G_list, restrict=True)      


####################################################################

	# Error: galList is modified ?
	for galaxy in galList:
	  if galaxy.inGroup == 0:
	    galaxy.mDist = 0
	    galaxy.setMeanDist(galaxy.dcf2, galaxy.ed, 0)
	    
	  
	
        group_moderator(G_list)
          
 
 


        if cluster == 'virgo':
	  groupWrite('virgo.iter.'+str(int(iter))+'.v24.group', G_list, galList)
        elif cluster == 'fornax': 
	  groupWrite('fornax.iter.'+str(int(iter))+'.v24.group', G_list, galList)
        elif cluster == 'north':
          groupWrite('north.iter.'+str(int(iter))+'.v24.group', G_list, galList)
        elif cluster == 'south':
	  groupWrite('south.iter.'+str(int(iter))+'.v24.group', G_list, galList)
        elif cluster == 'manual':
	  groupWrite('manual.iter.'+str(int(iter))+'.v24.group', G_list, galList)

	
	for galaxy in galList: galaxy.inGroup = 0
        

  
