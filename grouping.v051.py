#!/home/ehsan/Ureka/Ureka/variants/common/bin/python


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
from maxheap import *
import random
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py

from astropy import coordinates as coord
from astropy import units as unit
# **************************************
#    Global Variables
# **************************************
H0 = 75.           # hubble constant
sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude
G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.08567758E19
#t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]

# This is the val;use we are gonna adopt for the paper
t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr
# **************************************
def Mass(L_k):

  L = L_k / 1.E10
  
  if L < 1.:
    MtoL = 43.0
  elif L > 1000.:
    MtoL = 121.19
  else:
    MtoL = 43*(10**(0.15*log10(L)))
  
  Mass_out = L_k * MtoL
  
  return Mass_out

# **************************************

def Rho_cr(L_k, ML):
  
  L = L_k / 1.E10
  
  if L < 1.:
    MtoL = 43.0
  elif L > 1000.:
    MtoL = 121.19
  else:
    MtoL = 43*(10**(0.15*log10(L)))
  
  if ML <= 43: MtoL = ML
  
  Rho_cr = (Mpc_km*1000)**3/(4*pi*G*t0_age**2)/M_sun/MtoL
  Rho_cr_universe = 1.5 * Rho_cr
  Rho_200 = 200 * Rho_cr_universe
  
  Tx = (pi/2.) * Mpc_km / sqrt(2.5) / 350.
  
  Rho_X = (Mpc_km*1000)**3/(4*pi*G*Tx**2)/M_sun/MtoL
  
  return Rho_X
  
  

# **************************************
# returns radian
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


# **************************************
def distance(l1, b1, v1, d1, l2, b2, v2, d2):
   
   Vl = 300   # km/s  ... Limiting velocity
   
   theta = angle(l1, b1, l2, b2)
   R = 0.
   
   #if d1!=0 and d2!=0:
      #R = sqrt(d1**2+d2**2-2.*d1*d2*cos(theta))
   #elif d1!=0:
      #v1 = d1 * H0
   #elif d2!=0:
      #v2 = d2 * H0
     
   
   
   if R == 0:
	  V12 = abs (v1-v2)
	  D = 0.5 * (v1 + v2) / H0
	  
	  R = 0  # Separation distance
	  if V12 <= Vl:
	    R = (8./pi)*D*sin(theta/2)
	    #print "check point"
	    #print v1, v2, D, theta*90/pi
	  else:
	    R1 = (2*D*sin(theta/2))**2
	    R2 = (1 + (4/pi -1)*(Vl**2/V12**2))**2
	    R3 = (V12**2-Vl**2)/H0**2
	    R = sqrt(R1*R2 + R3)
	  
	  if R==0:
	    R = 1.E-6
	    
   
   return R
 
# **************************************

def sign(x):
  if x<0 : return -1
  if x>=0 : return 1
# **************************************
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
   
   l1 = galNode1.l
   l2 = galNode2.l
   b1 = galNode1.b
   b2 = galNode2.b
   
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
     
     
   l, b, d   = barycenter(L1,  l1,  b1, d1, L2,  l2,  b2, d2)
   gl, gb, d = barycenter(L1, gl1, gb1, d1, L2, gl2, gb2, d2)
   ra, dec, d = barycenter(L1, ra1, dec1, d1, L2, ra2, dec2, d2)
   
   if d_final==0:
      d = 0

   

   
   n1 = galNode1.subGalaxies
   n2 = galNode2.subGalaxies   
   Vls = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)

   newNode = GalxyNode(ID, gl, gb, l, b, Vls, 0, 0, d, 0)
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
   
   
   
   dist = distance(l1, b1, v1, d1, l2, b2, v2, d2)
   dist1 = distance(l, b, Vls, d, l1, b1, v1, d1)
   dist2 = distance(l, b, Vls, d, l2, b2, v2, d2)   
   
   


   
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
  sigma = 0.  # velocity dispersion
  dcf2 = 0.
  ed = 0.
  dcf2Copy = 0.
  edCopy = 0.
  mDist = 0.
  mDistErr = 0.
  Ks = -100000.0
  dend_x = 0.
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
  
  keyData = None  # pointing at cross calculated data (L/r**2)
  left = None
  right = None
  
  def __init__(self, id, gl, gb, l, b, Vls, Ks, Ty, dcf2, ed):
    
    
    #dcf2 = 0
    self.id = id
    self.gl = gl
    self.gb = gb    
    self.l = l
    self.b = b
    self.Vls = Vls
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
    self.keyData = None
    
    self.subGalaxies = 1
    self.level = 0
    self.Rgroup = 0.  
    #self.Rho_gr = 0.
    self.R_theta = 0.
    self.nest = id
    self.sigma = 0. 
    self.v_av = Vls
    self.v2_av = Vls**2
    
    
    self.dend_x = 0.
    self.inGroup = 0.

    
    self.logK = m_logK(Ks, Vls, 0, dcf2)
    
    #self.r_g = 0.
    #self.M_v = 0.
    #self.R_2t = 0.
    
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
    
    if self.dcf2 == 0 and self.ed == 0:
      self.mDist = Dist
      self.mDistErr = errDist
      if GRP_vel == 0 : 
	vel = self.Vls
      else: vel = GRP_vel
      self.logK = m_logK(self.Ks, vel, 0, 0)
      self.M_v2 = Mass(10**self.logK)
      self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
      
    
    
    
  
  
  
  def toString(self):
    
    print "\n"
    print "PGC: ", self.id
    print "ra: ", self.ra
    print "dec: ", self.dec   
    print "gl: ", self.gl
    print "gb: ", self.gb    
    print "sgl: ", self.l
    print "sgb: ", self.b
    print "Vls: ", self.Vls
    print "logK: ",  self.logK
    print "Ty: ", self.Ty
    print "R_theta: ",  self.R_theta
    print "Sigma: ", "%7.2f" % self.sigma
    print "# of subGal: ",  self.subGalaxies
    print "level: ",  self.level
    print "nest: ",  self.nest
    
    
    if self.left != None:
      print "left PGC ID: ", self.left.id
    
    if self.right != None:
      print "right PGC ID: ", self.right.id  
    
    print "\n"

  
  
# **************************************
def minDist(galList):
  
  length = len(galList)
  
  minDist = 10000000
  minIDi = 0
  minIDj = 0
  
  for i in range(0, length-1):
    for j in range(i+1, length):
       
       dist12 = distance(galList[i].l, galList[i].b, galList[i].Vls, galList[i].dcf2, galList[j].l, galList[j].b, galList[j].Vls, galList[j].dcf2)
       if dist12 < minDist:
          minIDi = i
          minIDj = j
          minDist = dist12
   
  return minIDi, minIDj, minDist, galList[minIDi].id, galList[minIDj].id
    
  
# **************************************
def maxKey(galList):
  
  length = len(galList)
  
  maxKey = 0
  maxIDi = 0
  maxIDj = 0
  
  for i in range(0, length-1):
    for j in range(i+1, length):
       
       key12 = pairKey(galList[i], galList[j])
       if key12 > maxKey:
          maxIDi = i
          maxIDj = j
          maxKey = key12
   
  return maxIDi, maxIDj, maxKey, galList[maxIDi].id, galList[maxIDj].id



# **************************************

def pairKey(galNode1, galNode2):
  
  dist12 = distance(galNode1.l, galNode1.b, galNode1.Vls, galNode1.dcf2, galNode2.l, galNode2.b, galNode2.Vls, galNode2.dcf2)
  L1 = 10**galNode1.logK
  L2 = 10**galNode2.logK
  L_k = L1+L2 # max([L1,L2])
  key = L_k/(dist12**2)

  return key
# **************************************
def m_logK(m, Vls, meanDist, distance):
    
    Mk_sun = 3.28   # K-band
    

    if distance==0 and meanDist==0:
       distance = Vls / H0  # Mpc
    elif meanDist!=0 :
       distance = meanDist  # Mpc
    
    
    M = m - 5*log10(distance) - 30 + 5
    logK = -0.4 * (M - Mk_sun)
  
    return logK
# **************************************
def groupID(root, ML):
  list = []
  groupIDCore(root, list, ML)
  
  for i in range(0, len(list)):
      for j in range(1, len(list[i])):
         list[i][j].inGroup = 1
  return list
  
  
def groupIDCore(root, list, ML):
  
  if root.left == None: return list
  

  Rg = root.R_2t2 / 1.22
  
  Dist = root.mDist
  if Dist == 0. :   Dist = root.Vls/H0
  #Rg_obs = root.R_theta*Dist
  angular_radii = Theta_max(root)
  Rg_obs = angular_radii*Dist
  

  if Rg_obs <= root.R_2t2 and angular_radii*180./pi < 10.:

    list.append(NodeLeaf(root))

  elif  root.left != None: 
    groupIDCore(root.left, list, ML)
    groupIDCore(root.right, list, ML)

# **************************************



def NodeLeaf(root):
  list = [root]
  NodeLeafCore(root, list)
  return list

def NodeLeafCore(root, list):  
  if root.left == None:
    #print "  ", root.id, root.l, root.b
    list.append(root)
  else:
    NodeLeafCore(root.left, list)
    NodeLeafCore(root.right, list)
    
    
    
# **************************************

def treeAppend(myTable, root):

    
    pgc = root.id 
    gl = root.gl  
    gb = root.gb  
    sgl = root.l  
    sgb = root.b  
    Vls = root.Vls  
    
    Ks =  root.Ks
    logK = root.logK  
    Ty = root.Ty  
    dcf2 = root.dcf2  
    ed = root.ed  
    
    dcf2Copy = root.dcf2Copy
    edCopy = root.edCopy
    
    
    ra  =  root.ra
    dec = root.dec
    coordinate_src = root.coordinate_src
    Ty_src  = root.Ty_src
    B_mag   = root.B_mag
    Ks_src  = root.Ks_src
    Vls_src = root.Vls_src
    objname = root.objname
    

    mDist = root.mDist
    mDistErr = root.mDistErr
    sumDist = root.sumDist
    sumError = root.sumError
    
    subGalaxies = root.subGalaxies  
    level = root.level  
    R_theta = root.R_theta  
    sigma = root.sigma  
    v_av = root.v_av  
    v2_av = root.v2_av  
    nest = root.nest  
    
    M_v2 = root.M_v2
    R_2t2 = root.R_2t2

    
    if root.left != None:
        leftID = root.left.id  
        rightID = root.right.id  
    else:
        leftID = 0  
        rightID = 0  
    

    myTable.add_row([pgc, ra, dec, gl, gb, sgl, sgb, Vls, Ks, logK, B_mag, Ty, dcf2, ed, \
      dcf2Copy, edCopy, mDist, mDistErr, sumDist, sumError, subGalaxies, level, \
	coordinate_src, Ty_src, Ks_src, Vls_src, objname, \
	  R_theta, sigma, v_av, v2_av, nest, M_v2, R_2t2, leftID, rightID])



# **************************************
def addSubTree(myTable, root):
  
  
  treeAppend(myTable, root)
  if root.left != None:
    addSubTree(myTable, root.left)
    addSubTree(myTable, root.right)
  
  
# **************************************

def writeTree (filename, root):
    
    
    myTable = Table()
    
    
    empty = []
    myTable.add_column(Column(data=empty,name='pgc'))
    myTable.add_column(Column(data=empty,name='ra'))
    myTable.add_column(Column(data=empty,name='dec'))
    myTable.add_column(Column(data=empty,name='gl'))
    myTable.add_column(Column(data=empty,name='gb'))    
    myTable.add_column(Column(data=empty,name='sgl'))
    myTable.add_column(Column(data=empty,name='sgb'))
    myTable.add_column(Column(data=empty,name='Vls'))
    myTable.add_column(Column(data=empty,name='Ks'))
    myTable.add_column(Column(data=empty,name='logK'))
    myTable.add_column(Column(data=empty,name='B_mag'))
    myTable.add_column(Column(data=empty,name='Ty'))
    myTable.add_column(Column(data=empty,name='dcf2'))
    myTable.add_column(Column(data=empty,name='ed'))
    
    
    myTable.add_column(Column(data=empty,name='dcf2Copy'))
    myTable.add_column(Column(data=empty,name='edCopy'))    
    myTable.add_column(Column(data=empty,name='mDist'))
    myTable.add_column(Column(data=empty,name='mDistErr'))
    myTable.add_column(Column(data=empty,name='sumDist'))
    myTable.add_column(Column(data=empty,name='sumError'))  
    myTable.add_column(Column(data=empty,name='No_Galaxies'))
    myTable.add_column(Column(data=empty,name='level'))
    
    myTable.add_column(Column(data=empty,name='coordinate_src', dtype='S15'))
    myTable.add_column(Column(data=empty,name='Ty_src', dtype='S10'))
    myTable.add_column(Column(data=empty,name='Ks_src', dtype='S15'))
    myTable.add_column(Column(data=empty,name='Vls_src', dtype='S15'))
    myTable.add_column(Column(data=empty,name='objname', dtype='S35'))
    

    myTable.add_column(Column(data=empty,name='R_theta'))
    myTable.add_column(Column(data=empty,name='sigma'))
    myTable.add_column(Column(data=empty,name='v_av'))
    myTable.add_column(Column(data=empty,name='v2_av'))
    myTable.add_column(Column(data=empty,name='nest'))
    myTable.add_column(Column(data=empty,name='M_v2'))
    myTable.add_column(Column(data=empty,name='R_2t2'))
    myTable.add_column(Column(data=empty,name='leftID'))
    myTable.add_column(Column(data=empty,name='rightID'))

    addSubTree(myTable, root)
 
    #myTable.write(filename+'.fits', format='fits', overwrite=True)
    myTable.write(filename, format='ascii',delimiter=',')
    


# **************************************
# input:  a row of the table which contain all properties of a node
# output: The genrated node based on the input properties 
def makeNode(table, index):
  
    id = table['pgc'][index]
    gl = table['gl'][index]
    gb = table['gb'][index]    
    sgl = table['sgl'][index]
    sgb = table['sgb'][index]
    Vls = table['Vls'][index]
    
    mDist = table['mDist'][index]
    mDistErr = table['mDistErr'][index]
    sumDist = table['sumDist'][index]
    sumError = table['sumError'][index]

    Ks  = table['Ks'][index]
    logK = table['logK'][index]
    Ty = table['Ty'][index]
    dcf2 = table['dcf2'][index]
    ed = table['ed'][index]
    dcf2Copy = table['dcf2Copy'][index]
    edCopy = table['edCopy'][index]
    
    ra  =  table['ra'][index]
    dec = table['dec'][index]
    coordinate_src = table['coordinate_src'][index]
    Ty_src  = table['Ty_src'][index]
    B_mag   = table['B_mag'][index]
    Ks_src  = table['Ks_src'][index]
    Vls_src = table['Vls_src'][index]
    objname = table['objname'][index]
    
    
    
    
    newNode = GalxyNode(id, gl, gb, sgl, sgb, Vls, Ks, Ty, dcf2, ed)
    newNode.mDist = mDist
    newNode.mDistErr = mDistErr
    newNode.sumDist = sumDist
    newNode.sumError = sumError
    
    newNode.dcf2Copy = dcf2Copy
    newNode.edCopy = edCopy
    
    newNode.ra = ra
    newNode.dec = dec
    newNode.coordinate_src = coordinate_src
    newNode.Ty_src = Ty_src
    newNode.B_mag =  B_mag
    newNode.Ks_src = Ks_src
    newNode.Vls_src = Vls_src
    newNode.objname = objname
    
    
    # middle node
    if Ks == -100000.0 : newNode.logK = logK
    
    
    newNode.subGalaxies = table['No_Galaxies'][index]
    newNode.level = table['level'][index]
    newNode.R_theta = table['R_theta'][index]
    newNode.nest = table['nest'][index]
    newNode.sigma = table['sigma'][index]
    newNode.v_av = table['v_av'][index]
    newNode.v2_av = table['v2_av'][index]
    
    newNode.M_v2 = table['M_v2'][index]
    newNode.R_2t2 = table['R_2t2'][index]
    
    
    
    
    
    
    
    return newNode
# **************************************
# Returns the root of the main tree
def makeTree(table, id):
  
  pgc   = table['pgc']
  leftID = table['leftID']
  rightID = table['rightID']
  
  index = 0 
  while pgc[index] != id: index+=1
  
  newNode = makeNode(table, index)
  
  left = leftID[index]
  right = rightID[index]
  
  #table.pop(index)
  
  if left != 0:
    newNode.left = makeTree (table, left)
    newNode.right = makeTree (table, right)
      
  
  return newNode

# **************************************


def readTree (filename):
  
    #skip_header=3
    table = np.genfromtxt( filename , delimiter=',', filling_values=0, names=True, dtype=None )
    
    
    return makeTree(table, table[0][0])

# **************************************
################################################################
################################################################
################################################################
################################################################
def Fornax_readgalList(table, R_min, R_max, n_gal, useDCF2=True):
  
 
  

  #print table.dtype
  id   = table['pgc']
  sgl  = table['sgl']
  sgb  = table['sgb']
  Vls  = table['Vls']
  Ks   = table['Ks'] 
  dcf2 = table['dcf2'] 
  ed   = table['ed'] 
  Ty   = table['Ty'] 
  
  gl = table['gl'] 
  gb = table['gb'] 
  
  
  ra  =  table['ra']
  dec = table['dec']
  coordinate_src = table['coordinate_src']
  Ty_src  = table['Ty_src']
  B_mag   = table['B_mag']
  Ks_src  = table['Ks_src']
  Vls_src = table['Vls_src']
  objname = table['objname']
  
  N_galaxies = len(id)

  
  
  print "No. of galaxies: ", N_galaxies;
  print "Making the main data structure .... "
  galList = []
  
  
  #RA = []
  #DEC = []
  #for i in range(0, N_galaxies):
      #point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
      #RA.append(point.fk5.ra.degree)
      #DEC.append(point.fk5.dec.degree)
    

  limit_counter = 0
  i = 0
  while(i < N_galaxies and limit_counter<n_gal ):
    
    if Vls[i] > 0 and Vls[i] < 4000 and Ks[i] > 0 and sgl[i]<303 and sgl[i]>240 and sgb[i]>-57 and sgb[i]<-20:  # -45
    #if Vls[i] > 0 and Vls[i] < 4000 and Ks[i] > 0 and sgl[i]<280 and sgl[i]>250 and sgb[i]>-40 and sgb[i]<-20:  # -45


       
       if useDCF2:
           node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[i], Vls[i], Ks[i], Ty[i], dcf2[i], ed[i])
	   node.ra = ra[i]
	   node.dec = dec[i]
	   node.coordinate_src = coordinate_src[i]
	   node.Ty_src = Ty_src[i]
	   node.B_mag =  B_mag[i]
	   node.Ks_src = Ks_src[i]
	   node.Vls_src = Vls_src[i]
	   node.objname = objname[i]
           galList.append(node)
       else:
	   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[i], Vls[i], Ks[i], Ty[i], 0, 0)
	   node.ra = ra[i]
	   node.dec = dec[i]
	   node.coordinate_src = coordinate_src[i]
	   node.Ty_src = Ty_src[i]
	   node.B_mag =  B_mag[i]
	   node.Ks_src = Ks_src[i]
	   node.Vls_src = Vls_src[i]
	   node.objname = objname[i]	   
	   node.dcf2Copy = dcf2[i]
	   node.edCopy = ed[i]
	   
	   galList.append(node)
       
       
       limit_counter+=1
    i+=1
  
  
  return galList 




################################################################
################################################################
def readgalList(table, R_min, R_max, n_gal, useDCF2=True):

  # names = true  ... takes the columns name from the header
  # dtype = true  ... recognizes the type of each column
  #inFile = 'brent.sorted.10.30.csv'
  #inFile = 'brent.notsorted.all.csv'
  

  #print table.dtype
  id   = table['pgc']
  sgl  = table['sgl']
  sgb  = table['sgb']
  Vls  = table['Vls']
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

  
  
  print "No. of galaxies: ", N_galaxies;
  print "Making the main data structure .... "
  galList = []
  
  
  


  limit_counter = 0
  i = 0
  while(i < N_galaxies and limit_counter<n_gal):
    ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)
    if Vls[i] > 0 and Vls[i] < 4000 and Ks[i] > 0 and ang12 > R_min and ang12 < R_max: # and sgl[i]<113 and sgl[i]>79  and sgb[i]<8 and sgb[i]>-35  :
       #Ks[i] = Median_k
       
       if useDCF2:
	   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[i], Vls[i], Ks[i], Ty[i], dcf2[i], ed[i])
	   node.ra = ra[i]
	   node.dec = dec[i]
	   node.coordinate_src = coordinate_src[i]
	   node.Ty_src = Ty_src[i]
	   node.B_mag =  B_mag[i]
	   node.Ks_src = Ks_src[i]
	   node.Vls_src = Vls_src[i]
	   node.objname = objname[i]
           galList.append(node)
           
       else:
	   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[i], Vls[i], Ks[i], Ty[i], 0, 0)
	   node.ra = ra[i]
	   node.dec = dec[i]
	   node.coordinate_src = coordinate_src[i]
	   node.Ty_src = Ty_src[i]
	   node.B_mag =  B_mag[i]
	   node.Ks_src = Ks_src[i]
	   node.Vls_src = Vls_src[i]
	   node.objname = objname[i]	   
	   node.dcf2Copy = dcf2[i]
	   node.edCopy = ed[i]
	   galList.append(node)

      
       
       limit_counter+=1
    i+=1
  
  
  return galList
  
################################################################
################################################################  
def galListJoint(galList, angle_limit = 10, grID = 1000000000):
  
  #angle_limit = 20   # degree
  list_dis_problem=[]
  
  t0 = time()
  numData = len(galList)
  galList[0].keyData = maxHeap()
  for i in range(1, numData):
       galList[i].keyData = maxHeap()
       j = i-1
       while j>=0:
	 
	    if angle(galList[i].l,galList[i].b,galList[j].l,galList[j].b)*180./pi < angle_limit :
	        key12 = pairKey(galList[i], galList[j])
                galList[i].keyData.push(key12, galList[j].id)
	    j-=1

       #galList[i].keyData.toString()

  subID = 1
  while  numData>1:
                  
                  
                  
		  maxKeyHeap = maxHeap()
		  
		  for i in range(1, numData):
		        if galList[i].keyData.size != 0:
			   # keyData is also a heap itself
			   maxKeyHeap.push(galList[i].keyData.peek().key, [i, galList[i].keyData.peek().ID])
		  
		  
		  if maxKeyHeap.size == 0: 
		      angle_limit = 1000
		      for i in range(1, numData):
		          key12 = pairKey(galList[i], newGalaxy)
		          galList[i].keyData.push(key12, newGalaxy.id)
		      continue

                  ID00 = maxKeyHeap.peek().ID[0]
		  ID11 = 0
		  for i in range(1, numData):
		      if galList[i].id == maxKeyHeap.peek().ID[1]: ID11 = i
		        
		  
		  
		  
		  
		  ID00_dist = ID00
		  ID11_dist = ID11
		  rank = 0
		  rank_dist = 0
		  
		  
		  sanitary = True
		  sanity_max = 10000
		  sanity_max_dist = 10000
		  no_pop = 0
		  
		  
		  ID0 = ID00
		  ID1 = ID11
		  
		  
###########################################################3
		  
		  while sanitary :

		    ID0 = maxKeyHeap.peek().ID[0]
		    ID1 = 0
		    for i in range(1, numData):
		        if galList[i].id == maxKeyHeap.peek().ID[1]: ID1 = i
		    
		    v0 = galList[ID0].Vls
		    v1 = galList[ID1].Vls
		    
		    Dist0 = galList[ID0].mDist
		    Dist1 = galList[ID1].mDist
		    eDist0 = Dist0 * galList[ID0].mDistErr
		    eDist1 = Dist1 * galList[ID1].mDistErr


		    sig_p0 = (galList[ID0].M_v2 / (2.0E6))**(1./3)
		    sig_p1 = (galList[ID1].M_v2 / (2.0E6))**(1./3)
		    ang12 = angle(galList[ID0].l,galList[ID0].b,galList[ID1].l,galList[ID1].b)*180./pi
		    sanity = (abs(v1-v0)) / (2.5*max(sig_p0,sig_p1))
		    
		    r1 = galList[ID0].R_theta
		    r2 = galList[ID1].R_theta
		    
		    if eDist0*eDist1 != 0 and Dist0*Dist1 != 0 :
		      if ang12 <= 0.7*(r1+r2):
		        sanityDist = abs(Dist0-Dist1) / (eDist0+eDist1)  #  
		      elif ang12 <= 1.05*(r1+r2):
			sanityDist = abs(Dist0-Dist1) / (1.0*sqrt(eDist0**2+eDist1**2))
		      else:
			sanityDist = abs(Dist0-Dist1) / max(eDist0,eDist1)
		    else:
		        sanityDist = 0.    
		    
		    if sanityDist!=0 and sanityDist <= 1.0:
		      sanitary = False
		    elif sanity <=1 and sanityDist == 0.:
		      sanitary = False
		    else:
		      
		      if sanityDist > 1.0:
			#print "sanity dist problem \n",
			if galList[ID0].level == 0: 
			  #print "pgc ID0: ", galList[ID0].id
			  mayI = True
			  for ip in range(0,len(list_dis_problem)):
			    if galList[ID0].id == list_dis_problem[ip][0]: mayI=False
			  if mayI: list_dis_problem.append([galList[ID0].id, galList[ID1].id])
			if galList[ID1].level == 0: 
			  #print "pgc ID1: ", galList[ID1].id
			  mayI = True
			  for ip in range(0,len(list_dis_problem)):
			    if galList[ID1].id == list_dis_problem[ip][0]: mayI=False
			  if mayI: list_dis_problem.append([galList[ID1].id, galList[ID0].id])
			#print "********** \n\n",
		      
		      
		      no_pop+=1
		      ###keeping track of the best velocity match
		      if sanity<=1 and sanityDist <= sanity_max_dist and sanityDist > 0:
			sanity_max_dist = sanityDist
			ID00_dist = ID0
			ID11_dist = ID1
			rank_dist = no_pop
			
		      if sanity <= sanity_max and sanityDist == 0:
			sanity_max = sanity
			ID00 = ID0
			ID11 = ID11
			rank = no_pop
			
		      if maxKeyHeap.size == 1:
			if rank_dist<=rank:
			   ID0 = ID00_dist
			   ID1 = ID11_dist
			else:
			   ID0 = ID00
			   ID1 = ID11
			sanitary = False
		      else:
			   temp = maxKeyHeap.pop()
			   ind = temp.ID[0]
			   galList[ind].keyData.pop()
			   if galList[ind].keyData.size > 0: 
			       maxKeyHeap.push(galList[ind].keyData.peek().key, [ind, galList[ind].keyData.peek().ID])

		       
###########################################################3
		  gal0 = galList[ID0].id
		  gal1 = galList[ID1].id
		  

		  
		  newGalaxy = galJoint(galList[ID0], galList[ID1], grID + subID)
		  
		  for i in range(ID0+1, numData):
		    galList[i].keyData.remove(galList[ID0].id)
		  for i in range(ID1+1, numData):
		    galList[i].keyData.remove(galList[ID1].id)    
		  
		  #if gal0 == 35043 or gal1 == 35043:
		  print "Process No.:", len(galList), newGalaxy.id, '    Serious pop No.: ', no_pop#, gal0, gal1
		  
		  galList.pop(ID0)
		  if ID0>ID1:
		    galList.pop(ID1)
		  else:
		    galList.pop(ID1-1)
		    
		  
		  
		  galList.insert(0, newGalaxy)
		  numData = len(galList)
		  galList[0].keyData = maxHeap()
		  for i in range(1, numData):
		            if angle(galList[i].l,galList[i].b,newGalaxy.l,newGalaxy.b)*180./pi < angle_limit :
			         key12 = pairKey(galList[i], newGalaxy)
			         galList[i].keyData.push(key12, newGalaxy.id)
	          subID += 1

  
  print "Process No.:", len(galList)
  
  print "Time: ", time()-t0
  #galList[0].toString()  # Root of tree
  
   
  
  return galList[0], list_dis_problem
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
def appendLeaf(root, leafList) :
  
  
  if root.left == None: 
    leafList.append([root.id, root.l, root.b])
  else: 
    appendLeaf(root.left, leafList)
    appendLeaf(root.right, leafList)
  
  return leafList
  
def allLeaves(root) :
  
  leafList = []
  appendLeaf(root, leafList)
  
  id = []
  sgl = []
  sgb = []
  
  
  for i in range(0, len(leafList)):
    id.append(leafList[i][0])
    sgl.append(leafList[i][1])
    sgb.append(leafList[i][2])
  
  return id, sgl, sgb
#################################################################
#################################################################
#################################################################
#################################################################
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

  plt.plot(X, Y, '-', markersize = 1)
  
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
  ##ax = fig.add_axes([0.13, 0.1, 0.85,  0.85]) 
  plt.ylim(-100,100)
  plt.xlim(380,-20)
  #plt.xlabel("SGL (deg)", fontsize=20)
  #plt.ylabel("SGB (deg)", fontsize=20)
  plt.yticks(fontsize=14)
  plt.xticks(fontsize=14)
  
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
           x = r*np.cos(u[q]) + G_list[i][0].l
           y = r*np.sin(u[q]) + G_list[i][0].b
           X[q], Y[q] = xymap(x,y, l0)

	
	if r <= 6.3:
	  line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Blue, Green)) 
	  line.set_dashes([8, 3]) 
	
	
	
	
	X  = []
	Y = []
	

	

	for j in range(1, len(G_list[i])):

	    x, y = xymap(G_list[i][j].l,G_list[i][j].b, l0)
	    X.append(x)
	    Y.append(y)
	    
	plt.plot(X, Y, 'o', markersize = 3, color=(Red, Blue, Green), markeredgecolor = 'none')
        

#################################################################
#################################################################
#################################################################
def dend_setX(root, n):
  
  if root.left == None: 
    root.dend_x = n
    return n+1
  else:
    l = dend_setX(root.left, n)
    r = dend_setX(root.right, l)
    root.dend_x = (root.left.dend_x+root.right.dend_x)/2.
    return r
 
#################################################################
def addGalaxyTable(galaxy, myTable, G_members, subGalaxies):
  
    pgc = galaxy.id 
    gl = galaxy.gl  
    gb = galaxy.gb      
    sgl = galaxy.l  
    sgb = galaxy.b  
    Vls = galaxy.Vls  
    logK = galaxy.logK  
    Ty = galaxy.Ty  
    dcf2 = galaxy.dcf2Copy
    ed = galaxy.edCopy
    
    mDist = galaxy.mDist
    mDistErr = galaxy.mDistErr
    
    level = galaxy.level  
    R_theta = galaxy.R_theta  
    sigma = galaxy.sigma  
    v_av = galaxy.v_av  
    v2_av = galaxy.v2_av  
    nest = galaxy.nest  
    
    myTable.add_row([pgc, gl, gb, sgl, sgb, Vls, logK, Ty, dcf2, ed, mDist, mDistErr, G_members, subGalaxies, level \
			 , R_theta, sigma, v_av, v2_av, nest])

#################################################################
def addGroupTable(root, myTable):
  
  n = root.subGalaxies
  addGalaxyTable(root, myTable, n, n)
  if root.right != None:
    addGroupTableCore(root.right, myTable, n, root.subGalaxies)
    addGroupTableCore(root.left, myTable, n, root.subGalaxies)
  
def addGroupTableCore(root, myTable, G_members, subGalaxies):
  
  if root.right != None:
    
    if root.left.right != None and root.right.dend_x - root.left.right.dend_x > 5:
	addGalaxyTable(root.right, myTable, G_members, root.right.subGalaxies)
        addGroupTableCore(root.right, myTable, G_members, root.right.subGalaxies)
        addGalaxyTable(root.left, myTable, G_members, root.left.subGalaxies)
        addGroupTableCore(root.left, myTable, G_members, root.left.subGalaxies)
    else:    
        addGroupTableCore(root.right, myTable, G_members, subGalaxies)
        addGroupTableCore(root.left, myTable, G_members, subGalaxies)
  
  if root.right == None:
    addGalaxyTable(root, myTable, G_members, subGalaxies)

  
  
  
  
  

################################################################# 
def groupWrite(outfile, G_list, root):
  
  
  # finding the leave nodes
  galList = NodeLeaf(root)
  galList.pop(0)
  
  
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
        sgl = galaxy.l  
        sgb = galaxy.b  
        Vls = galaxy.Vls  
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
       
 
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,Ks,logK,Vls, dcf2, ed, \
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_lum, \
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, tX_dyn, tX_lum, subGalaxies, nest, \
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname])
	

  
  # writing individual galaxies
  for galaxy in galList:


      if galaxy.inGroup == 0.:
	flag = 0
	galaxy.flag = 0
        pgc = galaxy.id  
        gl = galaxy.gl  
        gb = galaxy.gb
        sgl = galaxy.l  
        sgb = galaxy.b  
        Vls = galaxy.Vls  
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
       
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,Ks,logK,Vls, dcf2, ed, \
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_lum, \
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, tX_dyn, tX_lum, subGalaxies, nest, \
	            coordinate_src, Ty_src , Ks_src, Vls_src, objname])

	

  
  #myTable.write(outfile, format='ascii',delimiter=',')
  
  myTable.write(outfile, format='ascii.fixed_width',delimiter='|', bookend=False)
################################################################# 
def R_node(root, fraction):
  
  if root.left == None: return 0
  
  galList = NodeLeaf(root)
  
  L_tot = 10**root.logK
  L_lim = fraction * L_tot
  
  l0 = root.l
  b0 = root.b
  
  N = len(galList)
  theta = np.zeros(N-1) 
  logK = np.zeros(N-1) 
  Velocity = root.v_av
  
  for i in range(1,N):
      
      theta[i-1] = (angle(l0,b0, galList[i].l, galList[i].b))
      logK[i-1]  = (m_logK(galList[i].Ks, galList[i].Vls, Velocity, galList[i].dcf2))
  
  indices = np.argsort(theta)
  theta   =  theta[indices]
  logK  =    logK[indices]
  
  i = 0
  L_tot_radial = 0
  while L_lim >= L_tot_radial and i < N-1:
    L_tot_radial +=  10**logK[i]
    i+=1
  
  f = L_tot_radial / L_tot
  angular_radius = theta[i-1]
  
  
  for alfa in np.arange(0,1.0001,0.0001):
      if f <= 0.75*(alfa**2)*(1+sqrt(1-alfa**2))+0.25*((1-sqrt(1-alfa**2))**3): break
  
  radius = Velocity*tan(angular_radius/alfa)/H0
  
  return radius


################################################################# 


################################################################# 

def Theta_max(root):
  
  if root.left == None: return 0
  
  galList = NodeLeaf(root)
  
  N = len(galList)
  theta = np.zeros(N-1) 
  
  for i in range(1,N):
      
      theta[i-1] = angle(root.l, root.b, galList[i].l, galList[i].b)
  
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
      
      distij = angle(galList[i].l, galList[i].b, galList[j].l, galList[j].b)
      if distij !=0 : 
         sum += 1. / distij
      else: BOL = False    # If one of the didtnaces is zero, therre is change to have duplicates 
      
  if BOL == True and sum != 0:
    Rg = ((N)*(N)) / sum
  
  return Rg




################################################################# 
################################################################# 

def fornaxPlot(G_list, galList):
  
    
  
    #print table.dtype
    #id   = table['pgc']
    #gl  = table['gl']
    #gb  = table['gb']
    #dcf2 =  table['dcf2']
    #Vls  =  table['Vls']
    
    id   = []
    gl   = []
    gb   = []
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
  
    #col = 0
    
    
    
    #RA0 = []
    #DEC0 = []
    #for i in range(0, N_galaxies):
      #point = coord.Galactic(gl[i] , gb[i], unit=(unit.degree, unit.degree))
      #RA0.append(point.fk5.ra.degree)
      #DEC0.append(point.fk5.dec.degree)
    
    
    fig = plt.figure(figsize=(7, 7), dpi=100)
    
    #fig.add_axes([left,bottom,width,height]) 
    
    ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))  
  
  
    plt.minorticks_on()
    plt.tick_params(which='major', length=7, width=1.5)
    plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     
    
    plt.plot(RA0, DEC0, '.', markersize = 1, color='#696969')  # gray    color='#696969'
    plt.ylim(-45,-10)
    plt.xlim(70,35)
    plt.xlabel("RA (deg)", fontsize=20)
    plt.ylabel("DEC (deg)", fontsize=20)
    
    for i in range(0, NoGroups):  # for all groups
  
      print "ploting group #: ", i+1
      Key = True
      if Key:
	  random.seed(G_list[i][0].nest )
	  Red, Blue, Green = random.random(), random.random(), random.random()
	  
	  
	  #r = G_list[i][0].R_theta
	  
	  Dist = G_list[i][0].mDist
	  if Dist ==0 : 
	     Dist = G_list[i][0].Vls/H0
	  r = 180.*atan(G_list[i][0].R_2t2/Dist)/pi
	  #r = 180.*Theta_max(G_list[i][0])/pi
	  
	  
	  
	  d_theta = 0.01
	  u = np.arange(0,2*pi,d_theta)
	  X = u*0.
	  Y = u*0.
	  #point = coord.Galactic(G_list[i][0].gl , G_list[i][0].gb, unit=(unit.degree, unit.degree))
	  #RA0 = point.fk5.ra.degree
	  #DEC0 = point.fk5.dec.degree
	  RA0 = G_list[i][0].ra
	  DEC0 = G_list[i][0].dec
	  
	  for q in range(0,len(u)):
	    x = r*np.cos(u[q]) + RA0
	    y = r*np.sin(u[q]) + DEC0
	    X[q] = x
	    Y[q] = y

	  
	  if r <= 10000:
	    line, = plt.plot(X,Y, '-', markersize = 2, color=(Red, Blue, Green)) 
	    line.set_dashes([8, 3]) 
	    plt.text(RA0, DEC0+r, int(G_list[i][0].nest), fontsize=8, color=(Red, Blue, Green))
	  
	  
	  
	  X  = []
	  Y = []
	  

	  

	  for j in range(1, len(G_list[i])):

	      point = coord.Galactic(G_list[i][j].gl , G_list[i][j].gb, unit=(unit.degree, unit.degree))
	      RA0 = point.fk5.ra.degree
	      DEC0 = point.fk5.dec.degree
	      X.append(RA0)
	      Y.append(DEC0)
	      
	  plt.plot(X, Y, 'o', markersize = 3, color=(Red, Blue, Green), markeredgecolor = 'none')
	  

################################################################# 
################################################################# 
#################################################################
#################################################################

def virgoPlot(table, G_list, R_min, R_max):
  
  id   = table['pgc']
  sgl  = table['sgl']
  sgb  = table['sgb']
  Vls  = table['Vls']
  flag = np.where(Vls>0)
  sgl = sgl[flag]
  sgb = sgb[flag]
  Vls = Vls[flag]
  flag = np.where(Vls<=3500)
  sgl = sgl[flag]
  sgb = sgb[flag]
  Vls = Vls[flag]
  
  
  
  fig = plt.figure(figsize=(7, 7), dpi=100)
  ax = fig.add_axes([0.12, 0.1, 0.85,  0.85]) 
  ax.xaxis.set_major_locator(MultipleLocator(10))
  ax.yaxis.set_major_locator(MultipleLocator(10))
  ax.xaxis.set_minor_locator(MultipleLocator(1))
  ax.yaxis.set_minor_locator(MultipleLocator(1))  
  
  
  
  plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)  
  
  
  xmin = 5*floor((sglV+35)/5.)+3
  xmax = 5*ceil((sglV-35)/5.)-3
  #ymax = 5*floor((sgbV+35)/5.)+3
  ymin = 5*ceil((sgbV-35)/5.)-3
  ymax = ymin - (xmax-xmin)
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
  
  red = 620  # nm
  violet = 380 # nm
  NoGroups = len(G_list)
  print "VirgoPlot >> Number of groups: ", NoGroups
  #plt.plot(sgl, sgb, '.', markersize = 3, color='black')  # gray    color='#696969'
  plt.plot(sgl, sgb, '.', markersize = 1, color='#696969')  # gray    color='#696969'
  
  col = 0
  colors= ['#0000FF','#CC00CC','#999900','#006600','#FF0000','#FF8000', '#6600CC'] 
  for i in range(0, NoGroups):  # for all groups
  
    groupFromVirgo = angle(G_list[i][0].l, G_list[i][0].b, sglV, sgbV)*180./pi
    
    if groupFromVirgo>=6.8 and groupFromVirgo<=1000:

	random.seed(G_list[i][0].nest )
	Red, Blue, Green = random.random(), random.random(), random.random()
	
	
	#r = G_list[i][0].R_theta
        Dist = G_list[i][0].mDist
	if Dist==0.: Dist = G_list[i][0].Vls/H0
	r = 180.*atan(G_list[i][0].R_2t2/Dist)/pi
	#r = 180.*Theta_max(G_list[i][0])/pi
	
	d_theta = 0.001
	theta = np.arange(0,2*pi,d_theta)
	Circlx = r*np.cos(theta) + G_list[i][0].l
	Circly = r*np.sin(theta) + G_list[i][0].b
	
	if r <= 7:
	  line, = plt.plot(Circlx, Circly, '-', markersize = 2, color=(Red, Blue, Green)) 
	  line.set_dashes([8, 3]) 
	
	
	
	sRA  = []
	sDEC = []
	


	for j in range(1, len(G_list[i])):
 
	    sRA.append(G_list[i][j].l)
	    sDEC.append(G_list[i][j].b)
	    
	plt.plot(sRA, sDEC, 'o', markersize = 4, color=(Red, Blue, Green), markeredgecolor = 'none')
        

#################################################################
def mergGroup(G_list, restrict=False):
  
    NoGroups = len(G_list)
    print "Merge  >> N.Groups: ", NoGroups
    

    
    i = 0
    while i < NoGroups:
      j = i+1
      while j < NoGroups-1:
	
	#r1 = 180.*G_list[i][0].R_theta/pi
	#r2 = 180.*G_list[j][0].R_theta/pi#Theta_max(G_list[j][0])
	

	
	if G_list[i][0].mDist != 0 and G_list[j][0].mDist != 0 :
	  r1 = (180.*atan(G_list[i][0].R_2t2/(G_list[i][0].mDist))/pi)
	  r2 = (180.*atan(G_list[j][0].R_2t2/(G_list[j][0].mDist))/pi)
	else: 
	  r1 = (180.*atan(G_list[i][0].R_2t2/(G_list[i][0].Vls/H0))/pi)
	  r2 = (180.*atan(G_list[j][0].R_2t2/(G_list[j][0].Vls/H0))/pi)
	




	
	ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0].gb, G_list[j][0].gl, G_list[j][0].gb)
	
	
	d1 = G_list[i][0].mDist
	e1 = d1 * G_list[i][0].mDistErr
	d2 = G_list[j][0].mDist
	e2 = d2 * G_list[j][0].mDistErr
	delt = abs(d1-d2)
	
	v1 = G_list[i][0].Vls
	v2 = G_list[j][0].Vls
	#sig1 = G_list[i][0].sigma
	#sig2 = G_list[j][0].sigma
	
	
	# using dynamical velocity dispersions calculated based on the luminosities
        sig1 = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
	sig2 = (G_list[j][0].M_v2 / (2.0E6))**(1./3)
	
	
	
	
	n1 = G_list[i][0].subGalaxies
	n2 = G_list[j][0].subGalaxies
	#print ang12
	Bol = False
	


#######################################################################
        Sigquad = sqrt(sig1**2+sig2**2)
#######################################################################
        if restrict:
	  if ang12 <= 1.05*abs(G_list[i][0].R_theta-G_list[j][0].R_theta) and max(r1,r2)<6:
	    if d1!=0 and d2!=0 and delt <= 2.0*(e1+e2) and abs(v1-v2) <= 5.0*max(sig1,sig2):
		Bol = True
	    elif (d1==0 or d2==0) and abs(v1-v2) <= 5.0*max(sig1,sig2):
		Bol = True



#######################################################################
        
	if ang12 <= 0.2*(r1+r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt <= 2.0*(e1+e2) and abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True
	  elif (d1==0 or d2==0) and abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True

	      
	if ang12 <= 0.4*(r1+r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt <= 1.4*(e1+e2) and abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True
	  elif (d1==0 or d2==0) and abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True
	      
#######################################################################

	if ang12 <= 0.7*(r1+r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt <= (e1+e2) and abs(v1-v2) <= 1.*max(sig1,sig2):
	      Bol = True
	  elif (d1==0 or d2==0) and abs(v1-v2) <= 1.*max(sig1,sig2):
	      Bol = True


	if ang12 <= 0.8*(r1+r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt <= (e1+e2) and abs(v1-v2) <= 1.*max(sig1,sig2):
	      Bol = True
	  elif (d1==0 or d2==0) and abs(v1-v2) <= 1.*max(sig1,sig2):
	      Bol = True


	if ang12 <= 1.00*(r1+r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt <= 1.0*max(e1,e2) and abs(v1-v2) <= 2.*max(sig1,sig2):
	      Bol = True
	  #elif (d1==0 or d2==0) and abs(v1-v2) <= 1.0*min(sig1,sig2):
	      #Bol = True
	      
	      

#######################################################################


        if Bol:
	  
	  #print "YES"
	  newGroup = []
	  
	  for p in range(1, len(G_list[i])):
	    newGroup.append(G_list[i][p])
	  
	  for q in range(1, len(G_list[j])):
	    newGroup.append(G_list[j][q])
	  
	  #newGalaxy = galJoint(newGroup[0], newGroup[1], 10000000)
	  #p=2
	  #while p < len(newGroup):
	    #newGalaxy = galJoint(newGalaxy, newGroup[p], 10000000+p)
	    #p+=1
	  
	  
	  root = LoosgalListJoint(newGroup, grID = 2000000000)
	  #root = newGalaxy
	  #newGroup.insert(0, newGalaxy)
	  
	  G_list[i] = NodeLeaf(root)
	  #G_list[i] = newGroup
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
  
    print "AddGal >> N.Groups: ", len(G_list)
    # finding the leave nodes
    
    modified_G_list = []
    for group in G_list:
      if group[0].R_theta*180./pi > 10.:
	for galaxies in group[1:]:
	  galaxies.inGroup = 0
      else:
	modified_G_list.append(group)
    G_list = modified_G_list
	
      
    
    
    
    
    #galList = NodeLeaf(root)
    #galList.pop(0)
    
    
    LL=[]
    BB=[]
    N = len(galList)
    for i in range(0,N):
       LL.append(galList[i].gl)
       BB.append(galList[i].gb)
    LL = np.asarray(LL)
    BB = np.asarray(BB)
    
    
    # Lum_mass of all groups
    Lum_mass = []
    for i in range(0, len(G_list)):
        Lum_mass.append(G_list[i][0].M_v2)
    
    
    # We first start with larger Lum_mass groups
    Lum_mass = np.asarray(Lum_mass)
    indices = np.argsort(Lum_mass)
    

    #subgalList = galList
    for i in indices[::-1]:
            ChangeGroup = False
            p=0
	    r = 2.*(180.*atan(G_list[i][0].R_2t2*1.05*sqrt(1.5)/(G_list[i][0].Vls/H0))/pi)
	    
	    q1 = np.zeros((N,), dtype=np.int)
            q2 = np.zeros((N,), dtype=np.int)
            q3 = np.zeros((N,), dtype=np.int)
            q4 = np.zeros((N,), dtype=np.int)
	    q1[np.where(LL<=G_list[i][0].gl+3*r)] = 1
	    q2[np.where(LL>=G_list[i][0].gl-3*r)] = 1
	    q3[np.where(BB<=G_list[i][0].gb+3*r)] = 1
	    q4[np.where(BB>=G_list[i][0].gb-3*r)] = 1
	    qq = q1+q2+q3+q4
	    subgalList =[]
	    indices = np.where(qq==4)
	    for index in indices[0]:
	       subgalList.append(galList[index])
	    
	    while p < len(subgalList):

	      ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0].gb, subgalList[p].gl, subgalList[p].gb)
	      
	      
	      # if the galaxy is close enough to the center of the group
	      # and it's not already in any other group, then it would be added to the current group
	      if ang12<=10 and ang12 <= r:
		Bol = True
		if  subgalList[p].inGroup == 1:
		          Bol = False  #  It's already got groupped
	        
	        sig_p = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
	        Dist0 = subgalList[p].dcf2
		Dist1 = G_list[i][0].mDist
		eDist0 = Dist0 * subgalList[p].ed
		eDist1 = Dist1 * G_list[i][0].mDistErr
		
		if eDist0*eDist1 != 0 and Dist0*Dist1 != 0:
		    sanityDist = abs(Dist0-Dist1) / sqrt(eDist0**2+eDist1**2)
		else:
		    sanityDist = 0.
		if Bol==True and subgalList[p].dcf2==0 and Dist1 !=0 and ang12 <= 1.0*(180.*atan(G_list[i][0].R_2t2/(Dist1-eDist1))/pi) and abs(subgalList[p].Vls - G_list[i][0].Vls) < 2.*sig_p:
		  #print "Yaaay", subgalList[p].id, G_list[i][0].id
		  ChangeGroup = True
		  #G_list[i].pop(0)
		  
		  # Tha galaxy is absorbed
		  subgalList[p].inGroup = 1
		  G_list[i].append(subgalList[p])
		  subgalList.pop(p)
		  p-=1
	        elif Bol==True and subgalList[p].dcf2==0 and Dist1 ==0  and ang12 <= 1.0*(180.*atan(G_list[i][0].R_2t2/(G_list[i][0].Vls/H0))/pi) and abs(subgalList[p].Vls - G_list[i][0].Vls) < 2.*sig_p:
		  
		  ChangeGroup = True
		  
		  
		  # Tha galaxy is absorbed
		  subgalList[p].inGroup = 1
		  G_list[i].append(subgalList[p])
		  subgalList.pop(p)
		  p-=1
		
	        elif Bol==True and sanityDist > 0. and sanityDist<=1 and ang12 <= 1.0*(180.*atan(G_list[i][0].R_2t2/(Dist1-eDist1))/pi):
		  #print "Yaaay.Dist", subgalList[p].id, G_list[i][0].id
		  ChangeGroup = True
		  #G_list[i].pop(0)
		  
		  # Tha galaxy is absorbed
		  subgalList[p].inGroup = 1
		  G_list[i].append(subgalList[p])
		  subgalList.pop(p)
		  p-=1
              p+=1
              # end while
              
              
              
            if ChangeGroup:
	      # remove the root from the beginning of the list
	      G_list[i].pop(0)
	      # Group rearrangement, make an hirarchical binary tree
              root = LoosgalListJoint(G_list[i], grID = 3000000000)
              G_list[i] = NodeLeaf(root) 
    return G_list

#################################################################
def IndividualGal(galList):
  
  

  #galList = NodeLeaf(root)
  #galList.pop(0)
  
  
  
  p = 0 
  N = len(galList)
  while p < N:
    if galList[p].inGroup == 1: 
	   galList.pop(p)
	   p-=1
	   N-=1
    p+=1
  

  #p = 0 
  #while p < len(galList):
    #for i in range(0, len(G_list)):
      #for j in range(1, len(G_list[i])):
         #if galList[p].id == G_list[i][j].id:
	   #galList.pop(p)
	   #p-=1
    #p+=1




  
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
    individuals.pop(p)
    N-=1
    q = p
    piv = (p-1)%N
    grp = galaxies[0]
    Bol = False
    while q!=piv:
      
      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].gl, individuals[q].gb)
      if ang12 <= (180.*atan(grp.R_2t2/(grp.Vls/H0))/pi):
	sig_p = (grp.M_v2 / (2.0E6))**(1./3)
	Dist0 = grp.mDist
	Dist1 = individuals[q].dcf2
	eDist0 = Dist0 * grp.mDistErr
	eDist1 = Dist1 * individuals[q].ed
	if eDist0*eDist1 != 0 and Dist0*Dist1 != 0:
	    sanityDist = abs(Dist0-Dist1) / (eDist0+eDist1)
	else:
	    sanityDist = 0.
	if (sanityDist == 0 and abs(grp.Vls-individuals[q].Vls) <= 2.*sig_p) or (sanityDist > 0. and sanityDist<=1):
	    print " A new group is developing !!!"
	    Bol = True
	    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 4000000000)
	    individuals.pop(q)
	    N-=1
	    piv = q 
	    if p>=q:  p-=1
	    q-=1
	    
      q = (q+1)%N
    if Bol: #len(galaxies) > 1:
      NEWgpLIST.append(NodeLeaf(grp))
    p+=1
    
    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1
    
  return NEWgpLIST

#################################################################

def LoosgalListJoint(galList, grID = 5000000000):
   
   if len(galList)==0: return None
 
   Lum_mass = []
   for i in range(0, len(galList)):
        Lum_mass.append(galList[i].Ks)
   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass)
   
   root = None
   for i in indices:
       root = galJoint(root, galList[i], grID + galList[i].id)
   
   return root
   #return galListJoint(galList)






#################################################################

if __name__ == '__main__':
  
  
  R_min = 6.0
  R_max = 50.
  n_gal = 500000
  
  cluster = 'virgo'
  #cluster = 'fornax'
  
  #cluster  = 'north'
  #cluster = 'south'
  
  
  
  
  if cluster == 'virgo':
    inFile  = 'AllSky.north.v7.csv'
  elif cluster == 'fornax':
    inFile  = 'AllSky.south.v7.csv'
  elif cluster == 'north':
    inFile  = 'AllSky.north.v7.csv'
  else:
    inFile  = 'AllSky.south.v7.csv'   
  
  table = np.genfromtxt( inFile , delimiter=',', filling_values=0, names=True, dtype=None)
  
  
  
  
  
  if   cluster == 'virgo' :
     galList = readgalList(table, R_min, R_max, n_gal, useDCF2=True)
  elif cluster == 'fornax': 
     galList =  Fornax_readgalList(table, R_min, R_max, n_gal, useDCF2=True)
  elif cluster == 'north' :
     R_max = 100000000.
     galList = readgalList(table, R_min, R_max, n_gal, useDCF2=True)
  else:
     R_max = 100000000.
     galList = readgalList(table, R_min, R_max, n_gal, useDCF2=True)    
  
  

  
  
  
  ####### Continue the iteration
  #ML = 60
  #if cluster == 'virgo':
    #root = readTree('virgo.ML.'+str(int(ML))+'.v09.tree')
  #elif cluster == 'fornax':
    #root = readTree('fornax.ML.'+str(int(ML))+'.v09.tree')
  #elif cluster == 'north':
    #root   = readTree('north.ML.'+str(int(ML))+'.v09.tree')
  #else:
    #root = readTree('south.ML.'+str(int(ML))+'.v09.tree')
  #print "List making ..."
  #galList = NodeLeaf(root)
  #galList.pop(0)    # removing the root from the beginning of the list
  #print "Velocity modifying ..."
  #for i in range(0,len(galList)):
    #galList[i].v_av  = galList[i].Vls
    #galList[i].v2_av = galList[i].Vls**2
   
  for ML in [60]:
  #for ML in [60,61,62,63,64]:


	#root, badDist = galListJoint(galList)
	#if cluster == 'virgo':
	  #ehsan = Table()
	  #ehsan.add_column(Column(data=badDist,name='pgc'))
	  #ehsan.write('virgo_dist_error.ML.'+str(int(ML))+'.v09.txt', format='ascii',delimiter=',')
	#elif cluster == 'fornax':
	  #ehsan = Table()
	  #ehsan.add_column(Column(data=badDist,name='pgc'))
	  #ehsan.write('fornax_dist_error.ML.'+str(int(ML))+'.v09.txt', format='ascii',delimiter=',')
	  
	
	#if cluster == 'virgo':
	  #writeTree('virgo.ML.'+str(int(ML))+'.v09.tree', root) # write the tree
	#elif cluster == 'fornax':
	  #writeTree('fornax.ML.'+str(int(ML))+'.v09.tree', root) # write the tree
        #elif cluster == 'north':
          #writeTree('north.ML.'+str(int(ML))+'.v09.tree', root) # write the tree
        #else:
	  #writeTree('south.ML.'+str(int(ML))+'.v09.tree', root) # write the tree
	
	
	##if cluster == 'virgo':
	  ##root = readTree('virgo.ML.'+str(int(ML))+'.v09.tree')
	##elif cluster == 'fornax':
	  ##root = readTree('fornax.ML.'+str(int(ML))+'.v09.tree')
        ##else:
	  ##root   = readTree('north.ML.'+str(int(ML))+'.v09.tree')
	  ##root_1 = readTree('south.ML.'+str(int(ML))+'.v09.tree')



	##  Finding groups based on the criteria
	#G_list = groupID(root, ML)
	#G_list_1 = groupID(root_1, ML)


        #for pq in range(0,10): 
	  #G_list = addGalGroup(G_list, root)
	
        #for pq in range(0,10): 
	  #G_list = addGalGroup(G_list, root)	
	  #mergGroup(G_list)
	  #G_list = addGalGroup(G_list, root)	
  

        
        G_list = []
        for pq in range(0,20): 
           IndivG_lis = IndividualGal(galList)
           G_list = G_list + IndivG_lis  
        for pq in range(0,10): 
          mergGroup(G_list)
        
        # Don't comment it out
        for pq in range(0,10): 
          G_list = addGalGroup(G_list, galList)	
        for pq in range(0,10): 
          mergGroup(G_list)
        G_list = addGalGroup(G_list, galList, restrict=True)	





	
	
	for group in G_list:
	  groupHeader = group[0]
	  groupHeader.flag = 2
	  meanDist = groupHeader.mDist
	  meanDistErr = groupHeader.mDistErr
	  L_tot = 0
	  ID = groupHeader.id 
	  groupHeader.id = ID - ID % 1000000000 + groupHeader.nest
	  for galaxy in group[1:]:
	    galaxy.flag = 1
	    # and also modify the absolute luminosities
	    galaxy.setMeanDist(meanDist, meanDistErr, GRP_vel = groupHeader.Vls)
	    L_tot += 10**galaxy.logK
	 
          groupHeader.logK = log10(L_tot)
          if meanDist == 0 :  meanDist = groupHeader.Vls / H0
          Mk_sun = 3.28   # K-band
          M = Mk_sun - (groupHeader.logK / 0.4)
          groupHeader.Ks = M + 5*log10(meanDist) + 30 - 5
          
          
          
          







        #groupPlot(G_list, G_list_1)
        
        #galList = NodeLeaf(root)
	#galList.pop(0) 
        if cluster == 'virgo':
	  virgoPlot(table, G_list, R_min, R_max)
	else: fornaxPlot(G_list, galList)
	






        #if cluster == 'virgo':
	  #groupWrite('virgo.ML.'+str(int(ML))+'.v09.group', G_list, root)
        #elif cluster == 'fornax': 
	  #groupWrite('fornax.ML.'+str(int(ML))+'.v09.group', G_list, root)
        #elif cluster == 'north':
          #groupWrite('north.ML.'+str(int(ML))+'.v09.group', G_list, root)
        #else:
	  #groupWrite('south.ML.'+str(int(ML))+'.v09.group', G_list, root)

   

	#galList = NodeLeaf(root)
	#galList.pop(0)    # removing the root from the beginning of the list
	for galaxy in galList: galaxy.inGroup = 0
        
        
        #if cluster == 'virgo':
           #plt.savefig('virgo.ML.'+str(int(ML))+'.v09.eps', dpi=600)  # the best view and size happens in eps format
        #elif cluster == 'fornax':
	   #plt.savefig('fornax.ML.'+str(int(ML))+'.v09.eps', dpi=600)  # the best view and size happens in eps format
  
  

  plt.show()
  