#!/usr/bin/python						﻿#!/usr/bin/python
__author__ = "Ehsan Kourkchi"					__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2016"				__copyright__ = "Copyright 2016"
__credits__ = ["Ehsan Kourkchi"]				__credits__ = ["Ehsan Kourkchi"]
__version__ = "1.0"						__version__ = "1.0"
__maintainer__ = "Ehsan Kourkchi"				__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"				__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"					__status__ = "Production"
# As of April, 8, 2016						# As of April, 8, 2016
###								###
# Written by Ehsan Kourkchi (September 2015)			# Written by Ehsan Kourkchi (September 2015)
# email: ehsan@ifa.hawaii.edu					# email: ehsan@ifa.hawaii.edu
# This code, identifies groups of galaxies, given a 		# This code, identifies groups of galaxies, given a 
# a cataloge of galaxies. The radial velocity of galaxies wou	# a cataloge of galaxies. The radial velocity of galaxies wou
# be used to find galaxies with relatively the same radial ve	# be used to find galaxies with relatively the same radial ve
# and almost the same position on the sky			# and almost the same position on the sky
###								###

## Importing Important Python Libraries				## Importing Important Python Libraries
import sys							import sys
import os							import os
import random							import random
import matplotlib.pyplot as plt					import matplotlib.pyplot as plt
from matplotlib import rc, rcParams				from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrForma	from matplotlib.ticker import MultipleLocator, FormatStrForma
import numpy as np						import numpy as np
from math import *						from math import *
from time import time						from time import time
import wl_to_rgb as col						import wl_to_rgb as col
import random							import random
from astropy.io import ascii					from astropy.io import ascii
from astropy.table import Table, Column 			from astropy.table import Table, Column 
import pyfits							import pyfits
import pylab as py						import pylab as py

from astropy import coordinates as coord			from astropy import coordinates as coord
from astropy import units as unit				from astropy import units as unit
import subprocess						import subprocess

# **************************************			# **************************************
# Global Variables						# Global Variables
# Physical Constants						# Physical Constants
# **************************************			# **************************************
H0 = 75.           # hubble constant				H0 = 75.           # hubble constant
sglV = 102.8806    # M87 center - super galactic longitude	sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude	sgbV = -2.3479     # M87 center - super galactic latitude
G = 6.67384E-11   # Gravitational constant			G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]				M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.08567758E19						Mpc_km = 3.08567758E19
#t0_age = Mpc_km / H0   # Universe Age [sec]			#t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]				Yr= 365.25 * 24 * 3600    # [sec]
t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr		t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr
h = 0.75  # hubble constant					h = 0.75  # hubble constant
# **************************************			# **************************************

# Returns mass in solar unit, given the absolute K-band lumin	# Returns mass in solar unit, given the absolute K-band lumin
def Mass(L_k):							def Mass(L_k):

  L = L_k / 1.E10						  L = L_k / 1.E10
  								  
  if L < 1.:							  if L < 1.:
    #MtoL = 43.0#*(L**-0.3)					    #MtoL = 43.0#*(L**-0.3)
    MtoL = 32.0*(L**-0.5)					    MtoL = 32.0*(L**-0.5)
  elif L > 1000.:						  elif L > 1000.:
    MtoL = 91.						      |	    MtoL = 121.19
  else:								  else:
    MtoL = 32*(L**0.15)					      |	    MtoL = 43*(L**0.15)
  								  
  								  
  Mass_out = L_k * MtoL					      |	  Mass_out = h * L_k * MtoL
  								  
  return Mass_out						  return Mass_out


# **************************************			# **************************************
# returns angular separation of 				# returns angular separation of 
# two vectors in radian						# two vectors in radian
# **************************************			# **************************************
def angle(l1, b1, l2, b2):					def angle(l1, b1, l2, b2):
  								  
   cl1 = cos(l1*pi/180.)					   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)					   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)					   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)					   sb1 = sin(b1*pi/180.)
   								   
   x1 = cl1 * cb1						   x1 = cl1 * cb1
   y1 = sl1 * cb1						   y1 = sl1 * cb1
   z1 = sb1							   z1 = sb1
   								   
   cl2 = cos(l2*pi/180.)					   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)					   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)					   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)					   sb2 = sin(b2*pi/180.)
   								   
   x2 = cl2 * cb2						   x2 = cl2 * cb2
   y2 = sl2 * cb2						   y2 = sl2 * cb2
   z2 = sb2   							   z2 = sb2   
   								   
   XdotY = x1*x2 + y1*y2 + z1*z2				   XdotY = x1*x2 + y1*y2 + z1*z2
   #X2 = sqrt(x1**2 + y1**2 + z1**2)				   #X2 = sqrt(x1**2 + y1**2 + z1**2)
   #Y2 = sqrt(x2**2 + y2**2 + z2**2)				   #Y2 = sqrt(x2**2 + y2**2 + z2**2)
   								   
   if XdotY > 1 :						   if XdotY > 1 :
     theta12 = 0.						     theta12 = 0.
   elif XdotY < -1 :						   elif XdotY < -1 :
     theta12 = -1.*pi						     theta12 = -1.*pi
   else:							   else:
     theta12 = acos(XdotY)  					     theta12 = acos(XdotY)  
   return theta12   # radian					   return theta12   # radian
# **************************************			# **************************************
# returns sign of a number 					# returns sign of a number 
def sign(x):							def sign(x):
  if x<0 : return -1						  if x<0 : return -1
  if x>=0 : return 1						  if x>=0 : return 1
# **************************************			# **************************************
# L is luminosity						# L is luminosity
# l,b are galaxtic coordinates, and d is distance		# l,b are galaxtic coordinates, and d is distance
# returns the barycentric coordiantes of a pair 		# returns the barycentric coordiantes of a pair 
def barycenter(L1, l1, b1, d1, L2, l2, b2, d2):			def barycenter(L1, l1, b1, d1, L2, l2, b2, d2):
  								  
   if d1==0 or d2==0:						   if d1==0 or d2==0:
       								       
       dd1 = 1.							       dd1 = 1.
       dd2 = 1.							       dd2 = 1.
   else: 							   else: 
       dd1 = d1							       dd1 = d1
       dd2 = d2							       dd2 = d2
 								 
   								   
   cl1 = cos(l1*pi/180.)					   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)					   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)					   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)					   sb1 = sin(b1*pi/180.)
   								   
   x1 = dd1 * cl1 * cb1						   x1 = dd1 * cl1 * cb1
   y1 = dd1 * sl1 * cb1						   y1 = dd1 * sl1 * cb1
   z1 = dd1 * sb1						   z1 = dd1 * sb1
   								   
   cl2 = cos(l2*pi/180.)					   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)					   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)					   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)					   sb2 = sin(b2*pi/180.)
   								   
   x2 = dd2 * cl2 * cb2						   x2 = dd2 * cl2 * cb2
   y2 = dd2 * sl2 * cb2						   y2 = dd2 * sl2 * cb2
   z2 = dd2 * sb2   						   z2 = dd2 * sb2   


   L_tot = L1 + L2						   L_tot = L1 + L2
   x = (L1*x1+L2*x2)/L_tot					   x = (L1*x1+L2*x2)/L_tot
   y = (L1*y1+L2*y2)/L_tot					   y = (L1*y1+L2*y2)/L_tot
   z = (L1*z1+L2*z2)/L_tot					   z = (L1*z1+L2*z2)/L_tot


   r = sqrt(x**2+y**2)						   r = sqrt(x**2+y**2)
   b = atan(z/r)*180/pi						   b = atan(z/r)*180/pi
   l = atan(y/x)*180/pi						   l = atan(y/x)*180/pi
   								   

   if sign(x) < 0 : l+=180					   if sign(x) < 0 : l+=180
   if sign(x) > 0 and sign(y) < 0: l+=360   			   if sign(x) > 0 and sign(y) < 0: l+=360   
   								   
   if d1==0 or d2==0:						   if d1==0 or d2==0:
      d = 0.							      d = 0.
   else:							   else:
      d = sqrt(x**2+y**2+z**2)					      d = sqrt(x**2+y**2+z**2)


   								   
   return l, b, d						   return l, b, d

# **************************************			# **************************************

def galJoint(galNode1, galNode2, ID):				def galJoint(galNode1, galNode2, ID):
   								   
   if galNode1==None and galNode2!=None:			   if galNode1==None and galNode2!=None:
     return galNode2						     return galNode2
   if galNode2==None and galNode1!=None:			   if galNode2==None and galNode1!=None:
     return galNode1						     return galNode1
   if galNode2==None and galNode1==None:			   if galNode2==None and galNode1==None:
     return None						     return None
   								   

   								   
   L1 = 10**galNode1.logK					   L1 = 10**galNode1.logK
   L2 = 10**galNode2.logK					   L2 = 10**galNode2.logK
   								   
   if L1 == 1: L1 = 0						   if L1 == 1: L1 = 0
   if L2 == 1: L2 = 0						   if L2 == 1: L2 = 0
   								   
   L_tot = L1 + L2						   L_tot = L1 + L2
   								   
   logK_tot = log10(L_tot)					   logK_tot = log10(L_tot)
   								   
   								   
   d1 = galNode1.dcf2						   d1 = galNode1.dcf2
   d2 = galNode2.dcf2						   d2 = galNode2.dcf2
   								   
   sgl1 = galNode1.sgl						   sgl1 = galNode1.sgl
   sgl2 = galNode2.sgl						   sgl2 = galNode2.sgl
   sgb1 = galNode1.sgb						   sgb1 = galNode1.sgb
   sgb2 = galNode2.sgb						   sgb2 = galNode2.sgb
   								   
   gl1 = galNode1.gl						   gl1 = galNode1.gl
   gl2 = galNode2.gl						   gl2 = galNode2.gl
   gb1 = galNode1.gb						   gb1 = galNode1.gb
   gb2 = galNode2.gb 						   gb2 = galNode2.gb 
   								   
   ra1 = galNode1.ra						   ra1 = galNode1.ra
   ra2 = galNode2.ra						   ra2 = galNode2.ra
   dec1 = galNode1.dec						   dec1 = galNode1.dec
   dec2 = galNode2.dec   					   dec2 = galNode2.dec   
      								      
   								   
   								   
   								   
   v1 = galNode1.Vls						   v1 = galNode1.Vls
   v2 = galNode2.Vls						   v2 = galNode2.Vls
   								   
   if d1==0 and d2==0:						   if d1==0 and d2==0:
     d1 = 0.							     d1 = 0.
     d2 = 0							     d2 = 0
     d_final=0.							     d_final=0.
   elif d1!=0 and d2==0:					   elif d1!=0 and d2==0:
     d2 = galNode2.Vls/H0					     d2 = galNode2.Vls/H0
     d_final=0							     d_final=0
   elif d2!=0 and d1==0:					   elif d2!=0 and d1==0:
     d1 = galNode1.Vls/H0					     d1 = galNode1.Vls/H0
     d_final=0							     d_final=0
   else:							   else:
     d_final=1							     d_final=1
     								     
     								     
   sgl, sgb, d   = barycenter(L1,  sgl1,  sgb1, d1, L2,  sgl2	   sgl, sgb, d   = barycenter(L1,  sgl1,  sgb1, d1, L2,  sgl2
   gl, gb, d = barycenter(L1, gl1, gb1, d1, L2, gl2, gb2, d2)	   gl, gb, d = barycenter(L1, gl1, gb1, d1, L2, gl2, gb2, d2)
   ra, dec, d = barycenter(L1, ra1, dec1, d1, L2, ra2, dec2, 	   ra, dec, d = barycenter(L1, ra1, dec1, d1, L2, ra2, dec2, 
   								   
   if d_final==0:						   if d_final==0:
      d = 0							      d = 0
      dm = 0							      dm = 0
   else:							   else:
      dm = 5*log(d)+25						      dm = 5*log(d)+25


   n1 = galNode1.subGalaxies					   n1 = galNode1.subGalaxies
   n2 = galNode2.subGalaxies   					   n2 = galNode2.subGalaxies   
   Vls = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)	   Vls = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)
   if galNode1.v_av == 0 and galNode2.v_av != 0 :		   if galNode1.v_av == 0 and galNode2.v_av != 0 :
     Vls = galNode2.v_av					     Vls = galNode2.v_av
     if L1 == 0:						     if L1 == 0:
       galNode1.logK = m_logK(galNode1.Ks, Vls, d=galNode1.dc	       galNode1.logK = m_logK(galNode1.Ks, Vls, d=galNode1.dc
       								       
   if galNode1.v_av != 0 and galNode2.v_av == 0 :		   if galNode1.v_av != 0 and galNode2.v_av == 0 :
     Vls = galNode1.v_av   					     Vls = galNode1.v_av   
     if L2 == 0:						     if L2 == 0:
       galNode2.logK = m_logK(galNode2.Ks, Vls, d=galNode2.dc	       galNode2.logK = m_logK(galNode2.Ks, Vls, d=galNode2.dc

   newNode = GalxyNode(ID, gl, gb, sgl, sgb, Vls, 0, 0, d, 0,	   newNode = GalxyNode(ID, gl, gb, sgl, sgb, Vls, 0, 0, d, 0,
   newNode.logK = logK_tot					   newNode.logK = logK_tot
   newNode.ra = ra						   newNode.ra = ra
   newNode.dec = dec						   newNode.dec = dec
   newNode.coordinate_src = 'GRPhead'				   newNode.coordinate_src = 'GRPhead'
   newNode.Ty_src  =  'GRPhead'					   newNode.Ty_src  =  'GRPhead'
   newNode.Ks_src  =  'GRPhead'					   newNode.Ks_src  =  'GRPhead'
   newNode.Vls_src =  'GRPhead'					   newNode.Vls_src =  'GRPhead'
   newNode.objname =  'GRPhead'					   newNode.objname =  'GRPhead'
   								   
   newNode.subGalaxies = n1 + n2				   newNode.subGalaxies = n1 + n2
   newNode.level = max([galNode1.level, galNode2.level]) + 1	   newNode.level = max([galNode1.level, galNode2.level]) + 1
   								   
   galNode1.topLevel = newNode.level				   galNode1.topLevel = newNode.level
   galNode2.topLevel = newNode.level				   galNode2.topLevel = newNode.level
   								   
   newNode.sumDist = galNode1.sumDist + galNode2.sumDist	   newNode.sumDist = galNode1.sumDist + galNode2.sumDist
   newNode.sumError = galNode1.sumError + galNode2.sumError	   newNode.sumError = galNode1.sumError + galNode2.sumError
   								   
   if newNode.sumDist !=0 and newNode.sumError!=0:		   if newNode.sumDist !=0 and newNode.sumError!=0:
        meanDist = newNode.sumDist/newNode.sumError		        meanDist = newNode.sumDist/newNode.sumError
        meanDistErr = sqrt(1/newNode.sumError)/newNode.sumDis	        meanDistErr = sqrt(1/newNode.sumError)/newNode.sumDis
        mDM  = (log10(meanDist)+5.)*5.				        mDM  = (log10(meanDist)+5.)*5.
        meDM = meanDistErr / (0.2*log(10.))			        meDM = meanDistErr / (0.2*log(10.))
   else:							   else:
        mDM = 0.						        mDM = 0.
        meDM = 0.						        meDM = 0.
        meanDist = 0.						        meanDist = 0.
        meanDistErr = 0.					        meanDistErr = 0.
    								    
   								   
   								   
   if galNode1.mDist != 0 and galNode2.mDist != 0:		   if galNode1.mDist != 0 and galNode2.mDist != 0:
     if L1 > L2: 						     if L1 > L2: 
       newNode.mDist = galNode1.mDist				       newNode.mDist = galNode1.mDist
       newNode.mDistErr = galNode1.mDistErr			       newNode.mDistErr = galNode1.mDistErr
       newNode.mDM = galNode1.mDM				       newNode.mDM = galNode1.mDM
       newNode.meDM = galNode1.meDM				       newNode.meDM = galNode1.meDM
     else:							     else:
       newNode.mDist = galNode2.mDist				       newNode.mDist = galNode2.mDist
       newNode.mDistErr = galNode2.mDistErr			       newNode.mDistErr = galNode2.mDistErr
       newNode.mDM = galNode2.mDM				       newNode.mDM = galNode2.mDM
       newNode.meDM = galNode2.meDM				       newNode.meDM = galNode2.meDM
       								       
   elif galNode1.mDist != 0:					   elif galNode1.mDist != 0:
     newNode.mDist = galNode1.mDist				     newNode.mDist = galNode1.mDist
     newNode.mDistErr = galNode1.mDistErr			     newNode.mDistErr = galNode1.mDistErr
     newNode.mDM = galNode1.mDM					     newNode.mDM = galNode1.mDM
     newNode.meDM = galNode1.meDM				     newNode.meDM = galNode1.meDM
     								     
   elif galNode2.mDist != 0:					   elif galNode2.mDist != 0:
     newNode.mDist = galNode2.mDist				     newNode.mDist = galNode2.mDist
     newNode.mDistErr = galNode2.mDistErr			     newNode.mDistErr = galNode2.mDistErr
     newNode.mDM = galNode2.mDM					     newNode.mDM = galNode2.mDM
     newNode.meDM = galNode2.meDM				     newNode.meDM = galNode2.meDM
   else:							   else:
     newNode.mDist = meanDist					     newNode.mDist = meanDist
     newNode.mDistErr = meanDistErr				     newNode.mDistErr = meanDistErr
     newNode.mDM = mDM						     newNode.mDM = mDM
     newNode.meDM = meDM					     newNode.meDM = meDM

   								   
   # Brighter galaxy is the left child				   # Brighter galaxy is the left child
   if L1 >= L2:							   if L1 >= L2:
      newNode.left = galNode1					      newNode.left = galNode1
      newNode.right = galNode2					      newNode.right = galNode2
      newNode.nest = galNode1.nest				      newNode.nest = galNode1.nest
   else:							   else:
      newNode.left = galNode2					      newNode.left = galNode2
      newNode.right = galNode1					      newNode.right = galNode1
      newNode.nest = galNode2.nest				      newNode.nest = galNode2.nest
   								   

   newNode.v_av = Vls						   newNode.v_av = Vls
   newNode.v2_av = (n1*galNode1.v2_av + n2*galNode2.v2_av) / 	   newNode.v2_av = (n1*galNode1.v2_av + n2*galNode2.v2_av) / 
   if galNode1.v2_av == 0 and galNode2.v2_av != 0 :		   if galNode1.v2_av == 0 and galNode2.v2_av != 0 :
     newNode.v2_av = galNode2.v2_av				     newNode.v2_av = galNode2.v2_av
   if galNode1.v2_av != 0 and galNode2.v2_av == 0 :		   if galNode1.v2_av != 0 and galNode2.v2_av == 0 :
     newNode.v2_av = galNode1.v2_av   				     newNode.v2_av = galNode1.v2_av   
   								   
   								   
   								   
   if (newNode.v2_av - newNode.v_av**2) > 0:			   if (newNode.v2_av - newNode.v_av**2) > 0:
      newNode.sigma =  sqrt(newNode.v2_av - newNode.v_av**2) 	      newNode.sigma =  sqrt(newNode.v2_av - newNode.v_av**2) 
   								   
   								   
   newNode.R_theta = Theta_max(newNode)				   newNode.R_theta = Theta_max(newNode)
   								   

   if newNode.sigma == 0:					   if newNode.sigma == 0:
     sig = 1							     sig = 1
   else:							   else:
     sig = newNode.sigma					     sig = newNode.sigma
   								   

   mass = Mass(L_tot)						   mass = Mass(L_tot)
   newNode.M_v2 = mass						   newNode.M_v2 = mass
   newNode.R_2t2 = 0.215*((mass/1.E12)**(1./3))  # Mpc		   newNode.R_2t2 = 0.215*((mass/1.E12)**(1./3))  # Mpc

   return newNode						   return newNode
# **************************************			# **************************************
# The calls definition of a galaxy				# The calls definition of a galaxy
# each node contains all esseential property of a galaxy	# each node contains all esseential property of a galaxy
# when gaalxies get connected along a tree, the new nodes	# when gaalxies get connected along a tree, the new nodes
# are defnided as new entities, which are new galaxies in	# are defnided as new entities, which are new galaxies in
# the context of this code					# the context of this code
class GalxyNode:						class GalxyNode:
  								  
  # field variables						  # field variables
  id = 0.							  id = 0.
  l = 0.							  l = 0.
  b = 0.							  b = 0.
  gl = 0.							  gl = 0.
  gb = 0.							  gb = 0.
  Vls = 0.							  Vls = 0.
  logK = 0.							  logK = 0.
  subGalaxies = 1						  subGalaxies = 1
  level = 0							  level = 0
  R_theta = 0.							  R_theta = 0.
  nest = 0							  nest = 0
  Ty = 0.							  Ty = 0.
  v_av = 0.							  v_av = 0.
  v2_av = 0.							  v2_av = 0.
  sigma = 0.  # velocity dispersion (if a node contains sever	  sigma = 0.  # velocity dispersion (if a node contains sever
  dcf2 = 0.							  dcf2 = 0.
  ed = 0.							  ed = 0.
  dcf2Copy = 0.							  dcf2Copy = 0.
  edCopy = 0.							  edCopy = 0.
  mDist = 0.							  mDist = 0.
  mDistErr = 0.							  mDistErr = 0.
  Ks = -100000.0						  Ks = -100000.0
  inGroup = 0.							  inGroup = 0.
  								  
  M_v2 = 0.							  M_v2 = 0.
  R_2t2 = 0.							  R_2t2 = 0.
  								  
  sumDist = 0.							  sumDist = 0.
  sumError = 0.							  sumError = 0.
  								  
  								  
  ra  =  0.							  ra  =  0.
  dec = 0.							  dec = 0.
  coordinate_src = 'NotSet'					  coordinate_src = 'NotSet'
  Ty_src  = 'NotSet'						  Ty_src  = 'NotSet'
  B_mag   = -100000.0						  B_mag   = -100000.0
  Ks_src  = 'NotSet'						  Ks_src  = 'NotSet'
  Vls_src = 'NotSet'						  Vls_src = 'NotSet'
  objname = 'NotSet'						  objname = 'NotSet'
  								  
  flag = 0							  flag = 0
  								  
  left = None							  left = None
  right = None							  right = None
  								  
  # Class constructor						  # Class constructor
  def __init__(self, id, gl, gb, sgl, sgb, Vls, Ks, Ty, dcf2,	  def __init__(self, id, gl, gb, sgl, sgb, Vls, Ks, Ty, dcf2,
    								    
    								    
    #dcf2 = 0							    #dcf2 = 0
    self.id = id						    self.id = id
    self.gl = gl						    self.gl = gl
    self.gb = gb    						    self.gb = gb    
    self.sgl = sgl						    self.sgl = sgl
    self.sgb = sgb						    self.sgb = sgb
    self.Vls = Vls						    self.Vls = Vls
    self.Vhelio = 0.						    self.Vhelio = 0.
    self.Ks = Ks						    self.Ks = Ks
    self.Ty = Ty						    self.Ty = Ty
    								    
    								    
    								    
    								    
    								    
    self.dcf2 = dcf2						    self.dcf2 = dcf2
    self.ed = ed						    self.ed = ed
    self.dm = dm						    self.dm = dm
    self.edm = edm						    self.edm = edm
    								    
    								    
    								    
    								    
    self.dcf2Copy = dcf2					    self.dcf2Copy = dcf2
    self.edCopy = ed    					    self.edCopy = ed    
    self.dmCopy = dm						    self.dmCopy = dm
    self.edmCopy = edm        					    self.edmCopy = edm        
    								    
    self.mDist = dcf2						    self.mDist = dcf2
    self.mDistErr = ed						    self.mDistErr = ed
    self.mDM = dm						    self.mDM = dm
    self.meDM = edm    						    self.meDM = edm    
    								    
    self.left = None						    self.left = None
    self.right = None						    self.right = None
    								    
    self.subGalaxies = 1					    self.subGalaxies = 1
    self.level = 0						    self.level = 0
    self.Rgroup = 0.  						    self.Rgroup = 0.  
    self.R_theta = 0.						    self.R_theta = 0.
    self.nest = id						    self.nest = id
    self.sigma = 0. 						    self.sigma = 0. 
    self.v_av = Vls						    self.v_av = Vls
    self.v2_av = Vls**2						    self.v2_av = Vls**2
    								    
    # 0: if the galaxy is NOT in a group			    # 0: if the galaxy is NOT in a group
    # 1: if the galaxy falls in a group				    # 1: if the galaxy falls in a group
    self.inGroup = 0.						    self.inGroup = 0.
    								    
    self.logK = m_logK(Ks, Vls, d=dcf2)			      |	    self.logK = m_logK(Ks, Vls)
 								 
    self.M_v2 = Mass(10**self.logK)				    self.M_v2 = Mass(10**self.logK)
    self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc	    self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc


    								    
    if dcf2!=0 and ed!=0:					    if dcf2!=0 and ed!=0:
        							        
        err = dcf2*ed						        err = dcf2*ed
        self.sumDist = 1. * dcf2 / (err**2)			        self.sumDist = 1. * dcf2 / (err**2)
        self.sumError = 1. / (err**2)				        self.sumError = 1. / (err**2)
    else:							    else:
        self.sumDist = 0.					        self.sumDist = 0.
        self.sumError = 0.					        self.sumError = 0.
        							        
    								    
    								    
  def setMeanDist(self, Dist, errDist, GRP_vel = 0, dist=0):	  def setMeanDist(self, Dist, errDist, GRP_vel = 0, dist=0):
    								    
      self.mDist = Dist						      self.mDist = Dist
      self.mDistErr = errDist					      self.mDistErr = errDist
      if GRP_vel == 0 : 					      if GRP_vel == 0 : 
	vel = self.Vls							vel = self.Vls
      else: vel = GRP_vel					      else: vel = GRP_vel
      self.logK = m_logK(self.Ks, vel, d=dist)			      self.logK = m_logK(self.Ks, vel, d=dist)
      self.M_v2 = Mass(10**self.logK)				      self.M_v2 = Mass(10**self.logK)
      self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc	      self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
      								      
    								    


# **************************************			# **************************************
# Given the apparent magnitude of a galaxy, m			# Given the apparent magnitude of a galaxy, m
# it returns the absolute luminosity 				# it returns the absolute luminosity 
# Vls is the radial velocity of the galaxy			# Vls is the radial velocity of the galaxy
# So, the distance to a galaxy is estimated 			# So, the distance to a galaxy is estimated 
# based on its velocity						# based on its velocity
def m_logK(m, Vls, d=0):					def m_logK(m, Vls, d=0):
    								    
    Mk_sun = 3.28   # K-band					    Mk_sun = 3.28   # K-band
    distance = Vls / H0  # Mpc					    distance = Vls / H0  # Mpc
    								    
    if distance < 1:						    if distance < 1:
      distance = 1 # Mpc					      distance = 1 # Mpc
    								    
    if d!=0:							    if d!=0:
      distance = d						      distance = d
    								    
    M = m - 5*log10(distance) - 30 + 5				    M = m - 5*log10(distance) - 30 + 5
    logK = -0.4 * (M - Mk_sun)					    logK = -0.4 * (M - Mk_sun)
    								    
    if Vls == 0 and d==0:				      |	    if Vls == 0:
      logK = 0							      logK = 0
    								    
    								    
    return logK							    return logK
# **************************************			# **************************************
# It returns all leaf-nodes at the bottom of a gaalxy tree.	# It returns all leaf-nodes at the bottom of a gaalxy tree.
# root: is the root of the tree					# root: is the root of the tree
# It works recursively, so for large trees, it might 		# It works recursively, so for large trees, it might 
# exhaust the recursion stack, and the codes crashe		# exhaust the recursion stack, and the codes crashe
# so some exception handling is required in this case		# so some exception handling is required in this case
def NodeLeaf(root):						def NodeLeaf(root):
  list = [root]							  list = [root]
  NodeLeafCore(root, list)					  NodeLeafCore(root, list)
  return list							  return list

def NodeLeafCore(root, list):  					def NodeLeafCore(root, list):  
  if root.left == None:						  if root.left == None:

    list.append(root)						    list.append(root)
  else:								  else:
    NodeLeafCore(root.left, list)				    NodeLeafCore(root.left, list)
    NodeLeafCore(root.right, list)				    NodeLeafCore(root.right, list)
    								    
    								    
   								   


#############################################################	#############################################################
def readgalList(table, skyPatch='none', sgl_min=0, sgl_max=0,	def readgalList(table, skyPatch='none', sgl_min=0, sgl_max=0,

  id   = table['pgc']						  id   = table['pgc']
  sgl  = table['sgl']						  sgl  = table['sgl']
  sgb  = table['sgb']						  sgb  = table['sgb']
  Vls  = table['Vls']						  Vls  = table['Vls']
  Vhelio = table['Vhelio']					  Vhelio = table['Vhelio']
  Ks   = table['Ks'] 						  Ks   = table['Ks'] 
  dcf2 = table['dcf2'] 						  dcf2 = table['dcf2'] 
  ed   = table['ed'] 						  ed   = table['ed'] 
  dm = table['DM'] 						  dm = table['DM'] 
  edm   = table['eDM'] 						  edm   = table['eDM'] 
  								  
  Ty   = table['Ty'] 						  Ty   = table['Ty'] 
  								  
  gl  = table['gl']						  gl  = table['gl']
  gb  = table['gb']  						  gb  = table['gb']  
  								  
  								  
  ra  =  table['ra']						  ra  =  table['ra']
  dec = table['dec']						  dec = table['dec']
  coordinate_src = table['coordinate_src']			  coordinate_src = table['coordinate_src']
  Ty_src  = table['Ty_src']					  Ty_src  = table['Ty_src']
  B_mag   = table['B_mag']					  B_mag   = table['B_mag']
  Ks_src  = table['Ks_src']					  Ks_src  = table['Ks_src']
  Vls_src = table['Vls_src']					  Vls_src = table['Vls_src']
  objname = table['objname']					  objname = table['objname']
   								   
  N_galaxies = len(id)						  N_galaxies = len(id)


  print "Reading the data file .... "				  print "Reading the data file .... "
  galList = []							  galList = []
  								  

    								    

  i = 0								  i = 0
  								  
  for i in range(N_galaxies):					  for i in range(N_galaxies):
    								    
    								    
    								    
    								    
    ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)		    ang12 = (180./pi)*angle(sgl[i], sgb[i], sglV, sgbV)
    								    
    								    
    Bool = False						    Bool = False
    								    
    								    
    inHydra = False						    inHydra = False
    if (180./pi)*angle(sgl[i], sgb[i], 139.3599, -37.7051) < 	    if (180./pi)*angle(sgl[i], sgb[i], 139.3599, -37.7051) < 

    								    
    								    
    								    
    inCentaurus = False						    inCentaurus = False
    if (180./pi)*angle(sgl[i], sgb[i], 156.2981, -11.6740) < 	    if (180./pi)*angle(sgl[i], sgb[i], 156.2981, -11.6740) < 

    if  skyPatch=='manual':					    if  skyPatch=='manual':
      if (Vls[i] < Vls_max or inHydra or inCentaurus) and sgl	      if (Vls[i] < Vls_max or inHydra or inCentaurus) and sgl
        Bool = True						        Bool = True
    elif skyPatch=='local':					    elif skyPatch=='local':
      if (Vls[i] < Vls_max or inHydra or inCentaurus):		      if (Vls[i] < Vls_max or inHydra or inCentaurus):
	Bool = True							Bool = True
    else:							    else:
      if (Vls[i] < Vls_max or inHydra or inCentaurus): 		      if (Vls[i] < Vls_max or inHydra or inCentaurus): 
	  Bool = True							  Bool = True
	  								  
    if id[i] == 31599:						    if id[i] == 31599:
      print id[i], (180./pi)*angle(sgl[i], sgb[i], 139.3599, 	      print id[i], (180./pi)*angle(sgl[i], sgb[i], 139.3599, 
      								      
    if Bool:         						    if Bool:         

	   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[		   node = GalxyNode(id[i], gl[i], gb[i], sgl[i], sgb[
	   node.ra = ra[i]						   node.ra = ra[i]
	   node.dec = dec[i]						   node.dec = dec[i]
	   node.coordinate_src = coordinate_src[i]			   node.coordinate_src = coordinate_src[i]
	   node.Ty_src = Ty_src[i]					   node.Ty_src = Ty_src[i]
	   node.B_mag =  B_mag[i]					   node.B_mag =  B_mag[i]
	   node.Ks_src = Ks_src[i]					   node.Ks_src = Ks_src[i]
	   node.Vls_src = Vls_src[i]					   node.Vls_src = Vls_src[i]
	   node.objname = objname[i]					   node.objname = objname[i]
	   node.Vhelio = Vhelio[i]					   node.Vhelio = Vhelio[i]
           galList.append(node)					           galList.append(node)
           							           
           if id[i] in ignore_list or Ks[i] < 0 or Vls[i]==0:	           if id[i] in ignore_list or Ks[i] < 0 or Vls[i]==0:
	     #print "dis_include:", id[i]				     #print "dis_include:", id[i]
	     node.inGroup = -1						     node.inGroup = -1
           							           
  print "No. of galaxies: ", len(galList)			  print "No. of galaxies: ", len(galList)
  print "Data file loaded .... "				  print "Data file loaded .... "
  								  

  								  
  								  
  return galList						  return galList
#############################################################	#############################################################
#############################################################	#############################################################
# b is galactic latitude [deg]					# b is galactic latitude [deg]
# 2*theta + sin(2*theta) = Pi * sin(b)				# 2*theta + sin(2*theta) = Pi * sin(b)

def theta(b):							def theta(b):
  								  
  if b == 90.: return pi/2.					  if b == 90.: return pi/2.
  if b == -90.: return -pi/2.					  if b == -90.: return -pi/2.
  								  
  b = b*pi/180.							  b = b*pi/180.
  theta0 = b							  theta0 = b
  theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b))/	  theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b))/
  								  
  while abs(theta1-theta0) > 0.01:				  while abs(theta1-theta0) > 0.01:
    theta0 = theta1						    theta0 = theta1
    theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b)	    theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b)
  								  
  return theta1							  return theta1
  								  
#############################################################	#############################################################

def xymap(l,b,l0):						def xymap(l,b,l0):
  								  
  t = theta(b)							  t = theta(b)
  x = 180.-((2*90.)/180.)*(l0-l)*cos(t)				  x = 180.-((2*90.)/180.)*(l0-l)*cos(t)
  y = 90.*sin(t)						  y = 90.*sin(t)
  								  
  return x, y 							  return x, y 
#############################################################	#############################################################

def plot_border(l0):						def plot_border(l0):
  								  
  								  
  X = []							  X = []
  Y = []							  Y = []

  								  
  l = 0.							  l = 0.
  for b in np.arange(-90,90,0.01):				  for b in np.arange(-90,90,0.01):
    x, y = xymap(l,b,l0)					    x, y = xymap(l,b,l0)
    X.append(x)							    X.append(x)
    Y.append(y)							    Y.append(y)
  l = 360							  l = 360
  for b in np.arange(90,-90,-0.01):				  for b in np.arange(90,-90,-0.01):
    x, y = xymap(l,b,l0)					    x, y = xymap(l,b,l0)
    X.append(x)							    X.append(x)
    Y.append(y)							    Y.append(y)
  								  
  X.append(X[0])						  X.append(X[0])
  Y.append(Y[0])						  Y.append(Y[0])

  plt.plot(X, Y, '-', markersize = 1, linewidth=2., color='#0	  plt.plot(X, Y, '-', markersize = 1, linewidth=2., color='#0
  								  
#############################################################	#############################################################

def plot_galaxies(inFile, l0):					def plot_galaxies(inFile, l0):
  								  
  table = np.genfromtxt( inFile , delimiter=',', filling_valu	  table = np.genfromtxt( inFile , delimiter=',', filling_valu
  								  
  #print table.dtype						  #print table.dtype
  id   = table['pgc']						  id   = table['pgc']
  gl  = table['sgl']						  gl  = table['sgl']
  gb  = table['sgb']						  gb  = table['sgb']
  N_galaxies = len(id)						  N_galaxies = len(id)
  								  
  X0 = []							  X0 = []
  Y0 = []							  Y0 = []
  for i in range(0, N_galaxies):				  for i in range(0, N_galaxies):
    x, y = xymap(gl[i],gb[i],l0)				    x, y = xymap(gl[i],gb[i],l0)
    X0.append(x)						    X0.append(x)
    Y0.append(y)						    Y0.append(y)
  plt.plot(X0, Y0, '.', markersize = 1, color='#696969')  # g	  plt.plot(X0, Y0, '.', markersize = 1, color='#696969')  # g
  								  
  								  
#############################################################	#############################################################
def groupPlot(northGList, southGList):				def groupPlot(northGList, southGList):
  								  
  l0 = 180							  l0 = 180
  								  
  fig = plt.figure(figsize=(12, 12*0.5), dpi=100)		  fig = plt.figure(figsize=(12, 12*0.5), dpi=100)
  ax = fig.add_axes([0.13, 0.13, 0.80,  0.80]) 			  ax = fig.add_axes([0.13, 0.13, 0.80,  0.80]) 
  plt.ylim(-100,100)						  plt.ylim(-100,100)
  plt.xlim(380,-20)						  plt.xlim(380,-20)
  plt.xlabel("SGL (deg)", fontsize=20)				  plt.xlabel("SGL (deg)", fontsize=20)
  plt.ylabel("SGB (deg)", fontsize=20)				  plt.ylabel("SGB (deg)", fontsize=20)
  plt.yticks(fontsize=16)					  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)					  plt.xticks(fontsize=16)
  								  
  plot_border(l0)						  plot_border(l0)
  plot_galaxies('AllSky.north.csv', l0)				  plot_galaxies('AllSky.north.csv', l0)
  plot_galaxies('AllSky.south.csv', l0)				  plot_galaxies('AllSky.south.csv', l0)
  								  
  # Virgo Border (6.8 deg)					  # Virgo Border (6.8 deg)
  u = np.arange(0,1,0.010)					  u = np.arange(0,1,0.010)
  X = u*0.							  X = u*0.
  Y = u*0.							  Y = u*0.
  for i in range(0,len(u)):					  for i in range(0,len(u)):
     x = 6.8*np.cos(u[i]*2*pi) + sglV				     x = 6.8*np.cos(u[i]*2*pi) + sglV
     y = 6.8*np.sin(u[i]*2*pi) + sgbV				     y = 6.8*np.sin(u[i]*2*pi) + sgbV
     X[i], Y[i] = xymap(x,y, l0)				     X[i], Y[i] = xymap(x,y, l0)
  								  
  line00, = plt.plot(X,Y, 'r.', markersize = 2)   		  line00, = plt.plot(X,Y, 'r.', markersize = 2)   
  line00.set_dashes([2, 2]) 					  line00.set_dashes([2, 2]) 
  								  
  groupPlotak(northGList, l0)					  groupPlotak(northGList, l0)
  groupPlotak(southGList, l0)					  groupPlotak(southGList, l0)
##########							##########

def groupPlotak(G_list, l0):					def groupPlotak(G_list, l0):
  								  
 								 
  NoGroups = len(G_list)					  NoGroups = len(G_list)
  print "Number of groups: ", NoGroups				  print "Number of groups: ", NoGroups
  								  
  col = 0							  col = 0

  for i in range(0, NoGroups):  # for all groups		  for i in range(0, NoGroups):  # for all groups

    								    
    Key = True							    Key = True
    if Key:							    if Key:
	random.seed(G_list[i][0].nest )					random.seed(G_list[i][0].nest )
	Red, Blue, Green = random.random(), random.random(), 		Red, Blue, Green = random.random(), random.random(), 
									
									
	r = G_list[i][0].R_theta					r = G_list[i][0].R_theta
									
	d_theta = 0.001							d_theta = 0.001
	u = np.arange(0,2*pi,d_theta)					u = np.arange(0,2*pi,d_theta)
	X = u*0.							X = u*0.
        Y = u*0.						        Y = u*0.
        for q in range(0,len(u)):				        for q in range(0,len(u)):
           x = r*np.cos(u[q]) + G_list[i][0].sgl		           x = r*np.cos(u[q]) + G_list[i][0].sgl
           y = r*np.sin(u[q]) + G_list[i][0].sgb		           y = r*np.sin(u[q]) + G_list[i][0].sgb
           X[q], Y[q] = xymap(x,y, l0)				           X[q], Y[q] = xymap(x,y, l0)

									
	if r <= 6.3:							if r <= 6.3:
	  line, = plt.plot(X,Y, '-', markersize = 2, color=(R		  line, = plt.plot(X,Y, '-', markersize = 2, color=(R
	  line.set_dashes([8, 3]) 					  line.set_dashes([8, 3]) 
									

	X  = []								X  = []
	Y = []								Y = []
									

	for j in range(1, len(G_list[i])):				for j in range(1, len(G_list[i])):

	    x, y = xymap(G_list[i][j].sgl,G_list[i][j].sgb, l		    x, y = xymap(G_list[i][j].sgl,G_list[i][j].sgb, l
	    X.append(x)							    X.append(x)
	    Y.append(y)							    Y.append(y)
	    								    
	plt.plot(X, Y, 'o', markersize = 3, color=(Red, Blue,		plt.plot(X, Y, 'o', markersize = 3, color=(Red, Blue,
        							        

#############################################################	#############################################################
#############################################################	#############################################################
def groupWrite(outfile, G_list, galList):			def groupWrite(outfile, G_list, galList):
  								  
  								  
  # finding the leave nodes					  # finding the leave nodes
  #galList = NodeLeaf(root)					  #galList = NodeLeaf(root)
  #galList.pop(0)						  #galList.pop(0)
  								  
  #for grp in G_list:						  #for grp in G_list:
    #for gal in grp:						    #for gal in grp:
      #if gal.id == 31599:					      #if gal.id == 31599:
        #print 'I have in my catalog'				        #print 'I have in my catalog'
      								      
      								      
  NoGroups = len(G_list)					  NoGroups = len(G_list)
  #print "Number of groups: ", NoGroups  			  #print "Number of groups: ", NoGroups  

  myTable = Table()						  myTable = Table()
    								    
    								    
  empty = []							  empty = []
  myTable.add_column(Column(data=empty,name='pgc', dtype=np.d	  myTable.add_column(Column(data=empty,name='pgc', dtype=np.d
  myTable.add_column(Column(data=empty,name='flag', dtype=np.	  myTable.add_column(Column(data=empty,name='flag', dtype=np.
  myTable.add_column(Column(data=empty,name='ra', format='%0.	  myTable.add_column(Column(data=empty,name='ra', format='%0.
  myTable.add_column(Column(data=empty,name='dec', format='%0	  myTable.add_column(Column(data=empty,name='dec', format='%0
  myTable.add_column(Column(data=empty,name='gl', format='%0.	  myTable.add_column(Column(data=empty,name='gl', format='%0.
  myTable.add_column(Column(data=empty,name='gb', format='%0.	  myTable.add_column(Column(data=empty,name='gb', format='%0.
  myTable.add_column(Column(data=empty,name='sgl', format='%0	  myTable.add_column(Column(data=empty,name='sgl', format='%0
  myTable.add_column(Column(data=empty,name='sgb', format='%0	  myTable.add_column(Column(data=empty,name='sgb', format='%0
  myTable.add_column(Column(data=empty,name='Ty'))		  myTable.add_column(Column(data=empty,name='Ty'))
  myTable.add_column(Column(data=empty,name='B_mag', format='	  myTable.add_column(Column(data=empty,name='B_mag', format='
  myTable.add_column(Column(data=empty,name='Ks', format='%0.	  myTable.add_column(Column(data=empty,name='Ks', format='%0.
  myTable.add_column(Column(data=empty,name='logK', format='%	  myTable.add_column(Column(data=empty,name='logK', format='%
  myTable.add_column(Column(data=empty,name='Vls', format='%0	  myTable.add_column(Column(data=empty,name='Vls', format='%0
  myTable.add_column(Column(data=empty,name='Vhelio', format=	  myTable.add_column(Column(data=empty,name='Vhelio', format=
  myTable.add_column(Column(data=empty,name='dcf2', format='%	  myTable.add_column(Column(data=empty,name='dcf2', format='%
  myTable.add_column(Column(data=empty,name='ed', format='%0.	  myTable.add_column(Column(data=empty,name='ed', format='%0.

  myTable.add_column(Column(data=empty,name='mDist', format='	  myTable.add_column(Column(data=empty,name='mDist', format='
  myTable.add_column(Column(data=empty,name='mDistErr', forma	  myTable.add_column(Column(data=empty,name='mDistErr', forma
  myTable.add_column(Column(data=empty,name='R_theta', format	  myTable.add_column(Column(data=empty,name='R_theta', format
  myTable.add_column(Column(data=empty,name='sigmaP_dyn', for	  myTable.add_column(Column(data=empty,name='sigmaP_dyn', for
  myTable.add_column(Column(data=empty,name='sigmaP_lum', for	  myTable.add_column(Column(data=empty,name='sigmaP_lum', for
  								  
  myTable.add_column(Column(data=empty,name='Mv_dyn', format=	  myTable.add_column(Column(data=empty,name='Mv_dyn', format=
  myTable.add_column(Column(data=empty,name='Mv_lum', format=	  myTable.add_column(Column(data=empty,name='Mv_lum', format=
  myTable.add_column(Column(data=empty,name='Rg_angular', for	  myTable.add_column(Column(data=empty,name='Rg_angular', for
  myTable.add_column(Column(data=empty,name='Rg_dyn', format=	  myTable.add_column(Column(data=empty,name='Rg_dyn', format=
  myTable.add_column(Column(data=empty,name='R2t_dyn', format	  myTable.add_column(Column(data=empty,name='R2t_dyn', format
  myTable.add_column(Column(data=empty,name='R2t_lum', format	  myTable.add_column(Column(data=empty,name='R2t_lum', format
  myTable.add_column(Column(data=empty,name='tX_dyn', format=	  myTable.add_column(Column(data=empty,name='tX_dyn', format=
  myTable.add_column(Column(data=empty,name='tX_lum', format=	  myTable.add_column(Column(data=empty,name='tX_lum', format=
  myTable.add_column(Column(data=empty,name='No_Galaxies',dty	  myTable.add_column(Column(data=empty,name='No_Galaxies',dty
  myTable.add_column(Column(data=empty,name='nest', dtype=np.	  myTable.add_column(Column(data=empty,name='nest', dtype=np.

  myTable.add_column(Column(data=empty,name='coordinate_src',	  myTable.add_column(Column(data=empty,name='coordinate_src',
  myTable.add_column(Column(data=empty,name='Ty_src', dtype='	  myTable.add_column(Column(data=empty,name='Ty_src', dtype='
  myTable.add_column(Column(data=empty,name='Ks_src', dtype='	  myTable.add_column(Column(data=empty,name='Ks_src', dtype='
  myTable.add_column(Column(data=empty,name='Vls_src', dtype=	  myTable.add_column(Column(data=empty,name='Vls_src', dtype=
  myTable.add_column(Column(data=empty,name='objname', dtype=	  myTable.add_column(Column(data=empty,name='objname', dtype=

  myTable.add_column(Column(data=empty,name='DM', format='%0.	  myTable.add_column(Column(data=empty,name='DM', format='%0.
  myTable.add_column(Column(data=empty,name='eDM', format='%0	  myTable.add_column(Column(data=empty,name='eDM', format='%0
  myTable.add_column(Column(data=empty,name='mDM', format='%0	  myTable.add_column(Column(data=empty,name='mDM', format='%0
  myTable.add_column(Column(data=empty,name='meDM', format='%	  myTable.add_column(Column(data=empty,name='meDM', format='%
  								  
  								  
  #print "# of all groups: ", NoGroups				  #print "# of all groups: ", NoGroups
  								  
  for i in range(0, NoGroups):  # for all groups		  for i in range(0, NoGroups):  # for all groups
    								    
   								   
    if G_list[i][0].Vls <= 3500:				    if G_list[i][0].Vls <= 3500:
       dist = G_list[i][0].mDist				       dist = G_list[i][0].mDist
       if dist == 0 : dist = G_list[i][0].Vls / H0		       if dist == 0 : dist = G_list[i][0].Vls / H0
       if dist<1: dist = 1					       if dist<1: dist = 1
       Rg =  dist * Rg_radian(G_list[i][1:])			       Rg =  dist * Rg_radian(G_list[i][1:])
       Rg_angular = dist * G_list[i][0].R_theta			       Rg_angular = dist * G_list[i][0].R_theta
       for j in range(0, len(G_list[i])):  # for all galaxies	       for j in range(0, len(G_list[i])):  # for all galaxies
        							        
        							        
        galaxy = G_list[i][j]					        galaxy = G_list[i][j]
        flag = galaxy.flag					        flag = galaxy.flag
        pgc = galaxy.id  					        pgc = galaxy.id  
        							        
        							        
        							        
        							        
        gl = galaxy.gl  					        gl = galaxy.gl  
        gb = galaxy.gb						        gb = galaxy.gb
        sgl = galaxy.sgl  					        sgl = galaxy.sgl  
        sgb = galaxy.sgb  					        sgb = galaxy.sgb  
        Vls = galaxy.Vls  					        Vls = galaxy.Vls  
        Vhelio = galaxy.Vhelio					        Vhelio = galaxy.Vhelio
        logK = galaxy.logK  					        logK = galaxy.logK  
        Ks = galaxy.Ks						        Ks = galaxy.Ks
        Ty = galaxy.Ty  					        Ty = galaxy.Ty  
        dcf2 = galaxy.dcf2Copy					        dcf2 = galaxy.dcf2Copy
        ed = galaxy.edCopy					        ed = galaxy.edCopy
        dm = galaxy.dmCopy					        dm = galaxy.dmCopy
        edm = galaxy.edmCopy					        edm = galaxy.edmCopy
        mdm = galaxy.mDM					        mdm = galaxy.mDM
        medm = galaxy.meDM					        medm = galaxy.meDM
        							        
        							        
        ra = galaxy.ra						        ra = galaxy.ra
        dec = galaxy.dec					        dec = galaxy.dec
        coordinate_src = galaxy.coordinate_src			        coordinate_src = galaxy.coordinate_src
        Ty_src = galaxy.Ty_src					        Ty_src = galaxy.Ty_src
        B_mag = galaxy.B_mag					        B_mag = galaxy.B_mag
        Ks_src = galaxy.Ks_src					        Ks_src = galaxy.Ks_src
        Vls_src = galaxy.Vls_src				        Vls_src = galaxy.Vls_src
        objname = galaxy.objname				        objname = galaxy.objname
        							        
        mDist = galaxy.mDist					        mDist = galaxy.mDist
        mDistErr = galaxy.mDistErr				        mDistErr = galaxy.mDistErr
        							        

        							        
        subGalaxies = G_list[i][0].subGalaxies  		        subGalaxies = G_list[i][0].subGalaxies  
        R_theta = G_list[i][0].R_theta  			        R_theta = G_list[i][0].R_theta  
        sigmaP_dyn = G_list[i][0].sigma  			        sigmaP_dyn = G_list[i][0].sigma  
        nest = G_list[i][0].nest  				        nest = G_list[i][0].nest  
        Mv_lum = galaxy.M_v2					        Mv_lum = galaxy.M_v2
        R2t_lum = galaxy.R_2t2					        R2t_lum = galaxy.R_2t2
        							        
        							        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)			        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqr	        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqr
        if j != 0: 						        if j != 0: 
	  Rg = 0  # dynamical virial radius				  Rg = 0  # dynamical virial radius
	  R_theta = 0							  R_theta = 0
	  Rg_angular = 0						  Rg_angular = 0
        Mv_dyn = (1.E9 * (2.5*pi/G/2.) * (Mpc_km*Rg) * sigmaP	        Mv_dyn = (1.E9 * (2.5*pi/G/2.) * (Mpc_km*Rg) * sigmaP
        							        
        if sigmaP_dyn == 0: 					        if sigmaP_dyn == 0: 
	  tX_dyn = 0 							  tX_dyn = 0 
	else:								else:
          tX_dyn = Mpc_km*0.5*pi*Rg/sigmaP_dyn/sqrt(2.5)	          tX_dyn = Mpc_km*0.5*pi*Rg/sigmaP_dyn/sqrt(2.5)
        R2t_dyn = (Rg*pi*0.5)/1.05/sqrt(1.5)			        R2t_dyn = (Rg*pi*0.5)/1.05/sqrt(1.5)
       								       
 								 
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B	        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_l		       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_l
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2		          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2
	            coordinate_src, Ty_src , Ks_src, Vls_src,		            coordinate_src, Ty_src , Ks_src, Vls_src,
									

  								  
  # writing individual galaxies					  # writing individual galaxies
  for galaxy in galList:					  for galaxy in galList:
    								    


      if galaxy.inGroup <= 0. and galaxy.Vls<=3500:		      if galaxy.inGroup <= 0. and galaxy.Vls<=3500:
	flag = galaxy.inGroup						flag = galaxy.inGroup
	galaxy.flag = galaxy.inGroup					galaxy.flag = galaxy.inGroup
        pgc = galaxy.id  					        pgc = galaxy.id  
        gl = galaxy.gl  					        gl = galaxy.gl  
        gb = galaxy.gb						        gb = galaxy.gb
        sgl = galaxy.sgl  					        sgl = galaxy.sgl  
        sgb = galaxy.sgb  					        sgb = galaxy.sgb  
        Vls = galaxy.Vls  					        Vls = galaxy.Vls  
        Vhelio = galaxy.Vhelio					        Vhelio = galaxy.Vhelio
        logK = galaxy.logK  					        logK = galaxy.logK  
        Ks = galaxy.Ks						        Ks = galaxy.Ks
        Ty = galaxy.Ty  					        Ty = galaxy.Ty  
        dcf2 = galaxy.dcf2Copy					        dcf2 = galaxy.dcf2Copy
        ed = galaxy.edCopy					        ed = galaxy.edCopy
        dm = galaxy.dmCopy					        dm = galaxy.dmCopy
        edm = galaxy.edmCopy					        edm = galaxy.edmCopy
        mdm = galaxy.mDM					        mdm = galaxy.mDM
        medm = galaxy.meDM					        medm = galaxy.meDM
        							        
        ra = galaxy.ra						        ra = galaxy.ra
        dec = galaxy.dec					        dec = galaxy.dec
        coordinate_src = galaxy.coordinate_src			        coordinate_src = galaxy.coordinate_src
        Ty_src = galaxy.Ty_src					        Ty_src = galaxy.Ty_src
        B_mag = galaxy.B_mag					        B_mag = galaxy.B_mag
        Ks_src = galaxy.Ks_src					        Ks_src = galaxy.Ks_src
        Vls_src = galaxy.Vls_src				        Vls_src = galaxy.Vls_src
        objname = galaxy.objname				        objname = galaxy.objname
        							        
        mDist = galaxy.mDist					        mDist = galaxy.mDist
        mDistErr = galaxy.mDistErr				        mDistErr = galaxy.mDistErr

        							        
        subGalaxies = galaxy.subGalaxies  			        subGalaxies = galaxy.subGalaxies  
        R_theta = galaxy.R_theta  				        R_theta = galaxy.R_theta  
        sigmaP_dyn = galaxy.sigma  				        sigmaP_dyn = galaxy.sigma  
        nest = galaxy.nest  					        nest = galaxy.nest  
        Mv_lum = galaxy.M_v2					        Mv_lum = galaxy.M_v2
        R2t_lum = galaxy.R_2t2					        R2t_lum = galaxy.R_2t2
        							        
        							        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)			        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqr	        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqr
        Rg = 0  # dynamical virial radius			        Rg = 0  # dynamical virial radius
        Mv_dyn = 0  # solar mass				        Mv_dyn = 0  # solar mass
        tX_dyn = 0						        tX_dyn = 0
        R2t_dyn = 0						        R2t_dyn = 0
        Rg_angular = 0						        Rg_angular = 0
        R_theta = 0 						        R_theta = 0 
       								       
        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B	        myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_l		       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_l
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2		          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2
	            coordinate_src, Ty_src , Ks_src, Vls_src,		            coordinate_src, Ty_src , Ks_src, Vls_src,

									
  								  
  pgc = 999999999; 						  pgc = 999999999; 
  ra = 999.9999; dec=-99.99;					  ra = 999.9999; dec=-99.99;
  gl = ra; gb = dec						  gl = ra; gb = dec
  sgl = ra; sgb=dec						  sgl = ra; sgb=dec
  Ty = -100000.00; B_mag=Ty					  Ty = -100000.00; B_mag=Ty
  Ks= 99.99							  Ks= 99.99
  logK = 99.9999						  logK = 99.9999
  Vls = 9999							  Vls = 9999
  dcf2 = 99.99							  dcf2 = 99.99
  ed = 9.99							  ed = 9.99
  Mv_dyn = 9.99E99; Mv_lum = Mv_dyn				  Mv_dyn = 9.99E99; Mv_lum = Mv_dyn
  tX_dyn = Mv_lum; tX_lum=Mv_lum				  tX_dyn = Mv_lum; tX_lum=Mv_lum
  nest = 9999999						  nest = 9999999
  dm = 99.99							  dm = 99.99
  edm = 99.99							  edm = 99.99
  mdm = 99.99							  mdm = 99.99
  medm = 99.99							  medm = 99.99

  								  
  myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,K	  myTable.add_row([pgc,flag,ra,dec,gl, gb, sgl,sgb,Ty,B_mag,K
	       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_l		       mDist, mDistErr, R_theta, sigmaP_dyn, sigmaP_l
	          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2		          Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2
	            coordinate_src, Ty_src , Ks_src, Vls_src,		            coordinate_src, Ty_src , Ks_src, Vls_src,
  								  
  								  
  myTable.write(outfile, format='ascii.fixed_width',delimiter	  myTable.write(outfile, format='ascii.fixed_width',delimiter
  								  
  # removing the last line, (it sits o adjust the column wodt	  # removing the last line, (it sits o adjust the column wodt
  command =  ["csh", "remove_lastline.csh", outfile]		  command =  ["csh", "remove_lastline.csh", outfile]
  subprocess.call(command)					  subprocess.call(command)
#############################################################	#############################################################
#############################################################	#############################################################

def Theta_max(root):						def Theta_max(root):
  								  
  if root.left == None: return 0				  if root.left == None: return 0
  								  
  galList = NodeLeaf(root)					  galList = NodeLeaf(root)
  								  
  N = len(galList)						  N = len(galList)
  theta = np.zeros(N-1) 					  theta = np.zeros(N-1) 
  								  
  for i in range(1,N):						  for i in range(1,N):
      								      
      theta[i-1] = angle(root.sgl, root.sgb, galList[i].sgl, 	      theta[i-1] = angle(root.sgl, root.sgb, galList[i].sgl, 
  								  
  return np.max(theta)						  return np.max(theta)


#############################################################	#############################################################
# gets a list of galaxies and returns Rg in terms of radian	# gets a list of galaxies and returns Rg in terms of radian
def Rg_radian(galList):						def Rg_radian(galList):
  								  
  								  
  								  
  								  
  								  
  N = len(galList)						  N = len(galList)
  sum = 0.							  sum = 0.
  BOL = True							  BOL = True
  Rg = 0 							  Rg = 0 
  for i in range(0,N-1):					  for i in range(0,N-1):
    for j in range(i+1,N):					    for j in range(i+1,N):
      								      
      distij = angle(galList[i].sgl, galList[i].sgb, galList[	      distij = angle(galList[i].sgl, galList[i].sgb, galList[
      if distij !=0 : 						      if distij !=0 : 
         sum += 1. / distij					         sum += 1. / distij
      else: BOL = False    # If one of the didtnaces is zero,	      else: BOL = False    # If one of the didtnaces is zero,
      								      
  if BOL == True and sum != 0:					  if BOL == True and sum != 0:
    Rg = ((N)*(N)) / sum					    Rg = ((N)*(N)) / sum
  								  
  return Rg							  return Rg




#############################################################	#############################################################
#############################################################	#############################################################

        							        
#############################################################	#############################################################
def mergGroup(G_list, restrict=False):				def mergGroup(G_list, restrict=False):
    								    

    								    
    NoGroups = len(G_list)					    NoGroups = len(G_list)
 								 
    i = 0							    i = 0
    while i < NoGroups:						    while i < NoGroups:
      j = i+1							      j = i+1
      while j < NoGroups-1:					      while j < NoGroups-1:
									
        n1 = G_list[i][0].subGalaxies				        n1 = G_list[i][0].subGalaxies
        n2 = G_list[j][0].subGalaxies				        n2 = G_list[j][0].subGalaxies

	d1 = G_list[i][0].Vls/H0					d1 = G_list[i][0].Vls/H0
	d2 = G_list[j][0].Vls/H0					d2 = G_list[j][0].Vls/H0
	if d1 < 1 : d1 = 1 						if d1 < 1 : d1 = 1 
	if d2 < 1 : d2 = 1						if d2 < 1 : d2 = 1
	r1 = (180.*atan(G_list[i][0].R_2t2/d1)/pi)			r1 = (180.*atan(G_list[i][0].R_2t2/d1)/pi)
	r2 = (180.*atan(G_list[j][0].R_2t2/d2)/pi)			r2 = (180.*atan(G_list[j][0].R_2t2/d2)/pi)


	ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0]		ang12 = (180./pi)*angle(G_list[i][0].gl, G_list[i][0]
									
									
	d1 = G_list[i][0].mDist						d1 = G_list[i][0].mDist
	e1 = d1 * G_list[i][0].mDistErr					e1 = d1 * G_list[i][0].mDistErr
	d2 = G_list[j][0].mDist						d2 = G_list[j][0].mDist
	e2 = d2 * G_list[j][0].mDistErr					e2 = d2 * G_list[j][0].mDistErr
	delt = abs(d1-d2)						delt = abs(d1-d2)
									
	v1 = G_list[i][0].Vls						v1 = G_list[i][0].Vls
	v2 = G_list[j][0].Vls						v2 = G_list[j][0].Vls

									
									
	# using dynamical velocity dispersions calculated bas		# using dynamical velocity dispersions calculated bas
        sig1 = (G_list[i][0].M_v2 / (2.0E6))**(1./3)		        sig1 = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
	sig2 = (G_list[j][0].M_v2 / (2.0E6))**(1./3)			sig2 = (G_list[j][0].M_v2 / (2.0E6))**(1./3)
									

	n1 = G_list[i][0].subGalaxies					n1 = G_list[i][0].subGalaxies
	n2 = G_list[j][0].subGalaxies					n2 = G_list[j][0].subGalaxies
	Bol = False							Bol = False
									

#############################################################	#############################################################
        Sigquad = sqrt(sig1**2+sig2**2)				        Sigquad = sqrt(sig1**2+sig2**2)
#############################################################	#############################################################

	idd1 = G_list[i][0].nest					idd1 = G_list[i][0].nest
	idd2 = G_list[j][0].nest					idd2 = G_list[j][0].nest
        							        
        if restrict:						        if restrict:

	  if (idd1==13418 and  idd2==14077):				  if (idd1==13418 and  idd2==14077):
	      Bol = True						      Bol = True
	  if (idd1==36188 and  idd2==36136) or (idd1==36136 a		  if (idd1==36188 and  idd2==36136) or (idd1==36136 a
	      Bol = True						      Bol = True
	  if (idd1==13324 and  idd2==13108) or (idd1==13108 a		  if (idd1==13324 and  idd2==13108) or (idd1==13108 a
	      Bol = True    						      Bol = True    
	  if (idd1==36699 and  idd2==36875) or (idd1==36875 a		  if (idd1==36699 and  idd2==36875) or (idd1==36875 a
	      Bol = True  						      Bol = True  
	      								      

	      								      
#############################################################	#############################################################
      								      
	if ang12 <= 1.1*(r1+r2) and max(r1,r2)<6:			if ang12 <= 1.1*(r1+r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(sig1,sig2) and (min(r1,r2))**3		  if abs(v1-v2) <= max(sig1,sig2) and (min(r1,r2))**3
	       Bol = True    						       Bol = True    
	       								       
	if ang12 <= 1.0*(r1+r2) and max(r1,r2)<6 and min(n1,n		if ang12 <= 1.0*(r1+r2) and max(r1,r2)<6 and min(n1,n
	  if abs(v1-v2) <= 2.0*max(sig1,sig2):				  if abs(v1-v2) <= 2.0*max(sig1,sig2):
	      Bol = True						      Bol = True
	      								      
	if ang12 <= 0.6*(r1+r2) and max(r1,r2)<6:			if ang12 <= 0.6*(r1+r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):		  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True						      Bol = True

	      								      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:			if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):		  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True						      Bol = True
									
	      								      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:			if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt < min(e1, e2):			  if d1!=0 and d2!=0 and delt < min(e1, e2):
	      Bol = True						      Bol = True
	      								      
									
	  								  
	# one group completely projected on another one (e.g.		# one group completely projected on another one (e.g.
	if ang12 <= 1.1*(max(r1,r2)-min(r1,r2)) and max(r1,r2		if ang12 <= 1.1*(max(r1,r2)-min(r1,r2)) and max(r1,r2
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):		  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True	      					      Bol = True	      
#############################################################	#############################################################

        if (idd1==39241 and  idd2==39600) or (idd1==39600 and	        if (idd1==39241 and  idd2==39600) or (idd1==39600 and
	  Bol = False							  Bol = False
	if (idd1==13324 and  idd2==13108) or (idd1==13108 and		if (idd1==13324 and  idd2==13108) or (idd1==13108 and
	  Bol = True    						  Bol = True    
	  								  

	if (idd1==12651 and  idd2==12651):				if (idd1==12651 and  idd2==12651):
	  Bol = False							  Bol = False
	if (idd1==13418 and  idd2==12651) or (idd1==12651 and		if (idd1==13418 and  idd2==12651) or (idd1==12651 and
	  Bol = False  							  Bol = False  

	  								  
        if Bol:							        if Bol:
	  								  
	  newGroup = []							  newGroup = []
	  								  
	  for p in range(1, len(G_list[i])):				  for p in range(1, len(G_list[i])):
	    newGroup.append(G_list[i][p])				    newGroup.append(G_list[i][p])
	  								  
	  for q in range(1, len(G_list[j])):				  for q in range(1, len(G_list[j])):
	    newGroup.append(G_list[j][q])				    newGroup.append(G_list[j][q])
	  								  
 								 
	  root = LoosgalListJoint(newGroup, grID = 200000000)		  root = LoosgalListJoint(newGroup, grID = 200000000)
	  G_list[i] = NodeLeaf(root)					  G_list[i] = NodeLeaf(root)

	  								  
									
	  G_list.pop(j) 						  G_list.pop(j) 
	  NoGroups-=1							  NoGroups-=1
	  i=0								  i=0
	  j=0								  j=0
        j+=1							        j+=1
      i+=1							      i+=1
    								    
    return G_list						    return G_list

 ############################################################	 ############################################################
 # Function name: addGalGroup					 # Function name: addGalGroup
 #								 #
 # This function tries to find those galaxies which do not fa	 # This function tries to find those galaxies which do not fa
 # into any Group. If the galaxy is close enough to the cente	 # into any Group. If the galaxy is close enough to the cente
 # according to its radial velocity (i.e. within 2*sigma of t	 # according to its radial velocity (i.e. within 2*sigma of t
 # it might get absorbed by that group...			 # it might get absorbed by that group...
 ############################################################	 ############################################################

def addGalGroup(G_list, galList):				def addGalGroup(G_list, galList):
  								  

    singles = []						    singles = []
    for galaxy in galList:					    for galaxy in galList:
       if galaxy.inGroup == 0:					       if galaxy.inGroup == 0:
	  singles.append([galaxy, -1, -10000])  # [galaxy, ne		  singles.append([galaxy, -1, -10000])  # [galaxy, ne
    								    
    N = len(G_list)						    N = len(G_list)
    #BB = []							    #BB = []
    #LL = []							    #LL = []
    #for i in range(N):						    #for i in range(N):
      #LL.append(G_list[i][0].gl)				      #LL.append(G_list[i][0].gl)
      #BB.append(G_list[i][0].gb)				      #BB.append(G_list[i][0].gb)
    								    
    for entity in singles:					    for entity in singles:
        if  entity[0].id != 40001:				        if  entity[0].id != 40001:
            gl = entity[0].gl					            gl = entity[0].gl
            gb = entity[0].gb					            gb = entity[0].gb
    								    
	    #q1 = np.zeros((N,), dtype=np.int)				    #q1 = np.zeros((N,), dtype=np.int)
            #q2 = np.zeros((N,), dtype=np.int)			            #q2 = np.zeros((N,), dtype=np.int)
            #q3 = np.zeros((N,), dtype=np.int)			            #q3 = np.zeros((N,), dtype=np.int)
            #q4 = np.zeros((N,), dtype=np.int)			            #q4 = np.zeros((N,), dtype=np.int)
	    #q1[np.where(LL<=gl+10)] = 1				    #q1[np.where(LL<=gl+10)] = 1
	    #q2[np.where(LL>=gl-10)] = 1				    #q2[np.where(LL>=gl-10)] = 1
	    #q3[np.where(BB<=gb+10)] = 1				    #q3[np.where(BB<=gb+10)] = 1
	    #q4[np.where(BB>=gb-10)] = 1				    #q4[np.where(BB>=gb-10)] = 1
	    #qq = q1+q2+q3+q4						    #qq = q1+q2+q3+q4
	    #group_indices = np.where(qq==4)				    #group_indices = np.where(qq==4)
	       								       
	    for i in range(len(G_list)): # group_indices[0]:		    for i in range(len(G_list)): # group_indices[0]:
	      if  G_list[i][0].nest != 41220:  # Virgo			      if  G_list[i][0].nest != 41220:  # Virgo
		    		    						    		    
		    ang12 = (180./pi)*angle(G_list[i][0].gl, 			    ang12 = (180./pi)*angle(G_list[i][0].gl, 
		    d = G_list[i][0].Vls / H0					    d = G_list[i][0].Vls / H0
		    if d < 1: d = 1 						    if d < 1: d = 1 
		    r = (180.*atan(G_list[i][0].R_2t2/(d))/pi			    r = (180.*atan(G_list[i][0].R_2t2/(d))/pi
		    								    
		    # if the galaxy is close enough to the ce			    # if the galaxy is close enough to the ce
		    # and it's not already in any other group			    # and it's not already in any other group
		    sig_p = (G_list[i][0].M_v2 / (2.0E6))**(1			    sig_p = (G_list[i][0].M_v2 / (2.0E6))**(1
		    								    
		    								    
		    #if G_list[i][0].nest == 40692: 				    #if G_list[i][0].nest == 40692: 
		      #print "Esn TEst: ", G_list[i][0].Vls, 			      #print "Esn TEst: ", G_list[i][0].Vls, 
		    								    
		    if ang12 <= 2*r:						    if ang12 <= 2*r:

		      								      
		      								      
											
											
		      if ang12 <= 1.05*(180./pi)*atan(1.3*G_l			      if ang12 <= 1.05*(180./pi)*atan(1.3*G_l
			  d1 = G_list[i][0].mDist					  d1 = G_list[i][0].mDist
			  e1 = d1 * G_list[i][0].mDistErr				  e1 = d1 * G_list[i][0].mDistErr
			  d2 = entity[0].dcf2						  d2 = entity[0].dcf2
			  e2 = d2 * entity[0].ed					  e2 = d2 * entity[0].ed
			  delt = abs(d1-d2)						  delt = abs(d1-d2)
			  join = False							  join = False
			  if ang12 > 1.05*(180./pi)*atan(1.0*				  if ang12 > 1.05*(180./pi)*atan(1.0*
			    join = True							    join = True
			  if d1==0 or d2==0:						  if d1==0 or d2==0:
			    join = True							    join = True
			  if ang12 <= 1.05*(180./pi)*atan(1.0				  if ang12 <= 1.05*(180./pi)*atan(1.0
			    join = True							    join = True
			    								    

			  if join==True and  sig_p > entity[2				  if join==True and  sig_p > entity[2
			    entity[1] = i						    entity[1] = i
			    entity[2] = sig_p						    entity[2] = sig_p
		  								  
            if entity[1] > -1: # if the galaxy is absorbed (c	            if entity[1] > -1: # if the galaxy is absorbed (c
              entity[0].inGroup = 1				              entity[0].inGroup = 1
              G_list[entity[1]].append(entity[0])		              G_list[entity[1]].append(entity[0])
    for i in range(N):						    for i in range(N):
       G_list[i].pop(0)						       G_list[i].pop(0)
       root  = LoosgalListJoint(G_list[i], grID = 300000000)	       root  = LoosgalListJoint(G_list[i], grID = 300000000)
       G_list[i] = NodeLeaf(root) 				       G_list[i] = NodeLeaf(root) 
    return G_list						    return G_list


#############################################################	#############################################################
# Sometimes when groups grow the order of galaxy absorbation 	# Sometimes when groups grow the order of galaxy absorbation 
def IndividualGal_modifier(G_list, galList):			def IndividualGal_modifier(G_list, galList):
  								  

  								  

  tmp_galList = []						  tmp_galList = []
  N = len(galList)						  N = len(galList)
  for galaxy in galList:					  for galaxy in galList:
    if galaxy.inGroup > 0:					    if galaxy.inGroup > 0:
      galaxy.inGroup = 0					      galaxy.inGroup = 0
      tmp_galList.append(galaxy)				      tmp_galList.append(galaxy)
    								    

  new_G_list=[]							  new_G_list=[]
  for Group in G_list:						  for Group in G_list:
    Group[1].inGroup = 1  # first galaxy of group, dominant o	    Group[1].inGroup = 1  # first galaxy of group, dominant o
    new_G_list.append([Group[0], Group[1]])			    new_G_list.append([Group[0], Group[1]])
  								  

  								  
  								  
  Lum_mass = []							  Lum_mass = []
  for i in range(0, len(new_G_list)): Lum_mass.append(new_G_l	  for i in range(0, len(new_G_list)): Lum_mass.append(new_G_l
  Lum_mass = np.asarray(Lum_mass)				  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)				  indices = np.argsort(Lum_mass)
  								  
  NEWgpLIST = []						  NEWgpLIST = []
  for i in indices[::-1]:					  for i in indices[::-1]:
            							            
            							            
            grp = new_G_list[i][0]				            grp = new_G_list[i][0]
            gl = grp.gl						            gl = grp.gl
            gb = grp.gb						            gb = grp.gb
    								    

	    j = 0 							    j = 0 
	    N = len(tmp_galList)					    N = len(tmp_galList)
	    while j < N: 						    while j < N: 
                     						                     
                  gal = tmp_galList[j]				                  gal = tmp_galList[j]

	          ang12 = (180./pi)*angle(gl, gb, gal.gl, gal		          ang12 = (180./pi)*angle(gl, gb, gal.gl, gal
                  d = grp.Vls / H0				                  d = grp.Vls / H0
                  if d < 1: d = 1 				                  if d < 1: d = 1 
                  r = (180.*atan(0.95*grp.R_2t2/(d))/pi)	                  r = (180.*atan(0.95*grp.R_2t2/(d))/pi)
                  						                  
                  if ang12 <= r and gal.inGroup == 0:		                  if ang12 <= r and gal.inGroup == 0:
 								 
	              sig_p = (grp.M_v2 / (2.0E6))**(1./3)		              sig_p = (grp.M_v2 / (2.0E6))**(1./3)
	              if (abs(grp.Vls-gal.Vls) <= 2.0*sig_p):		              if (abs(grp.Vls-gal.Vls) <= 2.0*sig_p):
			new_G_list[i].append(gal)					new_G_list[i].append(gal)
			gal.inGroup = 1							gal.inGroup = 1
		   	tmp_galList.pop(j)					   	tmp_galList.pop(j)
			j -= 1								j -= 1
			N -= 1								N -= 1
	          j += 1						          j += 1

	    								    
    								    
	    new_G_list[i].pop(0)					    new_G_list[i].pop(0)
	    if len(new_G_list[i]) > 1: 					    if len(new_G_list[i]) > 1: 
	        grp =  LoosgalListJoint(new_G_list[i], grID =		        grp =  LoosgalListJoint(new_G_list[i], grID =
	        NEWgpLIST.append(NodeLeaf(grp))				        NEWgpLIST.append(NodeLeaf(grp))
	    else:							    else:
	        new_G_list[i][0].inGroup = 0				        new_G_list[i][0].inGroup = 0
	        							        
  								  
  virg_w = None							  virg_w = None
  for group in NEWgpLIST:					  for group in NEWgpLIST:
    for galaxy in group[1:]:					    for galaxy in group[1:]:
      if galaxy.id == 39659:					      if galaxy.id == 39659:
	virg_w = group							virg_w = group
        break							        break
  								  
  if virg_w != None:						  if virg_w != None:
    newG = []							    newG = []
    #print "Here I am"						    #print "Here I am"
    for galaxy in virg_w[1:]:					    for galaxy in virg_w[1:]:
      if galaxy.Vls < 1450 and galaxy.sgb > -7:			      if galaxy.Vls < 1450 and galaxy.sgb > -7:
	galaxy.inGroup = 0						galaxy.inGroup = 0
      else:							      else:
        newG.append(galaxy)					        newG.append(galaxy)
   								   
    grp = LoosgalListJoint(newG, grID = 400000000)		    grp = LoosgalListJoint(newG, grID = 400000000)
    virg_w = NodeLeaf(grp)					    virg_w = NodeLeaf(grp)
  								  
  								  
  								  
  # This is the output	    					  # This is the output	    
  #G_list = NEWgpLIST						  #G_list = NEWgpLIST
  return NEWgpLIST						  return NEWgpLIST


#############################################################	#############################################################
## Trying to find linkages when both galaxies have distances	## Trying to find linkages when both galaxies have distances

def find_pair_dist(galList_org):				def find_pair_dist(galList_org):
  								  
  galList = []							  galList = []
  for galaxy in galList_org:					  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 		    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 
	galList.append(galaxy)  					galList.append(galaxy)  
  								  
  individuals = []						  individuals = []
  Lum_mass = []							  Lum_mass = []
  for i in range(0, len(galList)): Lum_mass.append(galList[i]	  for i in range(0, len(galList)): Lum_mass.append(galList[i]
  Lum_mass = np.asarray(Lum_mass)				  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)				  indices = np.argsort(Lum_mass)
  for i in indices[::-1]: individuals.append(galList[i])  	  for i in indices[::-1]: individuals.append(galList[i])  
  								  
  NEWgpLIST =[]							  NEWgpLIST =[]
  p = 0 							  p = 0 
  N = len(individuals)						  N = len(individuals)
  while p < N-1:						  while p < N-1:
    #print p, N							    #print p, N
    galaxies = []						    galaxies = []
    galaxies.append(individuals[p])				    galaxies.append(individuals[p])

    q = (p+1)%N							    q = (p+1)%N
    pivot = p							    pivot = p
    grp = galaxies[0]						    grp = galaxies[0]
    Bol = False							    Bol = False
    								    
    if grp.dcf2==0:						    if grp.dcf2==0:
      p+=1							      p+=1
      continue							      continue
    								    
    while q!=pivot:						    while q!=pivot:
      								      
      if individuals[q].dcf2==0:				      if individuals[q].dcf2==0:
	q = (q+1)%N							q = (q+1)%N
	continue							continue
      								      
      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].	      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].
      								      
      d = grp.dcf2						      d = grp.dcf2
      if d == 0:						      if d == 0:
	q = (q+1)%N							q = (q+1)%N
	continue							continue
      								      
      								      
      mass = Mass(10**m_logK(grp.Ks, grp.Vls, d=d))		      mass = Mass(10**m_logK(grp.Ks, grp.Vls, d=d))
      R_2t2 = 0.215*((mass/1.E12)**(1./3))			      R_2t2 = 0.215*((mass/1.E12)**(1./3))

      thet = (180.*atan(1.0*R_2t2/d)/pi)			      thet = (180.*atan(1.0*R_2t2/d)/pi)

      coeff = 1.3						      coeff = 1.3

      if abs(individuals[q].dcf2 - d)/d < 0.1 and ang12 < coe	      if abs(individuals[q].dcf2 - d)/d < 0.1 and ang12 < coe
            							            
            mass = Mass(10**m_logK(individuals[q].Ks, individ	            mass = Mass(10**m_logK(individuals[q].Ks, individ
            R_2t2_2 = 0.215*((mass/1.E12)**(1./3))		            R_2t2_2 = 0.215*((mass/1.E12)**(1./3))
            							            
            Bol = True						            Bol = True

	    galaxies.append(individuals[q])				    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 40000000		    grp =  LoosgalListJoint(galaxies, grID = 40000000
	    individuals.pop(q)						    individuals.pop(q)

	    N-=1							    N-=1
	    								    
	    if p>=q:  p-=1						    if p>=q:  p-=1
	    								    
	    								    
	    if len(galaxies) == 2:					    if len(galaxies) == 2:
	      individuals.pop(p)					      individuals.pop(p)
	      N-=1							      N-=1	
	      if q>p: q = (q-1)%N					      if q>p: q = (q-1)%N
	      break							      break
            							            
            							            
	    q = (q-1)%N							    q = (q-1)%N
	    pivot = q 							    pivot = q 
	    								    
      q = (q+1)%N						      q = (q+1)%N
    if Bol: 							    if Bol: 
      NEWgpLIST.append(NodeLeaf(grp))				      NEWgpLIST.append(NodeLeaf(grp))
    p+=1							    p+=1
    								    
    for i in range(0, len(NEWgpLIST)):				    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):			      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1				         NEWgpLIST[i][j].inGroup = 1
  								  
  								  
  return NEWgpLIST  						  return NEWgpLIST  


#############################################################	#############################################################

def IndividualGal(galList_org, pairs=False):			def IndividualGal(galList_org, pairs=False):
  								  
  								  
  galList = []							  galList = []
  for galaxy in galList_org:					  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 		    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 
	galList.append(galaxy)   					galList.append(galaxy)   
  								  

  individuals = []						  individuals = []
  Lum_mass = []							  Lum_mass = []
  for i in range(0, len(galList)): Lum_mass.append(galList[i]	  for i in range(0, len(galList)): Lum_mass.append(galList[i]
  Lum_mass = np.asarray(Lum_mass)				  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)				  indices = np.argsort(Lum_mass)
  for i in indices[::-1]: individuals.append(galList[i])  	  for i in indices[::-1]: individuals.append(galList[i])  
  								  
  NEWgpLIST =[]							  NEWgpLIST =[]
  p = 0 							  p = 0 
  N = len(individuals)						  N = len(individuals)
  while p < N-1:						  while p < N-1:
    #print p, N							    #print p, N
    galaxies = []						    galaxies = []
    galaxies.append(individuals[p])				    galaxies.append(individuals[p])
    #individuals.pop(p)						    #individuals.pop(p)
    #N-=1							    #N-=1
    q = (p+1)%N							    q = (p+1)%N
    pivot = p							    pivot = p
    grp = galaxies[0]						    grp = galaxies[0]
    Bol = False							    Bol = False
    while q!=pivot:						    while q!=pivot:
      								      
      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].	      ang12 = (180./pi)*angle(grp.gl, grp.gb, individuals[q].
      								      
      d = grp.Vls/H0						      d = grp.Vls/H0
      if d <=1: d =1						      if d <=1: d =1

      thet = (180.*atan(1.0*grp.R_2t2/d)/pi)			      thet = (180.*atan(1.0*grp.R_2t2/d)/pi)
      #if grp.Vls<75:						      #if grp.Vls<75:
	#thet = (180.*atan(1.0*grp.R_2t2/(1.))/pi)			#thet = (180.*atan(1.0*grp.R_2t2/(1.))/pi)
      								      
      coeff = 1.05						      coeff = 1.05
      delta_v = abs(grp.Vls-individuals[q].Vls)			      delta_v = abs(grp.Vls-individuals[q].Vls)
      								      
      								      
      test = False						      test = False
      if pairs:							      if pairs:
	sig_sig = 100000						sig_sig = 100000
	if len(galaxies) == 1:						if len(galaxies) == 1:
	  d = individuals[q].Vls/H0					  d = individuals[q].Vls/H0
	  if d <=1: d =1						  if d <=1: d =1
	  thet +=  (180.*atan(individuals[q].R_2t2/d)/pi)		  thet +=  (180.*atan(individuals[q].R_2t2/d)/pi)
	  L1 = 10**grp.logK						  L1 = 10**grp.logK
	  L2 = 10**individuals[q].logK					  L2 = 10**individuals[q].logK
	  L_tot = L1 + L2						  L_tot = L1 + L2
	  mass = Mass(L_tot)						  mass = Mass(L_tot)
	  sig_sig = (mass / (2.0E6))**(1./3)				  sig_sig = (mass / (2.0E6))**(1./3)
	  if delta_v < 2.0*sig_sig: 					  if delta_v < 2.0*sig_sig: 
	    test = True							    test = True
									
      								      
      								      
      if ang12 <= coeff*thet:					      if ang12 <= coeff*thet:
	sig_p = (grp.M_v2 / (2.0E6))**(1./3)				sig_p = (grp.M_v2 / (2.0E6))**(1./3)

	if test == True or delta_v <= 2.0*sig_p:			if test == True or delta_v <= 2.0*sig_p:

	    Bol = True							    Bol = True
	    								    

	    galaxies.append(individuals[q])				    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 40000000		    grp =  LoosgalListJoint(galaxies, grID = 40000000
	    individuals.pop(q)						    individuals.pop(q)

	    N-=1							    N-=1
	    								    
	    if p>=q:  p-=1						    if p>=q:  p-=1
	    								    
	    								    
	    if len(galaxies) == 2:					    if len(galaxies) == 2:
	      individuals.pop(p)					      individuals.pop(p)
	      N-=1							      N-=1	
	      if q>p: q = (q-1)%N					      if q>p: q = (q-1)%N
            							            
            							            
	    q = (q-1)%N							    q = (q-1)%N
	    pivot = q 							    pivot = q 
	    								    
      #if grp.id == 40001: Bol = False				      #if grp.id == 40001: Bol = False
      q = (q+1)%N						      q = (q+1)%N
    if Bol: #len(galaxies) > 1:					    if Bol: #len(galaxies) > 1:
      NEWgpLIST.append(NodeLeaf(grp))				      NEWgpLIST.append(NodeLeaf(grp))
    p+=1							    p+=1
    								    
    for i in range(0, len(NEWgpLIST)):				    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):			      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1				         NEWgpLIST[i][j].inGroup = 1
  								  
  								  
  ##3 Dealingwith Virgo-W					  ##3 Dealingwith Virgo-W
  virg_w = None							  virg_w = None
  for group in NEWgpLIST:					  for group in NEWgpLIST:
    for galaxy in group[1:]:					    for galaxy in group[1:]:
      if galaxy.id == 39659:					      if galaxy.id == 39659:
	virg_w = group							virg_w = group
        break							        break
  								  
  if virg_w != None:						  if virg_w != None:
    newG = []							    newG = []
    #print "Here I am"						    #print "Here I am"
    for galaxy in virg_w[1:]:					    for galaxy in virg_w[1:]:
      if galaxy.Vls < 1450 and galaxy.sgb > -7:			      if galaxy.Vls < 1450 and galaxy.sgb > -7:
	galaxy.inGroup = 0						galaxy.inGroup = 0
      else:							      else:
        newG.append(galaxy)					        newG.append(galaxy)
   								   
    grp = LoosgalListJoint(newG, grID = 400000000)		    grp = LoosgalListJoint(newG, grID = 400000000)
    virg_w = NodeLeaf(grp)					    virg_w = NodeLeaf(grp)
  								  
  								  
  return NEWgpLIST						  return NEWgpLIST

#############################################################	#############################################################

def LoosgalListJoint(galList, grID = 500000000, noSort=False)	def LoosgalListJoint(galList, grID = 500000000, noSort=False)
   								   
   if len(galList)==0: return None				   if len(galList)==0: return None
   								   
   Lum_mass = []						   Lum_mass = []
   for i in range(0, len(galList)):				   for i in range(0, len(galList)):
        Lum_mass.append(galList[i].Ks)				        Lum_mass.append(galList[i].Ks)
   Lum_mass = np.asarray(Lum_mass)				   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass)				   indices = np.argsort(Lum_mass)
   								   
   if noSort:							   if noSort:
      Lum_mass = []						      Lum_mass = []
      for i in range(0, len(galList)):				      for i in range(0, len(galList)):
	    Lum_mass.append(galList[i].logK)				    Lum_mass.append(galList[i].logK)
      Lum_mass = np.asarray(Lum_mass)				      Lum_mass = np.asarray(Lum_mass)
      indices = np.argsort(Lum_mass)				      indices = np.argsort(Lum_mass)
      indices = indices[::-1]					      indices = indices[::-1]
   								   
   root = None							   root = None
   for i in indices:						   for i in indices:
       root = galJoint(root, galList[i], grID + galList[i].id	       root = galJoint(root, galList[i], grID + galList[i].id
   								   
   								   
   sumDist = 0							   sumDist = 0
   sumError = 0							   sumError = 0
   V_mean = 0 							   V_mean = 0 
   for galaxy in galList:					   for galaxy in galList:
    sumDist += galaxy.sumDist					    sumDist += galaxy.sumDist
    sumError += galaxy.sumError					    sumError += galaxy.sumError
    V_mean += galaxy.Vls					    V_mean += galaxy.Vls
   								   
   V_mean /= len(galList)					   V_mean /= len(galList)
   mDM = 0.							   mDM = 0.
   meDM = 0.							   meDM = 0.
   mDist = 0.							   mDist = 0.
   mDistErr = 0.						   mDistErr = 0.
   if sumDist !=0 and sumError!=0:				   if sumDist !=0 and sumError!=0:
      mDM = sumDist/sumError					      mDM = sumDist/sumError
      meDM = sqrt(1/sumError)					      meDM = sqrt(1/sumError)
      mDist = 10**(mDM/5.-5.)					      mDist = 10**(mDM/5.-5.)
      mDistErr = (0.2*log(10.))*meDM				      mDistErr = (0.2*log(10.))*meDM

   								   
   								   
   for gal in galList:						   for gal in galList:
      gal.logK = m_logK(gal.Ks, V_mean)				      gal.logK = m_logK(gal.Ks, V_mean)
      gal.mDist = mDist						      gal.mDist = mDist
      gal.mDistErr = mDistErr					      gal.mDistErr = mDistErr
      gal.mDM = mDM						      gal.mDM = mDM
      gal.meDM = meDM						      gal.meDM = meDM
     								     
   root = None							   root = None
   for i in indices:						   for i in indices:
      root = galJoint(root, galList[i], grID + galList[i].id)	      root = galJoint(root, galList[i], grID + galList[i].id)
   								   
   								   
   root.dcf2 = 0						   root.dcf2 = 0
   root.ed = 0							   root.ed = 0
   								   
   return root							   return root

#############################################################	#############################################################

def Theta_gr(head, galList):					def Theta_gr(head, galList):
  								  
  								  
  								  
  N = len(galList)						  N = len(galList)
  theta = np.zeros(N-1) 					  theta = np.zeros(N-1) 
  								  
  for i in range(1,N):						  for i in range(1,N):
      								      
      theta[i-1] = angle(head.l, head.b, galList[i].l, galLis	      theta[i-1] = angle(head.l, head.b, galList[i].l, galLis
  								  
  return np.max(theta)						  return np.max(theta)



#############################################################	#############################################################

def galList_moderator(galList):					def galList_moderator(galList):
  								  
  for galaxy in galList:					  for galaxy in galList:
    if galaxy.inGroup <= 0:					    if galaxy.inGroup <= 0:
      galaxy.mDist = galaxy.dcf2Copy				      galaxy.mDist = galaxy.dcf2Copy
      galaxy.mDistErr = galaxy.edCopy				      galaxy.mDistErr = galaxy.edCopy
      galaxy.mDM = galaxy.dmCopy				      galaxy.mDM = galaxy.dmCopy
      galaxy.meDM = galaxy.edmCopy				      galaxy.meDM = galaxy.edmCopy
      								      
     								     
      #if galaxy.sgl>110 and galaxy.sgl<125 and galaxy.sgb>-8	      #if galaxy.sgl>110 and galaxy.sgl<125 and galaxy.sgb>-8
	#if  galaxy.Vls<=1750: 						#if  galaxy.Vls<=1750: 
           #galaxy.setMeanDist(galaxy.mDist, galaxy.mDistErr,	           #galaxy.setMeanDist(galaxy.mDist, galaxy.mDistErr,
	#else:								#else:
	   #galaxy.setMeanDist(galaxy.mDist, galaxy.mDistErr,		   #galaxy.setMeanDist(galaxy.mDist, galaxy.mDistErr,

#############################################################	#############################################################


def group_moderator(G_list, use_dist=False):			def group_moderator(G_list, use_dist=False):

  for group in G_list:						  for group in G_list:
	  groupHeader = group[0]					  groupHeader = group[0]
	  groupHeader.flag = 2						  groupHeader.flag = 2
	  groupHeader.dcf2 = 0						  groupHeader.dcf2 = 0
	  groupHeader.ed = 0						  groupHeader.ed = 0
	  groupHeader.dcf2Copy = 0					  groupHeader.dcf2Copy = 0
          groupHeader.edCopy = 0				          groupHeader.edCopy = 0
          groupHeader.dm = 0					          groupHeader.dm = 0
          groupHeader.edm = 0					          groupHeader.edm = 0
          groupHeader.dmCopy = 0				          groupHeader.dmCopy = 0
          groupHeader.edmCopy = 0				          groupHeader.edmCopy = 0
          							          
          							          

          sum_v = 0.						          sum_v = 0.
          sum_v2 = 0.						          sum_v2 = 0.
          helio_sum_v = 0.					          helio_sum_v = 0.
          helio_sum_v2 = 0.          				          helio_sum_v2 = 0.          
          n_v = 0						          n_v = 0
          for galaxy in group[1:]:				          for galaxy in group[1:]:

	    if galaxy.Vhelio == 0:					    if galaxy.Vhelio == 0:
	      galaxy.Vls = 0						      galaxy.Vls = 0
	    if galaxy.Vls !=0:						    if galaxy.Vls !=0:
	      sum_v += galaxy.Vls					      sum_v += galaxy.Vls
	      sum_v2 += galaxy.Vls**2					      sum_v2 += galaxy.Vls**2
	      helio_sum_v += galaxy.Vhelio				      helio_sum_v += galaxy.Vhelio
	      helio_sum_v2 += galaxy.Vhelio**2				      helio_sum_v2 += galaxy.Vhelio**2
	      n_v += 1							      n_v += 1
	  								  
	  if n_v != 0:							  if n_v != 0:
	    v_av = sum_v / n_v						    v_av = sum_v / n_v
	    helio_v_av = helio_sum_v / n_v				    helio_v_av = helio_sum_v / n_v
	    								    
	    groupHeader.Vls = v_av					    groupHeader.Vls = v_av
	    groupHeader.sigma = sqrt(sum_v2/n_v - v_av**2) 		    groupHeader.sigma = sqrt(sum_v2/n_v - v_av**2) 
	    								    
	    groupHeader.Vhelio = helio_v_av				    groupHeader.Vhelio = helio_v_av

	  else:								  else:
	    groupHeader.Vls = 0						    groupHeader.Vls = 0
	    groupHeader.sigma = 0  # in Vls				    groupHeader.sigma = 0  # in Vls
	    groupHeader.Vhelio = 0					    groupHeader.Vhelio = 0
	  								  
	  								  
	  mDM = 0.							  mDM = 0.
	  meDM = 0.							  meDM = 0.
	  meanDist = 0.							  meanDist = 0.
	  meanDistErr = 0.						  meanDistErr = 0.
	  sumDist = 0							  sumDist = 0
          sumError = 0						          sumError = 0
	  for galaxy in group[1:]:					  for galaxy in group[1:]:
	    if galaxy.dcf2 != 0 and galaxy.ed != 0:			    if galaxy.dcf2 != 0 and galaxy.ed != 0:
	      if not galaxy.id in [43775, 42447, 13368]:		      if not galaxy.id in [43775, 42447, 13368]:
		    err =  galaxy.dcf2*galaxy.ed				    err =  galaxy.dcf2*galaxy.ed
		    sumDist += galaxy.dcf2/(err**2)				    sumDist += galaxy.dcf2/(err**2)
		    sumError += 1./(err**2)					    sumError += 1./(err**2)
	  								  
	  								  
	  if sumDist !=0 and sumError != 0:				  if sumDist !=0 and sumError != 0:
            meanDist = sumDist/sumError				            meanDist = sumDist/sumError
            meanDistErr = sqrt(1./sumError)/meanDist		            meanDistErr = sqrt(1./sumError)/meanDist
            							            
           							           
              							              
          if meanDist>0:					          if meanDist>0:
            groupHeader.mDM  = (log10(meanDist)+5.)*5.		            groupHeader.mDM  = (log10(meanDist)+5.)*5.
            groupHeader.meDM = meanDistErr / (0.2*log(10.))	            groupHeader.meDM = meanDistErr / (0.2*log(10.))
          else:							          else:
	    groupHeader.mDM = 0						    groupHeader.mDM = 0
	    groupHeader.meDM = 0					    groupHeader.meDM = 0
	  								  
	  groupHeader.dm  = groupHeader.mDM				  groupHeader.dm  = groupHeader.mDM
          groupHeader.edm = groupHeader.meDM			          groupHeader.edm = groupHeader.meDM
          groupHeader.mDist = meanDist				          groupHeader.mDist = meanDist
          groupHeader.mDistErr = meanDistErr			          groupHeader.mDistErr = meanDistErr
 								 
	  L_tot = 0							  L_tot = 0
	  ID = groupHeader.id 						  ID = groupHeader.id 
	  groupHeader.id = ID - ID % 100000000 + groupHeader.		  groupHeader.id = ID - ID % 100000000 + groupHeader.
	  								  
	  V_hel_sum = 0.						  V_hel_sum = 0.
	  								  
	  for galaxy in group[1:]:					  for galaxy in group[1:]:
	    galaxy.flag = 1						    galaxy.flag = 1
	    # and also modify the absolute luminosities			    # and also modify the absolute luminosities
	    if use_dist:						    if use_dist:
	       if galaxy.Vls<200: #galaxy.nest == 5064336 or 		       if galaxy.Vls<200: #galaxy.nest == 5064336 or 
	         galaxy.setMeanDist(meanDist, meanDistErr, GR		         galaxy.setMeanDist(meanDist, meanDistErr, GR
	         A= True						         A= True
	         							         
	       #elif groupHeader.sgl>110 and groupHeader.sgl<		       #elif groupHeader.sgl>110 and groupHeader.sgl<
		   #if groupHeader.mDist>0:					   #if groupHeader.mDist>0:
		     #galaxy.setMeanDist(meanDist, meanDistEr			     #galaxy.setMeanDist(meanDist, meanDistEr
		   #else:							   #else:
		     #if  groupHeader.Vls<=1750: 				     #if  groupHeader.Vls<=1750: 
		       #galaxy.setMeanDist(meanDist, meanDist			       #galaxy.setMeanDist(meanDist, meanDist
		     #else:							     #else:
		       #galaxy.setMeanDist(meanDist, meanDist			       #galaxy.setMeanDist(meanDist, meanDist
	       else: 							       else: 
		 galaxy.setMeanDist(meanDist, meanDistErr, GR			 galaxy.setMeanDist(meanDist, meanDistErr, GR
	    else:							    else:
	       galaxy.setMeanDist(meanDist, meanDistErr, GRP_		       galaxy.setMeanDist(meanDist, meanDistErr, GRP_
	   								   
	    galaxy.mDM = mDM						    galaxy.mDM = mDM
	    galaxy.meDM = meDM						    galaxy.meDM = meDM
	    L_tot += 10**galaxy.logK					    L_tot += 10**galaxy.logK
	    galaxy.mDM  = groupHeader.mDM				    galaxy.mDM  = groupHeader.mDM
	    galaxy.meDM = groupHeader.meDM				    galaxy.meDM = groupHeader.meDM
	    								    
	  								  
	  								  
          groupHeader.logK = log10(L_tot)			          groupHeader.logK = log10(L_tot)
          Dist_v =  groupHeader.Vls / H0			          Dist_v =  groupHeader.Vls / H0
          if Dist_v<1. : Dist_v=1				          if Dist_v<1. : Dist_v=1
          Mk_sun = 3.28   # K-band				          Mk_sun = 3.28   # K-band
          M = Mk_sun - (groupHeader.logK / 0.4)			          M = Mk_sun - (groupHeader.logK / 0.4)
          groupHeader.Ks = M + 5*log10(Dist_v) + 30 - 5		          groupHeader.Ks = M + 5*log10(Dist_v) + 30 - 5
          							          
          groupHeader.M_v2 = Mass(10**groupHeader.logK)		          groupHeader.M_v2 = Mass(10**groupHeader.logK)
          groupHeader.R_2t2 = 0.215*((groupHeader.M_v2/1.E12)	          groupHeader.R_2t2 = 0.215*((groupHeader.M_v2/1.E12)
  								  
  return G_list							  return G_list

#############################################################	#############################################################
#############################################################	#############################################################
#############################################################	#############################################################
# removing a galaxy from all group				# removing a galaxy from all group
def removeFrom(G_list, id):					def removeFrom(G_list, id):
  								  
  for i in range(len(G_list)):					  for i in range(len(G_list)):
    group = G_list[i]						    group = G_list[i]
    for galaxy in group[1:]:					    for galaxy in group[1:]:
      if galaxy.id == id:					      if galaxy.id == id:
        newG = []						        newG = []
        for galaxy in group[1:]:				        for galaxy in group[1:]:
           if galaxy.id == id:					           if galaxy.id == id:
        	galaxy.inGroup = 0				        	galaxy.inGroup = 0
           else:						           else:
                newG.append(galaxy)				                newG.append(galaxy)
        							        
	if len(newG) > 1 : 						if len(newG) > 1 : 
          grp = LoosgalListJoint(newG, grID = 500000000)	          grp = LoosgalListJoint(newG, grID = 500000000)
          group = NodeLeaf(grp)					          group = NodeLeaf(grp)
          G_list.pop(i)						          G_list.pop(i)
          G_list.append(group)					          G_list.append(group)
        elif  len(newG) == 1 :					        elif  len(newG) == 1 :
	  newG[0].inGroup = 0						  newG[0].inGroup = 0
	  G_list.pop(i)							  G_list.pop(i)
	return								return
#############################################################	#############################################################
def destroy_group(G_list, gal_id):				def destroy_group(G_list, gal_id):
  								  
  j = -1							  j = -1
  N_grp = len(G_list)						  N_grp = len(G_list)
  i = 0								  i = 0
  while i < N_grp:						  while i < N_grp:
    group = G_list[i]						    group = G_list[i]
    for galaxy in group[1:]:					    for galaxy in group[1:]:
      if galaxy.id == gal_id:					      if galaxy.id == gal_id:
	j = i								j = i
	i = N_grp+1							i = N_grp+1
	break								break
    i+=1							    i+=1
  								  
  if j>=0:							  if j>=0:
    for galaxy in G_list[j][1:]:				    for galaxy in G_list[j][1:]:
      galaxy.inGroup = 0					      galaxy.inGroup = 0
    								    
    G_list.pop(j)						    G_list.pop(j)
  								  
  return							  return
#############################################################	#############################################################
# removing a galaxy from a group				# removing a galaxy from a group
def removeFrom2(G_list, gr, gal_id):				def removeFrom2(G_list, gr, gal_id):

  id = None							  id = None
  for i in range(len(G_list)):					  for i in range(len(G_list)):
    group = G_list[i]						    group = G_list[i]
    if group[0].nest == gr:					    if group[0].nest == gr:
      id = i							      id = i
      break							      break
  								  
  if id!= None:							  if id!= None:
    group = G_list[id]						    group = G_list[id]
    for galaxy in group[1:]:					    for galaxy in group[1:]:
      if galaxy.id == gal_id:					      if galaxy.id == gal_id:
	  newG = []							  newG = []
	  for galaxy in group[1:]:					  for galaxy in group[1:]:
	    if galaxy.id == gal_id:					    if galaxy.id == gal_id:
		  galaxy.inGroup = 0						  galaxy.inGroup = 0
	    else:							    else:
		  newG.append(galaxy)						  newG.append(galaxy)
	  								  
	  if len(newG) > 1 : 						  if len(newG) > 1 : 
	    grp = LoosgalListJoint(newG, grID = 500000000)		    grp = LoosgalListJoint(newG, grID = 500000000)
	    group = NodeLeaf(grp)					    group = NodeLeaf(grp)
	    G_list.pop(i)						    G_list.pop(i)
	    G_list.append(group)					    G_list.append(group)
	  elif  len(newG) == 1 :					  elif  len(newG) == 1 :
	    newG[0].inGroup = 0						    newG[0].inGroup = 0
	    G_list.pop(i)						    G_list.pop(i)
	    								    
	  return							  return
#############################################################	#############################################################


#############################################################	#############################################################
# removing a galaxy from a group				# removing a galaxy from a group
def trimGroup_vel(G_list, gr, Vmin=-100000, Vmax=100000):	def trimGroup_vel(G_list, gr, Vmin=-100000, Vmax=100000):
  								  
  id = None							  id = None
  for i in range(len(G_list)):					  for i in range(len(G_list)):
     group = G_list[i]						     group = G_list[i]
     if group[0].nest == gr:					     if group[0].nest == gr:
       id = i							       id = i
       break							       break
  								  
  trimed_gal = []						  trimed_gal = []
  if id != None:						  if id != None:
     newGroup = []						     newGroup = []
     for p in range(1, len(G_list[id])):			     for p in range(1, len(G_list[id])):
        galaxy = G_list[id][p]					        galaxy = G_list[id][p]
        if galaxy.Vls < Vmax and galaxy.Vls > Vmin:		        if galaxy.Vls < Vmax and galaxy.Vls > Vmin:
	  newGroup.append(galaxy)					  newGroup.append(galaxy)
	else:								else:
	  galaxy.inGroup = 0						  galaxy.inGroup = 0
	  trimed_gal.append(galaxy)					  trimed_gal.append(galaxy)
 								 
     root = LoosgalListJoint(newGroup, grID = 700000000) 	     root = LoosgalListJoint(newGroup, grID = 700000000) 
     group = NodeLeaf(root)   					     group = NodeLeaf(root)   
     G_list.pop(id)						     G_list.pop(id)
     G_list.append(group)					     G_list.append(group)
     return trimed_gal						     return trimed_gal
  else:								  else:
    return None							    return None
#############################################################	#############################################################
# removing a galaxy from a group				# removing a galaxy from a group
# rad (radius) --> in degrees					# rad (radius) --> in degrees
def trimGroup_rad(G_list, gr, factor=1, rad=None):		def trimGroup_rad(G_list, gr, factor=1, rad=None):
  								  
  id = None							  id = None
  for i in range(len(G_list)):					  for i in range(len(G_list)):
     group = G_list[i]						     group = G_list[i]
     if group[0].nest == gr:					     if group[0].nest == gr:
       id = i							       id = i
       break							       break
     								     
  if id != None:						  if id != None:
     newGroup = []						     newGroup = []
     for p in range(1, len(G_list[id])):			     for p in range(1, len(G_list[id])):
        galaxy = G_list[id][p]					        galaxy = G_list[id][p]
        ang12 = (180./pi)*angle(G_list[id][0].gl, G_list[id][	        ang12 = (180./pi)*angle(G_list[id][0].gl, G_list[id][
        d = G_list[i][0].Vls / H0				        d = G_list[i][0].Vls / H0
        if d < 1: d = 1 					        if d < 1: d = 1 
        r = (180.*atan(G_list[i][0].R_2t2/(d))/pi)		        r = (180.*atan(G_list[i][0].R_2t2/(d))/pi)
        							        
        if rad != None:						        if rad != None:
	  radius = rad*pi/180.  # --> radius in radian			  radius = rad*pi/180.  # --> radius in radian
	else: 								else: 
	  radius = factor*G_list[id][0].R_2t2/d				  radius = factor*G_list[id][0].R_2t2/d
        							        
        if ang12 <= 1.05*(180./pi)*atan(radius):		        if ang12 <= 1.05*(180./pi)*atan(radius):
	  newGroup.append(galaxy)					  newGroup.append(galaxy)
	else:								else:
	  galaxy.inGroup = 0						  galaxy.inGroup = 0
 								 
     root = LoosgalListJoint(newGroup, grID = 700000000) 	     root = LoosgalListJoint(newGroup, grID = 700000000) 
     group = NodeLeaf(root)   					     group = NodeLeaf(root)   
     G_list.pop(id)						     G_list.pop(id)
     G_list.append(group)   					     G_list.append(group)   
        							        

#############################################################	#############################################################
# removing a galaxy from a group				# removing a galaxy from a group
def manualBlend(G_list, gr1, gr2):				def manualBlend(G_list, gr1, gr2):
  								  
   id1 = None; id2 = None					   id1 = None; id2 = None
   								   
   for i in range(len(G_list)):					   for i in range(len(G_list)):
     group = G_list[i]						     group = G_list[i]
     if group[0].nest == gr1:					     if group[0].nest == gr1:
       id1 = i							       id1 = i
       break							       break
   for j in range(len(G_list)):					   for j in range(len(G_list)):
     group = G_list[j]						     group = G_list[j]
     if group[0].nest == gr2:					     if group[0].nest == gr2:
       id2 = j							       id2 = j
       break   							       break   
   								   
   if id1 != None and id2 != None:				   if id1 != None and id2 != None:
     newGroup = []						     newGroup = []
     for p in range(1, len(G_list[id1])):			     for p in range(1, len(G_list[id1])):
	newGroup.append(G_list[id1][p])					newGroup.append(G_list[id1][p])
     for q in range(1, len(G_list[id2])):			     for q in range(1, len(G_list[id2])):
	newGroup.append(G_list[id2][q])   				newGroup.append(G_list[id2][q])   
     								     
     root = LoosgalListJoint(newGroup, grID = 600000000) 	     root = LoosgalListJoint(newGroup, grID = 600000000) 
     group = NodeLeaf(root)   					     group = NodeLeaf(root)   
     i = max(id1, id2); j = min(id1, id2)			     i = max(id1, id2); j = min(id1, id2)
     G_list.pop(i); G_list.pop(j)				     G_list.pop(i); G_list.pop(j)
     G_list.append(group)					     G_list.append(group)
   								   
    								    
#############################################################	#############################################################
# Note that dcf2 and ed would remain the same			# Note that dcf2 and ed would remain the same
# Bu their effect is discarded					# Bu their effect is discarded
#def distanceDiscard(G_list, id):				#def distanceDiscard(G_list, id):
  								  
  #for i in range(len(G_list)):					  #for i in range(len(G_list)):
    #group = G_list[i]						    #group = G_list[i]
    #for galaxy in group[1:]:					    #for galaxy in group[1:]:
      #if galaxy.id == id:					      #if galaxy.id == id:
        #newG = []						        #newG = []
        #for galaxy in group[1:]:				        #for galaxy in group[1:]:
           #if galaxy.id == id:					           #if galaxy.id == id:
        	#galaxy.sumDist = 0				        	#galaxy.sumDist = 0
        	#galaxy.sumError = 0 				        	#galaxy.sumError = 0 
           #newG.append(galaxy)					           #newG.append(galaxy)
        							        
        #grp = LoosgalListJoint(newG, grID = 500000000)		        #grp = LoosgalListJoint(newG, grID = 500000000)
        #group = NodeLeaf(grp)					        #group = NodeLeaf(grp)
        #G_list.pop(i)						        #G_list.pop(i)
        #G_list.append(group)					        #G_list.append(group)
	#return  							#return  

#############################################################	#############################################################
def addTo(G_list, galList, gr_id, gal_id):			def addTo(G_list, galList, gr_id, gal_id):
  								  
  removeFrom(G_list, gal_id)					  removeFrom(G_list, gal_id)
  								  
  id = None							  id = None
  B = False							  B = False
  for i in range(len(G_list)):					  for i in range(len(G_list)):
     group = G_list[i]						     group = G_list[i]
     if group[0].nest == gr_id:					     if group[0].nest == gr_id:
       id = i							       id = i
       break							       break
     else:							     else:
       for gal in group[1:]:					       for gal in group[1:]:
	 if gal.id ==  gr_id:						 if gal.id ==  gr_id:
	   id = i							   id = i
	   B = True							   B = True
	   break							   break
       if B: break   						       if B: break   
     								     
  idd = None							  idd = None
  for p in range(len(galList)):					  for p in range(len(galList)):
     galaxy = galList[p]					     galaxy = galList[p]
     if galaxy.id == gal_id:					     if galaxy.id == gal_id:
       idd = p							       idd = p
       break							       break
     								     
  								  
  if id != None and idd != None:				  if id != None and idd != None:
     newGroup = G_list[id][1:]					     newGroup = G_list[id][1:]
     newGroup.append(galList[p])				     newGroup.append(galList[p])
     galList[p].inGroup = 1					     galList[p].inGroup = 1

 								 
     root = LoosgalListJoint(newGroup, grID = 700000000, noSo	     root = LoosgalListJoint(newGroup, grID = 700000000, noSo
     group = NodeLeaf(root)   					     group = NodeLeaf(root)   
     G_list.pop(id)						     G_list.pop(id)
     G_list.append(group)       				     G_list.append(group)       
#############################################################	#############################################################
def formGroup(G_list, galList, galaxies):			def formGroup(G_list, galList, galaxies):
  								  
  newGroup=[]							  newGroup=[]
  for gal_id in galaxies:					  for gal_id in galaxies:
    brake = False						    brake = False
    for group in G_list:					    for group in G_list:
      for p in range(1, len(group)):				      for p in range(1, len(group)):
	if group[p].id == gal_id:					if group[p].id == gal_id:
	  brake = True							  brake = True
	  break								  break
      if brake: break						      if brake: break
    if brake == False:						    if brake == False:
      for galaxy in galList:					      for galaxy in galList:
        if galaxy.id == gal_id:					        if galaxy.id == gal_id:
	    if galaxy.inGroup == 1: 					    if galaxy.inGroup == 1: 
	      break	  						      break	  
	    else: 							    else: 
	      newGroup.append(galaxy)					      newGroup.append(galaxy)
	      galaxy.inGroup = 1					      galaxy.inGroup = 1
	      break							      break

  if len(newGroup) > 1:						  if len(newGroup) > 1:
    root = LoosgalListJoint(newGroup, grID = 700000000) 	    root = LoosgalListJoint(newGroup, grID = 700000000) 
    group = NodeLeaf(root)   					    group = NodeLeaf(root)   
    G_list.append(group)  					    G_list.append(group)  
  								  

#############################################################	#############################################################

def brent_modifications(G_list, galList):			def brent_modifications(G_list, galList):
  								  
  ### South							  ### South
  removeFrom(G_list, 12662)  # gal_id				  removeFrom(G_list, 12662)  # gal_id
  removeFrom(G_list, 2578)   # gal_id				  removeFrom(G_list, 2578)   # gal_id
  manualBlend(G_list, 62836, 62453) # group_id1, group_id2	  manualBlend(G_list, 62836, 62453) # group_id1, group_id2
  trimGroup_vel(G_list, 13255,  Vmax=1600) # gr_id1, remove t	  trimGroup_vel(G_list, 13255,  Vmax=1600) # gr_id1, remove t
  #distanceDiscard(G_list, 13368) # gal_id			  #distanceDiscard(G_list, 13368) # gal_id
  								  
  removeFrom2(G_list, 2557, 4126)				  removeFrom2(G_list, 2557, 4126)
  removeFrom2(G_list, 2789, 2578)				  removeFrom2(G_list, 2789, 2578)
  addTo(G_list, galList, 13505, 13470)				  addTo(G_list, galList, 13505, 13470)
  								  
  removeFrom2(G_list, 70090, 70069)				  removeFrom2(G_list, 70090, 70069)
  removeFrom2(G_list, 70090, 70184)				  removeFrom2(G_list, 70090, 70184)
  From_To(G_list, galList, 70080, 70184, 70069)			  From_To(G_list, galList, 70080, 70184, 70069)
  								  
  addTo(G_list, galList, 70069, 70306)				  addTo(G_list, galList, 70069, 70306)
  addTo(G_list, galList, 70069, 135466)				  addTo(G_list, galList, 70069, 135466)
    								    
  addTo(G_list, galList, 13434, 13154)				  addTo(G_list, galList, 13434, 13154)
  								  
  addTo(G_list, galList, 13419, 13342)				  addTo(G_list, galList, 13419, 13342)
  removeFrom(G_list, 13368)					  removeFrom(G_list, 13368)
  addTo(G_list, galList, 13419, 13304)				  addTo(G_list, galList, 13419, 13304)
  								  
  								  
  ### North							  ### North
  trimGroup_vel(G_list, 34695,  Vmin=400, Vmax=1000)		  trimGroup_vel(G_list, 34695,  Vmin=400, Vmax=1000)
  trimGroup_rad(G_list, 34695,  factor = 1.05)			  trimGroup_rad(G_list, 34695,  factor = 1.05)
  								  
  addTo(G_list, galList, 34426, 86673)				  addTo(G_list, galList, 34426, 86673)
    								    
  trimGroup_vel(G_list, 32256,  Vmax=800) # 			  trimGroup_vel(G_list, 32256,  Vmax=800) # 
  								  
  trimGroup_vel(G_list, 39600,  Vmax=700) # 			  trimGroup_vel(G_list, 39600,  Vmax=700) # 
  addTo(G_list, galList, 39241, 39344) # gr_id, gal_id		  addTo(G_list, galList, 39241, 39344) # gr_id, gal_id
  addTo(G_list, galList, 39241, 40228)				  addTo(G_list, galList, 39241, 40228)
  addTo(G_list, galList, 39600, 40537)				  addTo(G_list, galList, 39600, 40537)
  addTo(G_list, galList, 39600, 39191)				  addTo(G_list, galList, 39600, 39191)
  								  
  addTo(G_list, galList, 39600, 166131)				  addTo(G_list, galList, 39600, 166131)
  addTo(G_list, galList, 39600, 2832112)			  addTo(G_list, galList, 39600, 2832112)
  addTo(G_list, galList, 39600, 5057019)			  addTo(G_list, galList, 39600, 5057019)
  addTo(G_list, galList, 39600, 5057020)			  addTo(G_list, galList, 39600, 5057020)
  addTo(G_list, galList, 39600, 5061328)			  addTo(G_list, galList, 39600, 5061328)
  								  
  addTo(G_list, galList, 39600, 2832111)			  addTo(G_list, galList, 39600, 2832111)
  removeFrom2(G_list, 39600, 40665)				  removeFrom2(G_list, 39600, 40665)
  								  
  addTo(G_list, galList, 39241, 39237)				  addTo(G_list, galList, 39241, 39237)
  addTo(G_list, galList, 39241, 39864)				  addTo(G_list, galList, 39241, 39864)
  addTo(G_list, galList, 39241, 166129)				  addTo(G_list, galList, 39241, 166129)

  								  
  trimGroup_vel(G_list, 38440,  Vmin=700) # 			  trimGroup_vel(G_list, 38440,  Vmin=700) # 
  addTo(G_list, galList, 38440, 38068)				  addTo(G_list, galList, 38440, 38068)
  								  
  manualBlend(G_list, 28630, 28655)				  manualBlend(G_list, 28630, 28655)
  								  
  removeFrom(G_list, 24213) 					  removeFrom(G_list, 24213) 
  addTo(G_list, galList, 23324, 24213)				  addTo(G_list, galList, 23324, 24213)
  								  
  manualBlend(G_list, 21396, 21102)				  manualBlend(G_list, 21396, 21102)
  								  
  manualBlend(G_list, 46957, 45279)				  manualBlend(G_list, 46957, 45279)
  removeFrom2(G_list, 46957, 47847)				  removeFrom2(G_list, 46957, 47847)
  removeFrom2(G_list, 46957, 166152)				  removeFrom2(G_list, 46957, 166152)
  removeFrom2(G_list, 46957, 592761)				  removeFrom2(G_list, 46957, 592761)
  addTo(G_list, galList, 46957, 46680)				  addTo(G_list, galList, 46957, 46680)
  addTo(G_list, galList, 46957, 166158)				  addTo(G_list, galList, 46957, 166158)
  addTo(G_list, galList, 46957, 166167)				  addTo(G_list, galList, 46957, 166167)
  addTo(G_list, galList, 46957, 166172)				  addTo(G_list, galList, 46957, 166172)
  addTo(G_list, galList, 46957, 166175)				  addTo(G_list, galList, 46957, 166175)
  addTo(G_list, galList, 46957, 2815820)			  addTo(G_list, galList, 46957, 2815820)
  addTo(G_list, galList, 46957, 2815822)			  addTo(G_list, galList, 46957, 2815822)
  addTo(G_list, galList, 46957, 2815823)			  addTo(G_list, galList, 46957, 2815823)
  addTo(G_list, galList, 46957, 4689187)			  addTo(G_list, galList, 46957, 4689187)
  								  
  removeFrom2(G_list, 48082, 48334)				  removeFrom2(G_list, 48082, 48334)
  								  
  								  
  								  
  trimGroup_vel(G_list, 43495,  Vmax=500)			  trimGroup_vel(G_list, 43495,  Vmax=500)
  trimGroup_rad(G_list, 43495,  factor = 1.1)			  trimGroup_rad(G_list, 43495,  factor = 1.1)
  removeFrom2(G_list, 40973, 40904)				  removeFrom2(G_list, 40973, 40904)
  addTo(G_list, galList, 43495, 40973)				  addTo(G_list, galList, 43495, 40973)
  addTo(G_list, galList, 43495, 45314)				  addTo(G_list, galList, 43495, 45314)
  addTo(G_list, galList, 43495, 166146)				  addTo(G_list, galList, 43495, 166146)
  								  
  addTo(G_list, galList, 42575, 166140)				  addTo(G_list, galList, 42575, 166140)
  addTo(G_list, galList, 42575, 42045)				  addTo(G_list, galList, 42575, 42045)
  removeFrom(G_list, 41902)					  removeFrom(G_list, 41902)
  addTo(G_list, galList, 42575, 41902)				  addTo(G_list, galList, 42575, 41902)
  								  
  addTo(G_list, galList, 39225, 38881)				  addTo(G_list, galList, 39225, 38881)
  								  
  manualBlend(G_list, 39225, 39023)				  manualBlend(G_list, 39225, 39023)
  								  
  								  
  trimGroup_vel(G_list, 47404,  Vmin=400)			  trimGroup_vel(G_list, 47404,  Vmin=400)
  								  
  formGroup(G_list, galList, [46039, 45939, 45506, 46127, 505	  formGroup(G_list, galList, [46039, 45939, 45506, 46127, 505
  								  
  removeFrom2(G_list, 29265, 29033)				  removeFrom2(G_list, 29265, 29033)
  								  
  manualBlend(G_list, 13826, 15345)				  manualBlend(G_list, 13826, 15345)
  								  
  manualBlend(G_list, 25950, 26259)				  manualBlend(G_list, 25950, 26259)
  								  
  addTo(G_list, galList, 24930, 166096)				  addTo(G_list, galList, 24930, 166096)
  								  

  								  
  addTo(G_list, galList, 39724, 40692)				  addTo(G_list, galList, 39724, 40692)
  addTo(G_list, galList, 40692, 38742)				  addTo(G_list, galList, 40692, 38742)
  addTo(G_list, galList, 40692, 39690)				  addTo(G_list, galList, 40692, 39690)
  addTo(G_list, galList, 40692, 38598)				  addTo(G_list, galList, 40692, 38598)
  addTo(G_list, galList, 40692, 40495)				  addTo(G_list, galList, 40692, 40495)
  addTo(G_list, galList, 40692, 4326021)			  addTo(G_list, galList, 40692, 4326021)
  #addTo(G_list, galList, 40692, )				  #addTo(G_list, galList, 40692, )
  								  
  From_To(G_list, galList, 30083, 30445, 30493)			  From_To(G_list, galList, 30083, 30445, 30493)

  From_To(G_list, galList, 30744, 31166, 31029)			  From_To(G_list, galList, 30744, 31166, 31029)
  								  
  removeFrom2(G_list, 39225, 38685)				  removeFrom2(G_list, 39225, 38685)
  								  
  removeFrom2(G_list, 41789, 1181814)				  removeFrom2(G_list, 41789, 1181814)
  #distanceDiscard(G_list, 43775)				  #distanceDiscard(G_list, 43775)
  #distanceDiscard(G_list, 42447)				  #distanceDiscard(G_list, 42447)
  								  
  addTo(G_list, galList, 38440, 38356)				  addTo(G_list, galList, 38440, 38356)
  								  
  ### NGC784 Group (7671)					  ### NGC784 Group (7671)
  #addTo(G_list, galList, 7671, 166064)				  #addTo(G_list, galList, 7671, 166064)
  #addTo(G_list, galList, 7671, 6699)				  #addTo(G_list, galList, 7671, 6699)
  #addTo(G_list, galList, 7671, 8484)				  #addTo(G_list, galList, 7671, 8484)
  								  
  ### Cen A Group (46957)					  ### Cen A Group (46957)
  addTo(G_list, galList, 46957, 44110)				  addTo(G_list, galList, 46957, 44110)
  								  
  ### NGC1313 Group (12286)					  ### NGC1313 Group (12286)
  addTo(G_list, galList, 12286, 166073)				  addTo(G_list, galList, 12286, 166073)
  								  
  ### M83 Group (48082)						  ### M83 Group (48082)
  addTo(G_list, galList, 48082, 166170)				  addTo(G_list, galList, 48082, 166170)
  addTo(G_list, galList, 48082, 166176)				  addTo(G_list, galList, 48082, 166176)
  								  
  ### M81 Group (28630)						  ### M81 Group (28630)
  addTo(G_list, galList, 28630, 29231)				  addTo(G_list, galList, 28630, 29231)
  addTo(G_list, galList, 28630, 31286)				  addTo(G_list, galList, 28630, 31286)
  addTo(G_list, galList, 28630, 166101)				  addTo(G_list, galList, 28630, 166101)
  								  
  ### Hydra group - V_ls < 6000					  ### Hydra group - V_ls < 6000
  addTo(G_list, galList, 31478, 31599)				  addTo(G_list, galList, 31478, 31599)
  addTo(G_list, galList, 31478, 31951)				  addTo(G_list, galList, 31478, 31951)

  ### Centaurus group - V_ls < 6000  				  ### Centaurus group - V_ls < 6000  
  addTo(G_list, galList, 43296, 43466)				  addTo(G_list, galList, 43296, 43466)
  							      <
  manual_group(G_list, galList, [39225, 39023, 38881, 39145]) <
#############################################################	#############################################################

def From_To(G_list, galList, gr1, gr2, gal_id):			def From_To(G_list, galList, gr1, gr2, gal_id):
  								  
  removeFrom2(G_list, gr1, gal_id)				  removeFrom2(G_list, gr1, gal_id)
  addTo(G_list, galList, gr2, gal_id)				  addTo(G_list, galList, gr2, gal_id)
  								  

def isAllowed(group_id, galaxy_id):				def isAllowed(group_id, galaxy_id):
  								  
  allowed = True						  allowed = True
  								  
  ## not allowed list						  ## not allowed list
  if group_id == 34695 and galaxy_id in [3471336, 4098734, 28	  if group_id == 34695 and galaxy_id in [3471336, 4098734, 28
    allowed = False						    allowed = False
  								  
  return allowed						  return allowed
  								  

#############################################################	#############################################################
def manual_group(G_list, galList, id_list, verbose=False):  #	def manual_group(G_list, galList, id_list, verbose=False):  #
  								  
  								  
    my_group = []						    my_group = []
    for galaxy in galList: 					    for galaxy in galList: 
      id = galaxy.id						      id = galaxy.id
      if id in id_list: 					      if id in id_list: 
	removeFrom(G_list, id)						removeFrom(G_list, id)
	my_group.append(galaxy)						my_group.append(galaxy)
	galaxy.inGroup = 1						galaxy.inGroup = 1
	galaxy.logK = m_logK(galaxy.Ks, galaxy.Vls, d=galaxy.		galaxy.logK = m_logK(galaxy.Ks, galaxy.Vls, d=galaxy.
							      |	    root = LoosgalListJoint(my_group, grID = 900000000)
    if len(my_group)>0: 				      |	    new_Group = NodeLeaf(root)
      root = LoosgalListJoint(my_group, grID = 900000000)     |	    G_list.append(new_Group)  
      new_Group = NodeLeaf(root)			      <
      G_list.append(new_Group)  			      <
    else: return None					      <
    								    
    if verbose:							    if verbose:
      for all in new_Group:					      for all in new_Group:
	print all.id							print all.id
      								      
							      |	    
							      >	    
    return my_group						    return my_group
  								  

#############################################################	#############################################################

if __name__ == '__main__':					if __name__ == '__main__':
  								  
  								  


  								  
  R_max=50. # R_max=50.						  R_max=50. # R_max=50.
  cluster = ''							  cluster = ''

  								  
  if len(sys.argv) < 2:						  if len(sys.argv) < 2:
    print "\nEnter the cluster name as the input ..." 		    print "\nEnter the cluster name as the input ..." 
    print "\nexample: \n > python " + str(sys.argv[0])+ " man	    print "\nexample: \n > python " + str(sys.argv[0])+ " man
    print "\nPossible options: all, manual, local" 		    print "\nPossible options: all, manual, local" 
    print "Use north/south for whole sky.\n"			    print "Use north/south for whole sky.\n"
    sys.exit(1)							    sys.exit(1)
  								  
  cluster = str(sys.argv[1])					  cluster = str(sys.argv[1])
  								  
  								  
  inFile = 'AllSky.v23.cf3.csv'					  inFile = 'AllSky.v23.cf3.csv'
  								  
  								  
  table = np.genfromtxt(inFile , delimiter=',', filling_value	  table = np.genfromtxt(inFile , delimiter=',', filling_value
  								  
  								  


  ignore_list = [40001, 42447, 40240, 40851, 3085, 70090]   #	  ignore_list = [40001, 42447, 40240, 40851, 3085, 70090]   #
  Vh_zero = [29231,  31286,  40747,  46680, 166101, 166146, 1	  Vh_zero = [29231,  31286,  40747,  46680, 166101, 166146, 1
  								  
 								 
  ignore_list = np.asarray(ignore_list)				  ignore_list = np.asarray(ignore_list)
  Vh_zero = np.asarray(Vh_zero)					  Vh_zero = np.asarray(Vh_zero)
  ignore_list = np.concatenate((ignore_list, Vh_zero))		  ignore_list = np.concatenate((ignore_list, Vh_zero))
  								  
  # added May 3, 2016						  # added May 3, 2016
  Vh_zero = [166167, 2815820, 2815822, 2815823, 4689187]	  Vh_zero = [166167, 2815820, 2815822, 2815823, 4689187]
  Vh_zero = np.asarray(Vh_zero)					  Vh_zero = np.asarray(Vh_zero)
  ignore_list = np.concatenate((ignore_list, Vh_zero))		  ignore_list = np.concatenate((ignore_list, Vh_zero))
  								  


  # <> ---------- Local groups					  # <> ---------- Local groups
  MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', f	  MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', f
  id_MW = MW_g['PGC']						  id_MW = MW_g['PGC']
  ignore_list = np.concatenate((ignore_list, id_MW))		  ignore_list = np.concatenate((ignore_list, id_MW))
  								  
  M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=',',	  M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=',',
  id_M31 = M31_g['PGC']						  id_M31 = M31_g['PGC']
  ignore_list = np.concatenate((ignore_list, id_M31))		  ignore_list = np.concatenate((ignore_list, id_M31))
  # </> ---------- Local groups					  # </> ---------- Local groups
  								  

  								  
  if cluster == 'all' :						  if cluster == 'all' :
     galList = readgalList(table, ignore_list=ignore_list, Vl	     galList = readgalList(table, ignore_list=ignore_list, Vl
  elif cluster == 'manual' :					  elif cluster == 'manual' :
     alpha =  333.4714					      |	     alpha = 156.2981
     delta = -7.0663					      |	     delta = -11.6740
     step = 10							     step = 10
     galList = readgalList(table, skyPatch='manual', sgl_min=	     galList = readgalList(table, skyPatch='manual', sgl_min=
  elif cluster == 'local':					  elif cluster == 'local':
     galList = readgalList(table, skyPatch='local', Vls_max=2	     galList = readgalList(table, skyPatch='local', Vls_max=2
  else:								  else:
    sys.exit(1)							    sys.exit(1)

  								  
  ### Virgo: Ignoring galaxies inside 6.8 deg			  ### Virgo: Ignoring galaxies inside 6.8 deg
  Virgo_ignore = []						  Virgo_ignore = []
  for galaxy in galList:					  for galaxy in galList:
    ang12 = (180./pi)*angle(galaxy.sgl, galaxy.sgb, sglV, sgb	    ang12 = (180./pi)*angle(galaxy.sgl, galaxy.sgb, sglV, sgb
    if ang12<= 6.8:						    if ang12<= 6.8:
      Virgo_ignore.append(galaxy.id)				      Virgo_ignore.append(galaxy.id)
  								  
  if len(Virgo_ignore)>0:					  if len(Virgo_ignore)>0:
    Virgo_ignore = np.asarray(Virgo_ignore)			    Virgo_ignore = np.asarray(Virgo_ignore)
    ignore_list = np.concatenate((ignore_list, Virgo_ignore))	    ignore_list = np.concatenate((ignore_list, Virgo_ignore))
  								  
  Virgo_W1 = [40375, 40122, 40886, 40516, 40033, 40119,41050]	  Virgo_W1 = [40375, 40122, 40886, 40516, 40033, 40119,41050]
  Virgo_W1 = np.asarray(Virgo_W1)				  Virgo_W1 = np.asarray(Virgo_W1)
  								  
  Virgo_W2 = [42741, 42619, 43001]				  Virgo_W2 = [42741, 42619, 43001]
  Virgo_W2 = np.asarray(Virgo_W2)				  Virgo_W2 = np.asarray(Virgo_W2)
  								  
  Virgo_M   = [38890, 38916, 39002, 39025, 39040, 39152, 3939	  Virgo_M   = [38890, 38916, 39002, 39025, 39040, 39152, 3939
  Virgo_M   = np.asarray(Virgo_M)				  Virgo_M   = np.asarray(Virgo_M)
  								  

  Virgo_W_g = np.genfromtxt('brent_virgoW.csv' , delimiter=',	  Virgo_W_g = np.genfromtxt('brent_virgoW.csv' , delimiter=',
  id_Virgo_w = Virgo_W_g['PGC'][0:41]				  id_Virgo_w = Virgo_W_g['PGC'][0:41]
  ignore_list = np.concatenate((ignore_list, id_Virgo_w))	  ignore_list = np.concatenate((ignore_list, id_Virgo_w))


  Virgo_foreback = [40045, 42081, 43072, 44491, 43601, 43254,	  Virgo_foreback = [40045, 42081, 43072, 44491, 43601, 43254,
  Virgo_foreback = np.asarray(Virgo_foreback)			  Virgo_foreback = np.asarray(Virgo_foreback)
  Virgo_foreback = np.concatenate((Virgo_foreback, Virgo_W_g[	  Virgo_foreback = np.concatenate((Virgo_foreback, Virgo_W_g[
  #ignore_list = np.concatenate((ignore_list, Virgo_foreback)	  #ignore_list = np.concatenate((ignore_list, Virgo_foreback)
  								  
  virgo_g1 = [39620,40004,40087,40240,40339,40469,40494,40801	  virgo_g1 = [39620,40004,40087,40240,40339,40469,40494,40801
  virgo_g1 = np.asarray(virgo_g1)				  virgo_g1 = np.asarray(virgo_g1)
  ignore_list = np.concatenate((ignore_list, virgo_g1))		  ignore_list = np.concatenate((ignore_list, virgo_g1))
  								  
  virgo_g2 = [39483,39537,39646,39809,39878,40109,40229,40321	  virgo_g2 = [39483,39537,39646,39809,39878,40109,40229,40321
  virgo_g2 = np.asarray(virgo_g2)				  virgo_g2 = np.asarray(virgo_g2)
  ignore_list = np.concatenate((ignore_list, virgo_g2))		  ignore_list = np.concatenate((ignore_list, virgo_g2))
  								  
# ===========================================================	# ===========================================================
  #id_all   = table['pgc']					  #id_all   = table['pgc']
  #Ks_all   = table['Ks'] 					  #Ks_all   = table['Ks'] 
  								  
  								  
  #MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', 	  #MW_g = np.genfromtxt('brent_MWgroup.csv' , delimiter=',', 
  #id = MW_g['PGC']						  #id = MW_g['PGC']
  #Vmag = MW_g['Vmag']						  #Vmag = MW_g['Vmag']
  #Ks0=[]							  #Ks0=[]
  #Vm0=[]							  #Vm0=[]
  #for i in range(len(id)):					  #for i in range(len(id)):
    #for j in range(len(id_all)):				    #for j in range(len(id_all)):
      #if id_all[j] == id[i] and Vmag[i]!=Ks_all[j] and Vmag[	      #if id_all[j] == id[i] and Vmag[i]!=Ks_all[j] and Vmag[
	#Ks0.append(Vmag[i]-Ks_all[j])					#Ks0.append(Vmag[i]-Ks_all[j])
	#Vm0.append(Vmag[i])						#Vm0.append(Vmag[i])
	#break								#break
    								    
  #plt.plot(Vm0, Ks0, '+')					  #plt.plot(Vm0, Ks0, '+')
  								  
  								  
  								  
  #M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=','	  #M31_g = np.genfromtxt('brent_M31group.csv' , delimiter=','
  #id = M31_g['PGC']						  #id = M31_g['PGC']
  #Vmag = M31_g['Vmag']						  #Vmag = M31_g['Vmag']
  #Ks=[]							  #Ks=[]
  #Vm=[]							  #Vm=[]
  #for i in range(len(id)):					  #for i in range(len(id)):
    #for j in range(len(id_all)):				    #for j in range(len(id_all)):
      #if id_all[j] == id[i] and Vmag[i]!=Ks_all[j] and Vmag[	      #if id_all[j] == id[i] and Vmag[i]!=Ks_all[j] and Vmag[
	#Ks.append(Vmag[i]-Ks_all[j])					#Ks.append(Vmag[i]-Ks_all[j])
	#Vm.append(Vmag[i])						#Vm.append(Vmag[i])
	#break								#break
    								    
  #plt.plot(Vm, Ks, '*', color='red')  				  #plt.plot(Vm, Ks, '*', color='red')  

  #Ks = np.asarray(Ks)						  #Ks = np.asarray(Ks)
  #Vm = np.asarray(Vm)						  #Vm = np.asarray(Vm)
  #Ks0 = np.asarray(Ks0)					  #Ks0 = np.asarray(Ks0)
  #Vm0 = np.asarray(Vm0)  					  #Vm0 = np.asarray(Vm0)  
  #Ks = np.concatenate((Ks0, Ks))				  #Ks = np.concatenate((Ks0, Ks))
  #Vm = np.concatenate((Vm0, Vm))				  #Vm = np.concatenate((Vm0, Vm))
  #sig = np.std(Ks)						  #sig = np.std(Ks)
  #med = np.median(Ks)						  #med = np.median(Ks)
  #print med, sig						  #print med, sig
  #flag = np.zeros(len(Ks))					  #flag = np.zeros(len(Ks))
  #flag[np.where(Ks>med+1.5*sig)] = 1.				  #flag[np.where(Ks>med+1.5*sig)] = 1.
  #flag[np.where(Ks<med-1.5*sig)] = 1.				  #flag[np.where(Ks<med-1.5*sig)] = 1.
  #ind = np.where(flag == 0)					  #ind = np.where(flag == 0)
  #Ks = Ks[ind]							  #Ks = Ks[ind]
  #Vm = Vm[ind]							  #Vm = Vm[ind]
  #med = np.median(Ks)						  #med = np.median(Ks)
  #plt.plot([10, 20],[med, med], '--')				  #plt.plot([10, 20],[med, med], '--')
  #print med							  #print med
  								  
  								  
  #plt.xlabel('Vmag')						  #plt.xlabel('Vmag')
  #plt.ylabel('Vmag-Ks')					  #plt.ylabel('Vmag-Ks')
  #plt.show()							  #plt.show()
  #sys.exit(1)							  #sys.exit(1)
# ===========================================================	# ===========================================================
  								  

  for iter in range(4):						  for iter in range(4):
   								   
        for gal in galList:					        for gal in galList:
	  if gal.id in ignore_list:					  if gal.id in ignore_list:
	    gal.inGroup = -1						    gal.inGroup = -1
	    #print gal.id						    #print gal.id

        # these are then added manually				        # these are then added manually
        # Their Vh=0 but the belong to gr39600			        # Their Vh=0 but the belong to gr39600
        for gal in galList:					        for gal in galList:
	  if gal.id in [166131, 2832112, 5057019, 5057020, 50		  if gal.id in [166131, 2832112, 5057019, 5057020, 50
	    gal.inGroup = -1						    gal.inGroup = -1
    								    
        print "working on grouping ... iter:", iter		        print "working on grouping ... iter:", iter
	##  Finding groups based on the criteria			##  Finding groups based on the criteria
        G_list = []						        G_list = []
        							        

        							        
        							        
        for qp in range(0,3):					        for qp in range(0,3):
	   G_list += IndividualGal(galList)				   G_list += IndividualGal(galList)
	   G_list += find_pair_dist(galList)				   G_list += find_pair_dist(galList)
									
	## fixed							## fixed
        G_list = addGalGroup(G_list, galList)			        G_list = addGalGroup(G_list, galList)
        							        
        group_moderator(G_list)					        group_moderator(G_list)
        							        
        inject_list=[40001, 40240, 40851, 42447, 70090]   	        inject_list=[40001, 40240, 40851, 42447, 70090]   
        for gal in galList:					        for gal in galList:
	  if gal.id in inject_list:					  if gal.id in inject_list:
	    gal.inGroup = 0						    gal.inGroup = 0
        							        
        							        
        if True:						        if True:
          for pq in range(0,10): 				          for pq in range(0,10): 
	    G_list = IndividualGal_modifier(G_list, galList)		    G_list = IndividualGal_modifier(G_list, galList)
	    G_list += find_pair_dist(galList)				    G_list += find_pair_dist(galList)
	     								     
	  for pq in range(0,10): 					  for pq in range(0,10): 
	    G_list = addGalGroup(G_list, galList)			    G_list = addGalGroup(G_list, galList)
	     								     
          for qp in range(0,10):				          for qp in range(0,10):
            mergGroup(G_list)					            mergGroup(G_list)
          							          
          group_moderator(G_list)				          group_moderator(G_list)
          mergGroup(G_list, restrict=True)			          mergGroup(G_list, restrict=True)
           							           
	  for pq in range(0,10): 					  for pq in range(0,10): 
	    G_list = addGalGroup(G_list, galList)  			    G_list = addGalGroup(G_list, galList)  
	     								     
          for qp in range(0,2):					          for qp in range(0,2):
            mergGroup(G_list)					            mergGroup(G_list)
             							             
          mergGroup(G_list, restrict=True)    			          mergGroup(G_list, restrict=True)    
           							           
          group_moderator(G_list)				          group_moderator(G_list)
          							          


#############################################################	#############################################################

	##Error: galList is modified ?					##Error: galList is modified ?
	for galaxy in galList:						for galaxy in galList:
	  if galaxy.inGroup == 0:					  if galaxy.inGroup == 0:
	    galaxy.mDist = 0						    galaxy.mDist = 0
	    galaxy.setMeanDist(galaxy.dcf2, galaxy.ed, 0)		    galaxy.setMeanDist(galaxy.dcf2, galaxy.ed, 0)
	    								    

	brent_modifications(G_list, galList)  				brent_modifications(G_list, galList)  
	group_moderator(G_list)						group_moderator(G_list)
									
	G_list += IndividualGal(galList, pairs=True)			G_list += IndividualGal(galList, pairs=True)
	G_list = addGalGroup(G_list, galList)				G_list = addGalGroup(G_list, galList)
	brent_modifications(G_list, galList) 				brent_modifications(G_list, galList) 
	group_moderator(G_list)						group_moderator(G_list)
	G_list = addGalGroup(G_list, galList)				G_list = addGalGroup(G_list, galList)
									
	for qp in range(0,2):						for qp in range(0,2):
	  mergGroup(G_list, restrict=True)    				  mergGroup(G_list, restrict=True)    
        group_moderator(G_list)					        group_moderator(G_list)
        							        
        brent_modifications(G_list, galList) 			        brent_modifications(G_list, galList) 





        if True: # cluster=='manual':				        if True: # cluster=='manual':

	    destroy_group(G_list, 40001)				    destroy_group(G_list, 40001)
	    								    
	    for galaxy in galList: 					    for galaxy in galList: 
	      if galaxy.inGroup == 1 and galaxy.sgl>106 and g		      if galaxy.inGroup == 1 and galaxy.sgl>106 and g
		destroy_group(G_list, galaxy.id)				destroy_group(G_list, galaxy.id)
	    								    
	        							        
	    								    
	    ##manual_group(G_list, galList, virgo_g1)			    ##manual_group(G_list, galList, virgo_g1)
	    g_lst = []							    g_lst = []
	    for galaxy in galList:					    for galaxy in galList:
	      if galaxy.id in virgo_g1:					      if galaxy.id in virgo_g1:
		galaxy.inGroup = 0						galaxy.inGroup = 0
		g_lst.append(galaxy)						g_lst.append(galaxy)
		galaxy.setMeanDist(16, 0.1, GRP_vel = galaxy.			galaxy.setMeanDist(16, 0.1, GRP_vel = galaxy.
	    								    
	    new_grp = IndividualGal(g_lst)				    new_grp = IndividualGal(g_lst)
	    G_list += new_grp						    G_list += new_grp
	    meanDist    = new_grp[0][0].mDist				    meanDist    = new_grp[0][0].mDist
	    meanDistErr = new_grp[0][0].mDistErr			    meanDistErr = new_grp[0][0].mDistErr
	    for galaxy in g_lst:					    for galaxy in g_lst:
	      if galaxy.inGroup == 0: 					      if galaxy.inGroup == 0: 
		#galaxy.inGroup=-1						#galaxy.inGroup=-1
		galaxy.setMeanDist(meanDist, meanDistErr, GRP			galaxy.setMeanDist(meanDist, meanDistErr, GRP
	      								      
	    g_lst = []							    g_lst = []
	    for galaxy in galList:					    for galaxy in galList:
	      if galaxy.id in virgo_g2:					      if galaxy.id in virgo_g2:
		galaxy.inGroup = 0						galaxy.inGroup = 0
		g_lst.append(galaxy)						g_lst.append(galaxy)
		galaxy.setMeanDist(16, 0.1, GRP_vel = galaxy.			galaxy.setMeanDist(16, 0.1, GRP_vel = galaxy.
	    								    
	    new_grp = IndividualGal(g_lst)				    new_grp = IndividualGal(g_lst)
	    G_list += new_grp						    G_list += new_grp
	    meanDist    = new_grp[0][0].mDist				    meanDist    = new_grp[0][0].mDist
	    meanDistErr = new_grp[0][0].mDistErr			    meanDistErr = new_grp[0][0].mDistErr
	    for galaxy in g_lst:					    for galaxy in g_lst:
	      if galaxy.inGroup == 0: 					      if galaxy.inGroup == 0: 
		#galaxy.inGroup=-1						#galaxy.inGroup=-1
		galaxy.setMeanDist(meanDist, meanDistErr, GRP			galaxy.setMeanDist(meanDist, meanDistErr, GRP
	    								    
	    								    
	    addTo(G_list, galList, 39483, 39537)			    addTo(G_list, galList, 39483, 39537)
	    addTo(G_list, galList, 39483, 39620)			    addTo(G_list, galList, 39483, 39620)
	    addTo(G_list, galList, 39483, 39809)			    addTo(G_list, galList, 39483, 39809)
	    #addTo(G_list, galList, 39483, 3394197)			    #addTo(G_list, galList, 39483, 3394197)
	    								    
	    G_list = addGalGroup(G_list, galList)			    G_list = addGalGroup(G_list, galList)
	    								    
	    								    
	    manual_group(G_list, galList, Virgo_W1)			    manual_group(G_list, galList, Virgo_W1)
	    manual_group(G_list, galList, Virgo_W2)			    manual_group(G_list, galList, Virgo_W2)
	    manual_group(G_list, galList, Virgo_M)			    manual_group(G_list, galList, Virgo_M)
	    								    
									
										
	    #print id_Virgo_w						    #print id_Virgo_w
	    VirW_grp = manual_group(G_list, galList, id_Virgo		    VirW_grp = manual_group(G_list, galList, id_Virgo
	    #print "Virgo W length: " ,len(VirW_grp)			    #print "Virgo W length: " ,len(VirW_grp)
	    								    

	    								    
	    								    
            ### VVVVVIIIIIRRRRGGGGOOOOO				            ### VVVVVIIIIIRRRRGGGGOOOOO
            Virgo_grp = []					            Virgo_grp = []
            for galaxy in galList: 				            for galaxy in galList: 
	      if galaxy.id in Virgo_ignore and galaxy.inGroup		      if galaxy.id in Virgo_ignore and galaxy.inGroup
		Virgo_grp.append(galaxy.id)					Virgo_grp.append(galaxy.id)
            if len(Virgo_grp)>0:				            if len(Virgo_grp)>0:
	      print "Working on generating the Virgo cluster 		      print "Working on generating the Virgo cluster 
	      manual_group(G_list, galList, Virgo_grp)			      manual_group(G_list, galList, Virgo_grp)
	    								    
	    								    
            ### Assumption: any galaxy without a known distan	            ### Assumption: any galaxy without a known distan
            ### and has V_ls > 1200 km/s is in 39659 (Virgo W	            ### and has V_ls > 1200 km/s is in 39659 (Virgo W
	    for galaxy in galList:					    for galaxy in galList:
	      if (180./pi)*angle(galaxy.sgl, galaxy.sgb, 108.		      if (180./pi)*angle(galaxy.sgl, galaxy.sgb, 108.
		addTo(G_list, galList, 39659, galaxy.id)			addTo(G_list, galList, 39659, galaxy.id)
              # Virgo negative velocities			              # Virgo negative velocities
	      if (180./pi)*angle(galaxy.sgl, galaxy.sgb, sglV		      if (180./pi)*angle(galaxy.sgl, galaxy.sgb, sglV
		addTo(G_list, galList, 41220, galaxy.id)			addTo(G_list, galList, 41220, galaxy.id)
              							              

							      >			 
							      >		     
							      >			
							      >			
	    G_list = addGalGroup(G_list, galList)			    G_list = addGalGroup(G_list, galList)
	    G_list += IndividualGal(galList)				    G_list += IndividualGal(galList)
	    removeFrom2(G_list, 39659, 40001)				    removeFrom2(G_list, 39659, 40001)
	    trimed_virgoW = trimGroup_vel(G_list, 39659,  Vmi		    trimed_virgoW = trimGroup_vel(G_list, 39659,  Vmi
	    								    
	  								  
	    manualBlend(G_list, 40001, 40240)				    manualBlend(G_list, 40001, 40240)
	    addTo(G_list, galList, 40001, 40851)			    addTo(G_list, galList, 40001, 40851)
	    								    
	    								    
	    Grp40001 = []						    Grp40001 = []
	    for i in range(len(G_list)):				    for i in range(len(G_list)):
	      if G_list[i][0].nest == 40001:				      if G_list[i][0].nest == 40001:
		print "group 40001 exists ..."					print "group 40001 exists ..."
		Grp40001.append(G_list[i])					Grp40001.append(G_list[i])
		G_list.pop(i)							G_list.pop(i)
		break								break
	      								      
      								      
	    if len(Grp40001) > 0 and trimed_virgoW!=None: 		    if len(Grp40001) > 0 and trimed_virgoW!=None: 
	       print "adding to 40001 group ..."			       print "adding to 40001 group ..."
	       Grp40001 = addGalGroup(Grp40001, trimed_virgoW		       Grp40001 = addGalGroup(Grp40001, trimed_virgoW
	       G_list.append(Grp40001[0])				       G_list.append(Grp40001[0])
	    								    
	    								    
	    brent_modifications(G_list, galList)			    brent_modifications(G_list, galList)
    								    
   								   
	    for id in Virgo_foreback:					    for id in Virgo_foreback:
	      removeFrom(G_list, id)					      removeFrom(G_list, id)

        ## Very last stage ... No group formation after this 	        ## Very last stage ... No group formation after this 
        ### Hard coding groups					        ### Hard coding groups
        if cluster!='manual':					        if cluster!='manual':
	    manual_group(G_list, galList, id_MW)			    manual_group(G_list, galList, id_MW)
	    manual_group(G_list, galList, id_M31) 			    manual_group(G_list, galList, id_M31) 
 								 
							      |			
							      <
							      <
        group_moderator(G_list, use_dist=True)			        group_moderator(G_list, use_dist=True)
        galList_moderator(galList)				        galList_moderator(galList)


        groupWrite(cluster+'.iter.'+str(int(iter))+'.v37.grou |	        groupWrite(cluster+'.iter.'+str(int(iter))+'.v36.grou
        #sys.exit(1)						        #sys.exit(1)

	for galaxy in galList: 						for galaxy in galList: 
	  if galaxy.inGroup > 0 : 					  if galaxy.inGroup > 0 : 
	     galaxy.inGroup = 0						     galaxy.inGroup = 0
        							        

  								  
