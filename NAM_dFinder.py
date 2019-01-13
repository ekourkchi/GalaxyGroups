

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

import spherematch as spherematch

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


def Vgsrm2Vls(Vgsrm, l, b):
  
  l*=pi/180
  b*=pi/180
  Vls = Vgsrm - 37.*cos(l)*cos(b) + 66.*sin(l)*cos(b) - 15.*sin(b)
  return Vls


# **************************************
def nest(nest_file, nam_file):

    NESTS = np.genfromtxt(nest_file , delimiter=',', filling_values="-100000", names=True, dtype=None )
    nest_nest = NESTS['Nest']
    nest_sgl  = NESTS['sgl']
    nest_sgb  = NESTS['sgb']
    nest_gl   = NESTS['l']
    nest_gb   = NESTS['b']
    nest_Vls  = NESTS['Vls']
    nest_pgc  = NESTS['PGC1']
    nest_d    = NESTS['d']

    indices   = np.argsort(nest_nest)
    nest_nest = nest_nest[indices]
    nest_sgl  = nest_sgl[indices]
    nest_sgb  = nest_sgb[indices]
    nest_gl   = nest_gl[indices]
    nest_gb   = nest_gb[indices]
    nest_Vls  = nest_Vls[indices]
    nest_pgc  = nest_pgc[indices]
    nest_d    = nest_d[indices]

    NAM = np.genfromtxt(nam_file , delimiter=',', filling_values="-100000", names=True, dtype=None )
    nam_name  = NAM['name']
    nam_Vgsrm = NAM['cz_cat']
    nam_dist  = NAM['d']

    indices   = np.argsort(nam_name)
    nam_name  = nam_name[indices]
    nam_Vgsrm = nam_Vgsrm[indices]
    nam_dist  = nam_dist[indices]

    M = len(nest_nest)
    N = len(nam_name)

    i = 0 
    j = 0 
    p = 0 

    nest_dist    = np.zeros((M,), dtype=float)
    nest_flag    = np.zeros((M,), dtype=np.int)
    nest_dsource = np.zeros((M,), dtype=np.int)

    nam_sgl  = np.zeros((N,), dtype=float)
    nam_sgb  = np.zeros((N,), dtype=float)
    nam_gl   = np.zeros((N,), dtype=float)
    nam_gb   = np.zeros((N,), dtype=float)
    nam_Vls  = np.zeros((N,), dtype=float)
    nam_Vls_  = np.zeros((N,), dtype=float)

    while i < N:     # NAM -loop
      while j < M:   # NEST-loop
	
	if nam_name[i] == nest_nest[j]:
	    nam_sgl[i] = nest_sgl[j]
	    nam_sgb[i] = nest_sgb[j]
	    nam_gl[i]  = nest_gl[j]
	    nam_gb[i]  = nest_gb[j]
	    nam_Vls[i] = Vgsrm2Vls(nam_Vgsrm[i], nam_gl[i], nam_gb[i])
	    nam_Vls_[i] = nest_Vls[j]
	    
	    nest_dist[j]    = nam_dist[i]
	    nest_flag[j]    = 0
	    nest_dsource[j] = nam_name[i]
	    
	    #print nest_d[j],nam_dist[i]
	    
	    j+=1
	    p+=1
	    break
	j+=1
      
      i+=1

    if p!= N: print"Warning: Not found NEST ID for all NAM objects ..."
    #print p, N


    ################################################################################
    ################################################################################
    ###########################11111111111##########################################
    ################################################################################


    flag_indices = np.where(nest_flag==0)[0]   # where there is not previous match
    nest_nest_flag = nest_nest[flag_indices]
    nest_sgl_flag  = nest_sgl[flag_indices]
    nest_sgb_flag  = nest_sgb[flag_indices]


    (m1, m2, d)  = spherematch.match_radec(nest_sgl_flag, nest_sgb_flag, nam_sgl, nam_sgb, 5., notself=False, nearest=False)

    del_v  = np.zeros((len(d),), dtype=float)
    for i in range(len(d)):
      v1 =  nest_Vls[flag_indices[m1[i]]] 
      v2 =  nam_Vls_[m2[i]]
      del_v[i] = abs(v2-v1)/v1


    indices = np.argsort(m1)
    m1 = m1[indices]
    m2 = m2[indices]
    del_v = del_v[indices]

    m1_old = -1

    del_v_list_f = []
    m1_list_f = []
    m2_list_f = []
	
    del_v_list = []
    m1_list = []
    m2_list = []
    for i in range(len(m1)):
      
      if m1[i] != m1_old or i==len(m1)-1:
	
	if len(del_v_list) != 0:
	  indices = np.argsort(del_v_list)
	  index = indices[0]
	
	  del_v_list_f.append(del_v_list[index])
	  m1_list_f.append(m1_list[index])
	  m2_list_f.append(m2_list[index])
	
	del_v_list = []
	m1_list = []
	m2_list = []
	
	if del_v[i] < 0.3:
	  del_v_list.append(del_v[i])
	  m1_list.append(m1[i])
	  m2_list.append(m2[i])
	
	m1_old = m1[i]
	
      else:
	
	if del_v[i] < 0.3:
	  del_v_list.append(del_v[i])
	  m1_list.append(m1[i])
	  m2_list.append(m2[i])    


    #print "1: ", len(m1_list_f), " --------- "



    for p in range(len(m1_list_f)):
      
      j = flag_indices[m1_list_f[p]]
      i = m2_list_f[p]
      
      v1 = nest_Vls[j]
      v2 = nam_Vls_[i]  
      
      nest_dist[j]    = (v1/v2)*nam_dist[i]
      nest_flag[j]    = 1
      nest_dsource[j] = nam_name[i]
      

      
      #print v1, v2, abs(v2-v1)/v1  
    ################################################################################
    ################################################################################
    ################################################################################
    ###############################2222222222222####################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################


    flag_indices = np.where(nest_flag==0)[0]   # where there is not previous match
    nest_nest_flag = nest_nest[flag_indices]
    nest_sgl_flag  = nest_sgl[flag_indices]
    nest_sgb_flag  = nest_sgb[flag_indices]


    (m1, m2, d)  = spherematch.match_radec(nest_sgl_flag, nest_sgb_flag, nam_sgl, nam_sgb, 7., notself=False, nearest=False)

    del_v  = np.zeros((len(d),), dtype=float)
    for i in range(len(d)):
      v1 =  nest_Vls[flag_indices[m1[i]]] 
      v2 =  nam_Vls_[m2[i]]
      del_v[i] = abs(v2-v1)/v1


    indices = np.argsort(m1)
    m1 = m1[indices]
    m2 = m2[indices]
    del_v = del_v[indices]

    m1_old = -1

    del_v_list_f = []
    m1_list_f = []
    m2_list_f = []
	
    del_v_list = []
    m1_list = []
    m2_list = []
    for i in range(len(m1)):
      
      if m1[i] != m1_old or i==len(m1)-1:
	
	if len(del_v_list) != 0:
	  indices = np.argsort(del_v_list)
	  index = indices[0]
	
	  del_v_list_f.append(del_v_list[index])
	  m1_list_f.append(m1_list[index])
	  m2_list_f.append(m2_list[index])
	
	del_v_list = []
	m1_list = []
	m2_list = []
	
	if del_v[i] < 1.:
	  del_v_list.append(del_v[i])
	  m1_list.append(m1[i])
	  m2_list.append(m2[i])
	
	m1_old = m1[i]
	
      else:
	
	if del_v[i] < 1.:
	  del_v_list.append(del_v[i])
	  m1_list.append(m1[i])
	  m2_list.append(m2[i])    


    #print "2: ", len(m1_list_f), " --------- "


    for p in range(len(m1_list_f)):
      
      j = flag_indices[m1_list_f[p]]
      i = m2_list_f[p]

      v1 = nest_Vls[j]
      v2 = nam_Vls_[i]
      
      nest_dist[j]    = (v1/v2)*nam_dist[i]
      nest_flag[j]    = 2
      nest_dsource[j] = nam_name[i]
      

      
      #print v1, v2, abs(v2-v1)/v1   
    ################################################################################
    ################################################################################
    ################################################################################
    ###############################3333333333333####################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################


    flag_indices = np.where(nest_flag==0)[0]   # where there is not previous match
    nest_nest_flag = nest_nest[flag_indices]
    nest_sgl_flag  = nest_sgl[flag_indices]
    nest_sgb_flag  = nest_sgb[flag_indices]


    (m1, m2, d)  = spherematch.match_radec(nest_sgl_flag, nest_sgb_flag, nam_sgl, nam_sgb, 10., notself=False, nearest=False)

    del_v  = np.zeros((len(d),), dtype=float)
    for i in range(len(d)):
      v1 =  nest_Vls[flag_indices[m1[i]]] 
      v2 =  nam_Vls_[m2[i]]
      del_v[i] = abs(v2-v1)/v1


    indices = np.argsort(m1)
    m1 = m1[indices]
    m2 = m2[indices]
    del_v = del_v[indices]

    m1_old = -1

    del_v_list_f = []
    m1_list_f = []
    m2_list_f = []
	
    del_v_list = []
    m1_list = []
    m2_list = []
    for i in range(len(m1)):
      
      if m1[i] != m1_old or i==len(m1)-1:
	
	if len(del_v_list) != 0:
	  indices = np.argsort(del_v_list)
	  index = indices[0]
	
	  del_v_list_f.append(del_v_list[index])
	  m1_list_f.append(m1_list[index])
	  m2_list_f.append(m2_list[index])
	
	del_v_list = []
	m1_list = []
	m2_list = []
	
	if del_v[i] < 1.:
	  del_v_list.append(del_v[i])
	  m1_list.append(m1[i])
	  m2_list.append(m2[i])
	
	m1_old = m1[i]
	
      else:
	
	if del_v[i] < 1.:
	  del_v_list.append(del_v[i])
	  m1_list.append(m1[i])
	  m2_list.append(m2[i])    


    #print "3: ", len(m1_list_f), " --------- "


    for p in range(len(m1_list_f)):
      
      j = flag_indices[m1_list_f[p]]
      i = m2_list_f[p]

      v1 = nest_Vls[j]
      v2 = nam_Vls_[i]
      
      nest_dist[j]    = (v1/v2)*nam_dist[i]
      nest_flag[j]    = 3
      nest_dsource[j] = nam_name[i]


      #print v1, v2, abs(v2-v1)/v1   

    ################################################################################
    ################################################################################
    ###############################     444444444444444        #####################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################


    flag_indices = np.where(nest_flag==0)[0]   # where there is not previous match
    nest_nest_flag = nest_nest[flag_indices]
    nest_sgl_flag  = nest_sgl[flag_indices]
    nest_sgb_flag  = nest_sgb[flag_indices]


    (m1, m2, d)  = spherematch.match_radec(nest_sgl_flag, nest_sgb_flag, nam_sgl, nam_sgb, 20., notself=False, nearest=False)

    del_v  = np.zeros((len(d),), dtype=float)
    for i in range(len(d)):
      v1 =  nest_Vls[flag_indices[m1[i]]] 
      v2 =  nam_Vls_[m2[i]]
      del_v[i] = abs(v2-v1)/v1


    indices = np.argsort(m1)
    m1 = m1[indices]
    m2 = m2[indices]
    del_v = del_v[indices]

    m1_old = -1

    del_v_list_f = []
    m1_list_f = []
    m2_list_f = []
	
    del_v_list = []
    m1_list = []
    m2_list = []
    for i in range(len(m1)):
      
      if m1[i] != m1_old or i==len(m1)-1:
	
	if len(del_v_list) != 0:
	  indices = np.argsort(del_v_list)
	  index = indices[0]
	
	  del_v_list_f.append(del_v_list[index])
	  m1_list_f.append(m1_list[index])
	  m2_list_f.append(m2_list[index])
	
	del_v_list = []
	m1_list = []
	m2_list = []
	
	del_v_list.append(del_v[i])
	m1_list.append(m1[i])
	m2_list.append(m2[i])
	
	m1_old = m1[i]
	
      else:
	
	del_v_list.append(del_v[i])
	m1_list.append(m1[i])
	m2_list.append(m2[i])    


    #print "4: ", len(m1_list_f), " --------- "


    for p in range(len(m1_list_f)):
      
      j = flag_indices[m1_list_f[p]]
      i = m2_list_f[p]
      
      v1 = nest_Vls[j]
      v2 = nam_Vls_[i]
      
      nest_dist[j]    = (v1/v2)*nam_dist[i]
      nest_flag[j]    = 4
      nest_dsource[j] = nam_name[i]


      #print v1, v2, abs(v2-v1)/v1

    ################################################################################



    ################################################################################
    flag_indices = np.where(nest_flag==0)[0]   # where there is not previous match
    if len(flag_indices) > 0 : 
      print "Warning:", len(flag_indices), "groups have No distance associations ... "

    ################################################################################
    indices = np.argsort(nest_pgc)
    nest_nest    = nest_nest[indices]
    nest_sgl     = nest_sgl[indices]
    nest_sgb     = nest_sgb[indices]
    nest_gl      = nest_gl[indices]
    nest_gb      = nest_gb[indices]
    nest_Vls     = nest_Vls[indices]
    nest_pgc     = nest_pgc[indices]
    nest_d       = nest_d[indices]
    nest_dist    = nest_dist[indices]
    nest_flag    = nest_flag[indices]
    nest_dsource = nest_dsource[indices]
    
    nest_table = {'Nest':nest_nest, 'sgl':nest_sgl, 'sgb':nest_sgb, 'gl':nest_gl, 'gb':nest_gb, \
                  'Vls':nest_Vls, 'PGC1':nest_pgc, 'd':nest_d, 'dist':nest_dist, 'nest_flag': nest_flag, 'nest_dsource':nest_dsource}
    
    return nest_table
#################################################################
#################################################################

if __name__ == '__main__':
  
  nest_file = 'nests_VI38.1799extended_lob.csv'
  nam_file  = 'NAMtestV5_08.dat.csv'
  nest(nest_file, nam_file)





