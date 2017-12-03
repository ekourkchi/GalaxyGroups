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
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py

from grouping_v06_8pc  import *
# **************************************
#    Global Variables
# **************************************

# **************************************
#### Don't Delete the following method
# **************************************
def readTree (filename):
  
    #skip_header=3
    table = np.genfromtxt(filename , delimiter=',', filling_values=0, names=True, dtype=None )
    
    
    return table 

# **************************************

#################################################################

#################################################################

def densityplot(root, ax, n): 
  
  densityaxCore(root, ax)
  ax.set_xlabel('PGC ID', fontsize=12)
  ax.set_ylabel(r'Inverse density ($L_{\odot}\/Mpc^{-3}$)', fontsize=12)
  #ax.set_title('ax1 title')
  
  yticks = [item for item in ax.get_yticks()]
  d = abs(yticks[1] - yticks[0])
  Rho = (10**root.logK) / ((4./3)*pi*((root.R_2t2*pi/2.)**3)) 
  y_max = max(ax.get_ylim()[1], 1./Rho+2*d)
  #y_max = 1./root.Rho_gr+d

  
  
  y = 1./Rho_cr(10**root.logK, 50)
  #line0, = ax.plot([0,n+1], [y,y], '-k', lw=2)
  #line0.set_dashes([8,3])
  
  
  Level = [1./root.Rho_gr, y_max]
  X = [root.dend_x, root.dend_x]
  line, = ax.plot(X, Level, '-r', lw=2)
  line.set_dashes([2, 2])
  ax.set_xlim([0, n+1])
  ax.set_ylim([0, y_max])
  
def densityaxCore(root, ax):
  
  
  if root.left!= None:
    
    

    
    if root.left.Rho_gr==0: left=0
    else: left = 1./root.left.Rho_gr
    
    if root.right.Rho_gr==0: right=0
    else: right = 1./root.right.Rho_gr
    
    if root.Rho_gr==0: center=0
    else: center = 1./root.Rho_gr
    
    
    
    
    X = [root.left.dend_x, root.left.dend_x, root.right.dend_x, root.right.dend_x]
    Level = [left, center, center, right]
    
    #print X
    #print Level
    ax.plot(X, Level, '-r')
    
    densityaxCore(root.left, ax)
    densityaxCore(root.right, ax)
    
#################################################################
def radiusplot(root, ax, n): 

 
  radiusaxCore(root, ax)
  ax.set_xlabel('PGC ID', fontsize=12)
  ax.set_ylabel(r'Radius ($Mpc$)', fontsize=12)
  #ax.set_title('ax1 title')
  
  yticks = [item for item in ax.get_yticks()]
  d = yticks[1] - yticks[0]
  y_max = root.Rgroup+d
  Level = [root.Rgroup, y_max]
  X = [root.dend_x, root.dend_x]
  line, = ax.plot(X, Level, '-b', lw=2)
  line.set_dashes([2, 2])
  ax.set_xlim([0, n+1])
  ax.set_ylim([0, y_max])
  
def radiusaxCore(root, ax):
  
  
  if root.left!= None:
    
    
    left = root.left.Rgroup
    right = root.right.Rgroup
    center = root.Rgroup
    

    X = [root.left.dend_x, root.left.dend_x, root.right.dend_x, root.right.dend_x]
    Level = [left, center, center, right]
    
    #print X
    #print Level
    ax.plot(X, Level, '-b')
    
    radiusaxCore(root.left, ax)
    radiusaxCore(root.right, ax)

#################################################################
#################################################################
def levelplot(root, ax, n): 
  
  levelaxCore(root, ax)

  ax.set_xlabel('PGC ID', fontsize=12)
  ax.set_ylabel("Level", fontsize=12)
  #ax.set_title('ax1 title')
  
  
  yticks = [item for item in ax.get_yticks()]
  d = yticks[1] - yticks[0]
  y_max = max(ax.get_ylim()[1], root.level+d)
  Level = [root.level, y_max]
  X = [root.dend_x, root.dend_x]
  line, = ax.plot(X, Level, '-k', lw=2)
  line.set_dashes([2, 2])
  ax.set_xlim([0, n+1])
  ax.set_ylim([0, y_max])

def levelaxCore(root, ax):
  
  
  if root.left!= None:
    
    
    left = root.left.level
    right = root.right.level
    center = root.level    

    X = [root.left.dend_x, root.left.dend_x, root.right.dend_x, root.right.dend_x]
    Level = [left, center, center, right]
    
    #print X
    #print Level
    ax.plot(X, Level, '-k')
    
    levelaxCore(root.left, ax)
    levelaxCore(root.right, ax)
    
    
 #################################################################

def Criticalradiusplot(root, ax): 
  
  X = []
  Level = []
  CriticalradiusCore(root, X, Level)
  X.extend([0])
  Level.extend([Level[len(Level)-1]])
  X.insert(0, X[0]+2)
  Level.insert(0, Level[0])
  
  line, = ax.plot(X, Level, '-g', linewidth=2)
  line.set_dashes([3,5]) 

  yticks = [item for item in ax.get_yticks()]
  d = yticks[1] - yticks[0]
  Rho = (10**root.logK) / ((4./3)*pi*((root.R_2t2*pi/2.)**3)) 
  y_max = max(ax.get_ylim()[1], root.R_2t2*pi/2.+d)
  ax.set_ylim([0, y_max])


  
def CriticalradiusCore(root, x, level):
  
  if root == None: return x, level
  
  elif root.left != None:
       
      x.append(root.dend_x)
      Rho = (10**root.logK) / ((4./3)*pi*((root.R_2t2*pi/2.)**3)) 
      #level.append(1./Rho_cr(10.**root.logK, 50))
      level.append(root.R_2t2*pi/2.)
      
      if root.left.right != None:# and root.right.dend_x - root.left.right.dend_x > 5:
	print root.right.nest, root.left.nest
      
      
      
      CriticalradiusCore(root.left, x, level)
  

#################################################################  
#################################################################

def Criticaldensityplot(root, ax): 
  
  X = []
  Level = []
  CriticaldensityaCore(root, X, Level)
  X.extend([0])
  Level.extend([Level[len(Level)-1]])
  X.insert(0, X[0]+2)
  Level.insert(0, Level[0])
  
  line, = ax.plot(X, Level, '-g', linewidth=2)
  line.set_dashes([3,5]) 

  yticks = [item for item in ax.get_yticks()]
  d = yticks[1] - yticks[0]
  Rho = (10**root.logK) / ((4./3)*pi*((root.R_2t2*pi/2.)**3)) 
  y_max = max(ax.get_ylim()[1], 1./Rho+d)
  ax.set_ylim([0, y_max])


  
def CriticaldensityaCore(root, x, level):
  
  if root == None: return x, level
  
  elif root.left != None:
       
      x.append(root.dend_x)
      Rho = (10**root.logK) / ((4./3)*pi*((root.R_2t2*pi/2.)**3)) 
      #level.append(1./Rho_cr(10.**root.logK, 50))
      level.append(1./Rho)
      
      if root.left.right != None :#and root.right.dend_x - root.left.right.dend_x > 5:
	print root.right.nest, root.left.nest
      
      
      
      CriticaldensityaCore(root.left, x, level)
  

#################################################################
 
def subGroup(root):
  
  list = []
  subGroupCore(root,list)
  return list

def subGroupCore(root,list):
  
  if root == None: 
    return list
  elif root.left != None:
  
    if root.left.right != None and root.right.dend_x - root.left.right.dend_x > 5:
      list.append(root)
  
    subGroupCore(root.left,list)
  

#################################################################

def dentplot(table, id, isHorizontal=True, filename=None):
   
   
   
   root = makeTree(table, id)
   #root.toString()
   dend_setX(root, 1)
   
   
   pgc, sgl, sgb = allLeaves(root)
   #print pgc 
   
   
   if isHorizontal:
      fig = py.figure(figsize=(15, 5), dpi=100)
      fig.subplots_adjust(wspace=0.2, top=0.95, bottom=0.15, left=0.05, right=0.98) 
      ax1 = fig.add_subplot(131)
      ax2 = fig.add_subplot(132)
      ax3 = fig.add_subplot(133)
      levelplot(root, ax1)
      radiusplot(root, ax2)
      densityplot(root, ax3)  
      
      labels = []
      for i in range(0,len(pgc)):
	  labels.append(str(int(pgc[i])))
      
      for ax in [ax1,ax2,ax3]:
         ax.set_xticks(range(1,len(pgc)+1))
         ax.set_xticklabels(labels, rotation=90, fontsize=6)
         ax.tick_params(top='off', bottom='off')
      
   else:
      fig = py.figure(figsize=(5, 8), dpi=100) # 8 -> 12
      fig.subplots_adjust(hspace=0.1, top=0.95, bottom=0.10, left=0.15, right=0.95)   
      #ax1 = fig.add_subplot(311)
      #ax2 = fig.add_subplot(312)
      #ax3 = fig.add_subplot(313)   
      ax2 = fig.add_subplot(211)
      ax3 = fig.add_subplot(212)
      
      levelplot(root, ax3, len(pgc))
      
      
      radiusplot(root, ax2, len(pgc))
      Criticalradiusplot(root, ax2)
      
      #densityplot(root, ax1,len(pgc) ) 
      #Criticaldensityplot(root, ax1)
      
      
      ax3.tick_params( bottom='off')
      ax2.tick_params( bottom='off')
      #ax1.tick_params( bottom='off')
      #ax1.set_xlabel('')
      ax2.set_xlabel('')
      py.setp(ax2.get_xticklabels(), visible=False)
      #py.setp(ax1.get_xticklabels(), visible=False)
      
      labels = [item.get_text() for item in ax3.get_xticklabels()]
      labels = []
      ax3.set_xticks(range(1,len(pgc)+1))

      for i in range(0,len(pgc)):
	  labels.append(str(int(pgc[i])))
      ax3.set_xticklabels(labels, rotation=90, fontsize=6)
      
      ax4 = ax3.twiny()
      xticks = [item for item in ax2.get_xticks()]
      ax4.set_xticks(xticks)
      ax4.tick_params(bottom='off')
      py.setp(ax4.get_xticklabels(), visible=False)
      


      # Turn off axis lines and ticks of the big subplot
      #ax.spines['top'].set_color('none')
      #ax.spines['bottom'].set_color('none')
      #ax.spines['left'].set_color('none')
      #ax.spines['right'].set_color('none')

   

   if filename == None:
     plt.show()
   else:
     plt.gcf() # get current figure
     plt.savefig(filename, dpi=600)  # the best view and size happens in eps format
   

#################################################################

#################################################################

def singleGroupPlot(inFile, treeTable, ID, R_min, R_max, New=True, filename=None):
  
  
  
  
  
  table = np.genfromtxt( inFile , delimiter=',', usecols=(0,1,2,3,4,5,6,7), \
			filling_values=0, names=True, dtype=None)
  id   = table['PGC']
  sgl  = table['sgl']
  sgb  = table['sgb']

  
  if New :
	    fig = plt.figure(figsize=(7.5, 7), dpi=100)
	    ax = fig.add_axes([0.13, 0.1, 0.85,  0.85]) 
	    ax.xaxis.set_major_locator(MultipleLocator(20))
	    ax.yaxis.set_major_locator(MultipleLocator(20))
	    ax.xaxis.set_minor_locator(MultipleLocator(5))
	    ax.yaxis.set_minor_locator(MultipleLocator(5))  
	    
	    
	    
	    plt.minorticks_on()
	    plt.tick_params(which='major', length=7, width=1.5)
	    plt.tick_params(which='minor', length=4, color='#000033', width=1.0)  
	    
	    
	    xmin = 5*floor((sglV+R_max)/5.)+3
	    xmax = 5*ceil((sglV-R_max)/5.)-3
	    ymax = 5*floor((sgbV+R_max)/5.)+3
	    ymin = 5*ceil((sgbV-R_max)/5.)-3
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
	    
	    
	    plt.plot(sgl, sgb, '.', markersize = 1, color='#696969')  # gray
	  
  
  G_list = []
  root = makeTree(treeTable, ID)
  G_list.append(NodeLeaf(root))
  
  dend_setX(root, 1)
  subGroupsRoot = subGroup(root)
  n = len(subGroupsRoot)
  
  #print n
  if n>0:
      G_list.append(NodeLeaf(subGroupsRoot[n-1].left))
      G_list.append(NodeLeaf(subGroupsRoot[n-1].right))
      
      for i in range(0, n-1): 
		      G_list.append(NodeLeaf(subGroupsRoot[i].right))
		      #print subGroupsRoot[i].id
  
  NoGroups = len(G_list)
  if NoGroups>1: 

	    
	    print "Number of groups: ", NoGroups
	    

	    for i in range(0, NoGroups):  # for all groups
	    
	      
	      if i == 0:
	         random.seed( G_list[0][0].nest )
	      Red, Blue, Green = random.random(), random.random(), random.random()
	      #Gcolor = i*(red-violet)/NoGroups+violet
	      #Red, Blue, Green = col.wavelength_to_rgb(Gcolor, 0.5)   ## 0 .. 255
	      
	      r = G_list[i][0].R_theta
	      #d_theta = 0.2*pi/r
	      d_theta = 0.001
	      theta = np.arange(0,2*pi,d_theta)
	      Circlx = r*np.cos(theta) + G_list[i][0].l
	      Circly = r*np.sin(theta) + G_list[i][0].b
	      line, = plt.plot(Circlx, Circly, '-', markersize = 2, color=(Red, Blue, Green)) 
	      #line.set_dashes([8, 4, 2, 4, 2, 4]) 
	      #[8, 4, 2, 4, 2, 4] means

	      #8 points on, (dash)
	      #4 points off,
	      #2 points on, (dot)
	      #4 points off,
	      #2 points on, (dot)
	      #4 points off.
	      line.set_dashes([8, 3]) 
	      
	      
	      
	      
	      
	      
	      sRA  = []
	      sDEC = []
	      
	      
	      #plt.text(G_list[i][0].l, G_list[i][0].b, G_list[i][0].nest, fontsize=8, color=(Red, Blue, Green))
	      for j in range(1, len(G_list[i])):
		  sRA.append(G_list[i][j].l)
		  sDEC.append(G_list[i][j].b)
	      plt.plot(sRA, sDEC, 'o', markersize = 5, color=(Red, Blue, Green), markeredgecolor = 'none')
	    
  
  
  #plt.show()
#################################################################

#################################################################

if __name__ == '__main__':
  
  


  
  
  #nameRoot = 'final.virgo'
  inTree  = 'south.ML.60.v14.8pc.tree'
  
  
  print "Reding the tree ... "
  table = readTree(inTree)
  print "The tree is read ... "

  
  R_min = 6.0
  R_max = 100000000.
  n_gal = 500000
  inFile  = 'AllSky.north.v9.csv'
  table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
  galList = readgalList(table, R_min, R_max, n_gal, useDCF2=True)
  
  id =  300043495
  #################################################################
  file = 'north.ML.64.v12.group'
  Gtable = np.genfromtxt(file , delimiter='|', filling_values="-100000", names=True, dtype=None )
  flag = Gtable['flag']
  No = Gtable['No_Galaxies']
  R2t_dyn = Gtable['R2t_dyn']
  R2t_lum = Gtable['R2t_lum']
  sigmaP_dyn = Gtable['sigmaP_dyn']
  sigmaP_lum = Gtable['sigmaP_lum']
  mDist = Gtable['mDist']
  Vls = Gtable['Vls']
  gl = Gtable['gl']
  gb = Gtable['gb']
  #################################################################
  
  
  
  
  dentplot(table, 1000000090, isHorizontal=False) 
  
  



  plt.show()
