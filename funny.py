
import matplotlib.pyplot as plt
import numpy as np
from math import *


# **************************************
if __name__ == '__main__':
  
  
  
  x = []
  y = []
  force = []

  p = 0 
  for i in np.arange(0.001,1,0.01):
    for j in np.arange(0.001,1,0.01):
      x.append(i)
      y.append(j)
      force.append(sin(atan(j/i))/(i**2+j**2))
      p += 1

  
  indices = np.argsort(force)
  x = np.asarray(x)
  y = np.asarray(y)

  m = 1500  # number of added blocks
  indices = indices[p-m:p]
  
  XX =  x[indices]
  YY =  y[indices]
  
      
  plt.plot(XX, YY, 'o', markersize = 1)
  plt.show()

  
  # The end
  
  