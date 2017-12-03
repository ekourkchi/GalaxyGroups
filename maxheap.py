#!/home/ehsan/Ureka/Ureka/variants/common/bin/python


import numpy as np
from math import *
from copy import *

class heapNode: 
  
  key = None
  ID = None
  flag = False
  
  
  def __init__(self, key, ID):
    
    self.key = key
    self.ID = ID
  
  def toString(self):
    print self.key, self.ID, self.flag
  
# *********************************************
class maxHeap:
  
  size = 0   # Number of current elements
  array = []
# *****************
  
  def __init__(self):
    
    self.size = 0
    self.array = []
# *****************
  
  def push(self, key, ID):
    #print "push:", key, ID, self.size
    
    
    
    newNode = heapNode(key, ID)
    self.array.append(newNode)
    child = self.size
    
    while child > 0:
      parent = (child+1)/2-1
      if self.array[child].key > self.array[parent].key:
	
	self.array[parent], self.array[child] = self.array[child], self.array[parent]
	child = parent
      else:
	break
    
    #for i in range(0,self.size+1):
      #print self.array[i].key
      
    
    
    self.size+=1
    return 0
# *****************
   
  def lrmax(self, left, right):
     
     if right <= self.size-1:
	  if self.array[left].key >= self.array[right].key:
	    return left
	  else:
	    return right
     elif left <= self.size-1:
       return left
     else:
       return 0
     
# *****************

  def pop(self):
     
     if self.size == 0 :
       print "\n[Error] No elements in the mean Heap ...\n"
       return None
     
     N = self.size
     output = self.array[0]
     self.array[0] = self.array[N-1]
     parent = 0
     
     while parent <= N-1:
	left = 2*parent+1
	right = 2*parent+2
	
	child = self.lrmax(left, right)

	if child != 0:
	  if self.array[child].key >= self.array[parent].key:
            self.array[parent], self.array[child] = self.array[child], self.array[parent]
	    parent = child
	  else:
	    break
	else:
	   break
	   
     self.array.pop(N-1)
     self.size -= 1
     return output  

# *****************
     
  def setFlag(self, key):
      if self.size == 0 :
       print "\n[Error] No elements in the mean Heap ...\n"
       return False  
      
      for i in range(0, self.size):
	if self.array[i].key == key:
	  self.array[i].flag = True
# *****************
  
  def peek(self):
      if self.size == 0 :
       print "\n[Error] No elements in the mean Heap ...\n"
       return None  
      
      else:
	return self.array[0]
     
# *****************
  """
  This method removes heap elements which have the same id as the input ID
  The number of removed elements would be returned
  """
  
  def remove(self, ID):
      
      boolean = 0
      if self.size == 0 :
       #print "\n[Error] No elements in the mean Heap ...\n"
       return boolean  
      
      else:

	  i = 0
	  while i < self.size:
		
		# ID would be the object ID
		if self.array[i].ID == ID:
		        
			parent = i
			N = self.size
			self.array[parent] = self.array[N-1] 
			while parent <= N-1:
			    left = 2*parent+1
			    right = 2*parent+2
		    
			    child = self.lrmax(left, right)

			    if child != 0:
			      if self.array[child].key >= self.array[parent].key:
				  self.array[parent], self.array[child] = self.array[child], self.array[parent]
				  parent = child
			      else:
				  break
			    else:
			      break
			
			self.array.pop(N-1)
			self.size -= 1
			boolean+=1
			i-=1  # The new item must be checked again
	    
		i+=1
	  return boolean
     
   
# ***************** 
  
  def Size(self):  return self.size
  
# *****************
  
  def toString(self):
      
      for i in range(0,self.size):
	self.array[i].toString();
  
# *********************************************
# *********************************************

if __name__ == '__main__':
  
  myHeap = maxHeap()
  
  myHeap.push(4, "e4")
  myHeap.push(7, "e7")
  myHeap.push(2, "e2")
  myHeap.push(6, "e6")
  myHeap.push(8, "e7")
  myHeap.push(5, "e5")
  myHeap.push(3, "e7")
  print "\n", myHeap.Size()
  print myHeap.remove("e5")


  print "\n", myHeap.Size()
  
  while myHeap.Size()>0:
     myHeap.pop().toString()
  #print myHeap.peek().key


