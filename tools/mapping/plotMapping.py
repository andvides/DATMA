#!/usr/bin/env python
"""
0	611283	
1	1861243	
2	1741633	1741828	
"""

import matplotlib.pyplot as plt
import numpy as np
import os,argparse

def getList(result,x,y,i):
  fileIter = open(result)
  #lista =set()
  for items in fileIter:
      	temp = items.split()
     	x.append(int(temp[1]))
     	y.append(i+1)
  fileIter.close()

# Fixing random state for reproducibility
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot mapping')
  parser.add_argument('-r', '--r', help='result file', required=True)

  args = parser.parse_args()

  files=args.r.split(',')
  x=[]
  y=[]
  runs=len(files)
  #colors=['b','r','g','y']	
  for i in range (runs):
  	getList(files[i],x,y,i)
        
  plt.scatter(x, y, marker='_')
  plt.show()

