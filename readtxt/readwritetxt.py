import string
from scipy import *

def readtxt():
      f=open("test.txt")
      l=f.readlines()
      for i in arange(len(l)):
            if string.find(l[i],"name")!=-1:
                  filename=l[i+1][:-1]
            if string.find(l[i],"start")!=-1:
                  projections_start=l[i+1][:-1]
            if string.find(l[i],"end")!=-1:
                  projections_end=l[i+1][:-1]
            if string.find(l[i],"element")!=-1:
                  if string.find(l[i],"number")!=-1:
                        element=int(l[i+1][:-1])
      return filename,projections_start,projections_end,element
