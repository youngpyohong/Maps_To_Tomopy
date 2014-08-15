import string
from scipy import *

def readtxt(textfile=None):
      if textfile==None:
          textfile="configuration.txt"
      f=open(textfile)
      l=f.readlines()
      for i in arange(len(l)):
            if string.find(l[i],"file")!=-1:
                  filen=l[i+1][:-1]
            if string.find(l[i],"folder")!=-1:
                  foldern=l[i+1][:-1]
            if string.find(l[i],"start")!=-1:
                  projections_start=l[i+1][:-1]
            if string.find(l[i],"end")!=-1:
                  projections_end=l[i+1][:-1]
            if string.find(l[i],"element")!=-1:
                  if string.find(l[i],"number")!=-1:
                        element=int(l[i+1][:-1])
            if string.find(l[i],"reconstruct")!=-1:
                  rtype=l[i+1][:-1]
      filename=foldern+filen
      return filename,projections_start,projections_end,element,rtype
