# -*- coding: utf-8 -*-
import tomopy
import numpy as np
from scipy import *
import matplotlib as plt
import pylab
from pylab import *
from Maps_To_Tomopy.readwritetxt import readtxt, writetxt

from Maps_To_Tomopy.sinogram import sinogram
from PIL import Image
import os

from Maps_To_Tomopy.xtomo_io_f import xtomo_reader_f, xtomo_writer_f, xtomo_reader_config
from tomopy.xtomo.xtomo_io import xtomo_writer
from Maps_To_Tomopy.xcor import xcor
import sys


# Read HDF5 file.
#data, white, dark, theta = tomopy.xtomo_reader('demo/demo3.h5',
#                                               slices_start=0,
#                                               slices_end=1)

### Reads argument so that one can run other configuration files
try:
    textfile=sys.argv[1]
except IndexError:
    textfile=None
textfile="configuration_jake.txt"

### choose either reading from configuration file or reading script
try:
      file_name_f,projections_start,projections_end,exclude_numbers,element,reconstructiontype,sinogramsave,cor_number=readtxt(textfile)
      print "reading txt file"
      data_orig,theta,channelname=xtomo_reader_config(textfile)
      print "yes"
except IOError:
      data_orig,theta,channelname=xtomo_reader_f("/Users/youngpyohong/Documents/Work/2014-1/glassrod_W_Au/h5files/tomodata/2xfm",projections_start=172,projections_end=190)
      print "reading script beep beep"

### create a text file with channel name
print textfile
writetxt(channelname)


###
if len(cor_number)!=0:
    data_orig,shift=xcor(data_orig,cor_number[0],cor_number[1])
    print "cross correlation has been applied"


projections=int(projections_end)-int(projections_start)+1-len(exclude_numbers)

###$$$ temporary: shifts
##diff=array([ -4,   5,   8,   4,   3,  -3, -10, -10,  -7,   3,   1,   6,   6,
##         4,   0,   1,   3,   0,  -9], dtype=int32)
##
##diff=diff[int(projections_start)-172:int(projections_end)-172+1]

### counts how many elements to be reconstructed
count=element.shape[0]

theta=zeros(theta.shape)
for i in arange(theta.shape[0]):
    theta[i]=i*3



### making a loop for each element
for numb in arange(count):
      elem=element[numb]
      print elem
      ### select specific element 3d array (projections, x pixel, y pixel)
      data=data_orig[elem,:,:,:]

      
##      for i in arange(projections):
##            olddata=data[i,:,:]
##
##            data[i,:,:]=np.roll(olddata[:,:],diff[i],axis=1)

      sino=sinogram(data,thickness=0)

### saves if one specified to save sinogram.
### The path of the sinogram will be printed on the screen
      if sinogramsave==True:
            j=Image.fromarray(sino.astype(np.float32))
            try:
                  j.save(os.getcwd()+"/sinogram/"+channelname[elem]+"_sinogram.tiff")
                  print "saving sinogram... "
                  print os.getcwd()+"/sinogram/"+channelname[elem]+"_sinogram.tiff"
            except IOError:
                  os.mkdir(os.getcwd()+"/sinogram")
                  j.save(os.getcwd()+"/sinogram/"+channelname[elem]+"_sinogram.tiff")
                  print "saving sinogram... "
                  print os.getcwd()+"/sinogram/"+channelname[elem]+"_sinogram.tiff"
#data=np.expand_dims(data[:,25,:],axis=1)
# Xtomo object creation and pipeline of methods.  
      d = tomopy.xtomo_dataset(log='debug')
      data[data == inf] = 1
      whereAreNaNs = np.isnan(data)
      data[whereAreNaNs] = 1.
      #d.dataset(data, theta=theta)
      d.dataset(data[30:,:,:], theta=theta[30:]*np.pi/180)

#print d.data.shape
#tomopy.xtomo_writer(d.data, 'tmp/test_', axis=1, overwrite=True)
#d.normalize()
#d.correct_drift()
#d.phase_retrieval()
#d.data = (np.exp(-0.0001*d.data)).astype('float32')
#d.correct_drift()
#tomopy.xtomo_writer(d.data, 'tmp/test_', axis=1, overwrite=True)
#d.data[d.data is np.nan]=1.

#print d.data
#d.center=982
#d.theta *= np.pi/180
#d.optimize_center()

#d.diagnose_center()
#d.mlem()
#d.art()

# Write to stack of TIFFs.
#tomopy.xtomo_writer(d.data_recon, 'tmp/gridrec/test_', axis=0)

      #d.center=959

      print d.theta
      print reconstructiontype, channelname[elem]
      if reconstructiontype=="mlem":
            d.mlem()
      if reconstructiontype=="art":
            d.art()
      if reconstructiontype=="gridrec":
            d.gridrec()

      xtomo_writer(d.data_recon, 'tmp/'+str(channelname[elem])+'/'
                            +reconstructiontype+'/test_'+str(d.center), axis=0)
