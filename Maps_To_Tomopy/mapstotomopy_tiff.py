# -*- coding: utf-8 -*-
import tomopy
import numpy as np
from scipy import *

from Maps_To_Tomopy.readwritetxt import readtxt, writetxt

from Maps_To_Tomopy.sinogram import sinogram
from PIL import Image
import os

from Maps_To_Tomopy.xtomo_io_f import xtomo_reader_f, xtomo_writer_f, xtomo_reader_config
from tomopy.xtomo.xtomo_io import xtomo_writer
import sys


# Read HDF5 file.
#data, white, dark, theta = tomopy.xtomo_reader('demo/demo3.h5',
#                                               slices_start=0,
#                                               slices_end=1)


RH=[u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted24.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted25.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted26.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted27.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted28.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted29.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted30.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted31.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted32.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted33.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted34.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted35.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted36.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted37.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted38.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted39.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted40.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted41.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted42.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted43.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted44.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted45.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted46.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted47.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted48.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted49.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted50.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted51.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted52.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted53.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted54.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted55.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted56.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted57.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted58.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted59.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted60.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted61.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted62.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted63.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted64.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted65.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted66.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted67.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted68.tiff', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted69.tiff']

#### temp
temp=np.asarray(Image.open(RH[0]),dtype=float32)
y,x = temp.shape
projections=len(RH)
data=zeros([projections,y,x])
theta=zeros(projections)
for i in arange(projections):
    print i
    data[i,:,:]=np.asarray(Image.open(RH[i]),dtype=float32)
    theta[i]=i*6

for i in arange(29):
    theta[i]+=3

sino=sinogram(data)

### saves if one specified to save sinogram.
### The path of the sinogram will be printed on the screen

d = tomopy.xtomo_dataset(log='debug')
data[data == inf] = 1
d.dataset(data, theta=theta*np.pi/180)

#d.dataset(data, theta)
#d.data = (np.exp(-0.0001*d.data)).astype('float32')

#print d.data.shape
#tomopy.xtomo_writer(d.data, 'tmp/test_', axis=1, overwrite=True)
#d.normalize()
#d.correct_drift()
#d.phase_retrieval()
#d.data = (np.exp(-0.0001*d.data)).astype('float32')
#d.correct_drift()
#tomopy.xtomo_writer(d.data, 'tmp/test_', axis=1, overwrite=True)
#d.data[d.data is np.nan]=1.
#whereAreNaNs = np.isnan(d.data)
#d.data[whereAreNaNs] = 1.
#print d.data
#d.center=982
#d.theta *= np.pi/180
#d.optimize_center()

#d.diagnose_center()
#d.mlem()
#d.art()

# Write to stack of TIFFs.
#tomopy.xtomo_writer(d.data_recon, 'tmp/gridrec/test_', axis=0)
d.optimize_center()
for i in arange(10):
      d.center=40+i
      d.mlem()
      xtomo_writer(d.data_recon, 'tmp/tmp/'+str(d.center)+'/test_'+str(d.center), axis=0)

