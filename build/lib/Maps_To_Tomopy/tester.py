# -*- coding: utf-8 -*-
import tomopy
import numpy as np
from scipy import *
from scipy.signal import correlate2d
from Maps_To_Tomopy.readwritetxt import readtxt, writetxt

from Maps_To_Tomopy.sinogram import sinogram
from PIL import Image
import os

from Maps_To_Tomopy.xtomo_io_f import xtomo_reader_f, xtomo_writer_f, xtomo_reader_config
from tomopy.xtomo.xtomo_io import xtomo_writer
from Maps_To_Tomopy.xcor import xcor
import sys

textfile="configuration_jake.txt"

data_orig,theta,channelname=xtomo_reader_config(textfile)

for i in arange(theta.shape[0]):
      theta[i]=i*3

data=data_orig[37,:,:,:]
data[data==inf]=1
data[np.isnan(data)]=1

for i in arange(data.shape[0]-1):
      temp=data[i+1,:,:]
      cor2d=correlate2d(data[i,:,:],temp,fillvalue=np.average(temp))
      position=np.where(cor2d==cor2d.max())

      x1=position[1][0]


      data[i+1,:,:]=np.roll(temp,x1-data.shape[2]+1,axis=1)


sino=sinogram(data)

d=tomopy.xtomo_dataset(log='debug')
d.dataset(data, theta=theta*np.pi/180)

d.optimize_center()

d.mlem()
xtomo_writer(d.data_recon, 'tmp/test/test_',axis=0)
