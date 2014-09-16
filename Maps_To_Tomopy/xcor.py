import h5py
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy.signal import correlate2d
from PIL import Image


######
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

def xcor(data,element1,element2):
      xmatrix=data[element1,:,:,:]
      ymatrix=data[element2,:,:,:]

      finalmatrix=zeros(data.shape,dtype=float32)

      shift=zeros([2,data.shape[1]],dtype=int)
      for i in arange(len(data[0,:,0,0])-1):
            cor2dx=correlate2d(xmatrix[i,:,:],xmatrix[i+1,:,:],fillvalue=np.average(xmatrix[i,:,:]))
            position1=np.where(cor2dx==cor2dx.max())
            cor2dy=correlate2d(ymatrix[i,:,:],ymatrix[i+1,:,:],fillvalue=np.average(ymatrix[i,:,:]))
            position2=np.where(cor2dy==cor2dy.max())

            x1=position1[1][0]
            y1=position2[0][0]

            shift[0,i+1]=y1
            shift[1,i+1]=x1

            x=len(xmatrix[i,0,:])
            y=len(ymatrix[i,:,0])

            for j in arange(len(data[0,:,0,0])-1-i):
                  finalmatrix[i+j+1,:,:,:]=np.roll(data[i+j+1,:,:,:],x1-x+1,axis=2)
                  finalmatrix[i+j+1,:,:,:]=np.roll(data[i+j+1,:,:,:],y1-y+1,axis=1)
            print "Xcor", i
      return finalmatrix, shift
