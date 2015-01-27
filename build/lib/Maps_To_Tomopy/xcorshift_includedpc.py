import h5py
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy.signal import correlate2d
from PIL import Image
import subpixelshift
from subpixelshift import *

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
#from Maps_To_Tomopy.xcorshift import xcor_shift
######

path=[u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0219.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0220.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0221.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0222.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0223.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0224.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0225.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0226.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0227.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0228.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0229.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0230.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0247.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0248.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0249.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0250.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0251.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0252.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0253.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0254.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0255.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0256.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0257.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0258.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0276.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0277.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0278.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0279.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0280.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0281.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0282.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0283.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0284.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0285.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0286.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0287.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0288.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0289.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0290.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0291.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0292.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0293.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0294.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0295.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0296.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0297.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0298.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0299.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0300.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0301.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0302.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0303.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0304.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0305.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0306.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0307.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0308.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0309.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0310.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0311.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0312.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0313.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0314.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0315.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0316.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0317.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0318.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0319.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0320.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0321.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0352.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0353.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0354.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0355.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0356.h5', u'/Users/youngpyohong/Documents/Work/2014-2/hong/hong/30ms/2xfm_0357.h5']

def xyshift(path):
      path1=path[24:70]
      path1=path1[0:29]
      f1=h5py.File(path1[0],'r')
      f2=f1["MAPS"]["XRF_roi"][...]


      ## element to use for crosscorrelation horizontally
      element1=15
      ## element to use for xcor vertically
      element2=18
      ## dpc
      element4=7
      ## element for output
      element3=18

      matrix1=zeros([len(path1),24,38,151])
      matrix4=zeros([len(path1),19,38,151])
      matrix2=zeros([len(path1),24,38,151])
      matrix3=zeros([len(path1),24,38,151])
      theta=zeros(len(path1))
      
      for i in arange(len(path1)):
            fp=h5py.File(path1[i],'r')["MAPS"]
            f=fp["XRF_roi"][...]
            matrix2[i,:,:,:]=f
            matrix3[i,:,:,:]=f
            matrix1[i,:,:,:]=f
            matrix4[i,:,:,:]=fp["scalers"][...]
            temp=fp["extra_pvs_as_csv"][99]
            theta[i]=temp[temp.rfind(",")+2:]

      matrix2[:,element3,:,:]=np.roll(matrix2[:,element3,:,:],40,axis=2)
      matrix3[:,element3,:,:]=np.roll(matrix1[:,element3,:,:],40,axis=2)     

      shift=zeros([2,len(path1)],dtype=int)

      for i in arange(len(path1)-1):
            cor2d=correlate2d(matrix4[i,element4,:,:],matrix4[i+1,element4,:,:],fillvalue=np.average(matrix4[i,element4,:,:]))
            b1=np.where(cor2d==cor2d.max())
            cor2d=correlate2d(matrix1[i,element2,:,:],matrix1[i+1,element2,:,:],fillvalue=np.average(matrix1[i,element2,:,:]))
            b2=np.where(cor2d==cor2d.max())
            x1=b1[1][0]
            y1=b2[0][0]

            x=151
            y=38

            for j in arange(len(path1)-1-i):
                  matrix2[i+j+1,element3,:,:]=np.roll(matrix2[i+j+1,element3,:,:],x1-x+1,axis=1)
                  matrix2[i+j+1,element3,:,:]=np.roll(matrix2[i+j+1,element3,:,:],y1-y+1,axis=0)

            print y1,x1
 
            shift[0,i+1]=y1
            shift[1,i+1]=x1
##      for i in arange(len(path1)):
##            j=Image.fromarray(matrix2[i,element2,8:30,:].astype(np.float32))
##            if i>=10:
##                  j.save("/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted"+str(i)+".tiff")
##            else:
##                  j.save("/Users/youngpyohong/Documents/Work/2014-2/hong/hong/projections/shifted0"+str(i)+".tiff")
##      ## wkatl
      return shift,matrix2[:,element3,8:30,:],theta
###### Addition
a,b,theta=xyshift(path)
##for i in arange(b.shape[0]):
##      plt.imshow(b[i,:,:])
##      show()
theta[0]=-177
theta[6]=-135
theta[11]=-105
d = tomopy.xtomo_dataset(log='debug')
data=b
data[data == inf] = 1
#d.dataset(data, theta=theta*np.pi/180)
d.dataset(data,theta=theta)
d.data = (np.exp(-0.0001*d.data)).astype('float32')

#d.diagnose_center()

for i in arange(10):
      d.center=40+i
      d.gridrec()
      xtomo_writer(d.data_recon, 'tmp/tmp/'+str(d.center)+'/test_'+str(d.center), axis=0)

##      try:
##            j.save(os.getcwd()+"/tmp/tmp"+channel+"/projections"+"/"+channel+"_projection"+str(i)+".tiff")
##      except IOError:
##            os.makedirs(os.getcwd()+"/tmp/"+channel+"/projections")
##            j.save(os.getcwd()+"/tmp/"+channel+"/projections"+"/"+channel+"_projection"+str(i)+".tiff")
