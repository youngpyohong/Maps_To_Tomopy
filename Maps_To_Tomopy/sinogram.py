import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy import ndimage
from PIL import Image

def sinogram(data,thickness=10):
      projections=data.shape[0]
      sinogram=zeros([data.shape[0]*thickness,data.shape[2]],dtype=float32)
      for i in arange(projections):
            for j in arange(thickness):
                  sinogram[i*thickness+j,:]=data[i,data.shape[1]/2,:]
      plt.imshow(sinogram)
      show()
      return sinogram
