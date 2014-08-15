# -*- coding: utf-8 -*-
import numpy as np
from scipy import misc
import shutil
import os
import h5py
from PIL import Image
import matplotlib.pyplot as plt
from pylab import *
from pylab import *
import readwritetxt
from readwritetxt import readtxt

def xtomo_reader_f(file_name_f,
                 projections_start=None,
                 projections_end=None,
                 projections_step=None,
                 slices_start=None,
                 slices_end=None,
                 slices_step=None,
                 pixels_start=None,
                 pixels_end=None,
                 pixels_step=None,
                 white_start=None,
                 white_end=None,
                 dark_start=None,
                 dark_end=None):
    """
    Read Data Exchange HDF5 file.
    
    Parameters
    ----------
    file_name : str
        Input file.
    
    projections_start, projections_end, projections_step : scalar, optional
        Values of the start, end and step of the projections to
        be used for slicing for the whole data.
    
    slices_start, slices_end, slices_step : scalar, optional
        Values of the start, end and step of the slices to
        be used for slicing for the whole data.
    
    pixels_start, pixels_end, pixels_step : scalar, optional
        Values of the start, end and step of the pixels to
        be used for slicing for the whole data.
    
    white_start, white_end : scalar, optional
        Values of the start and end of the
        slicing for the whole white field shots.
    
    dark_start, dark_end : scalar, optional
        Values of the start and end of the
        slicing for the whole dark field shots.
        
    Examples
    --------
    - Import data, white-field, dark-field and projection angles
      from HDF5 file:
        
        >>> import tomopy
        >>> 
        >>> # Load data
        >>> myfile = 'demo/data.h5'
        >>> data, white, dark, theta = tomopy.xtomo_reader(myfile)
        >>>
        >>> # Image data
        >>> import pylab as plt
        >>> plt.figure()
        >>> plt.imshow(data[:, 0, :])
        >>> plt.show()
    
    - Import only 4th slice from HDF5 file:

        >>> import tomopy
        >>> 
        >>> # Load data
        >>> myfile = 'demo/data.h5'
        >>> data, white, dark, theta = tomopy.xtomo_reader(myfile,  slices_start=4, slices_end=5)
        >>> 
        >>> # Image data
        >>> import pylab as plt
        >>> plt.figure()
        >>> plt.imshow(data[:, 0, :])
        >>> plt.show()
    """

    # Start working on checks and stuff.
    print "working fine"
    try:
      file_name_f,projections_start,projections_end,elemen,dummy,dummy2=readtxt()
      print "reading txt"
    except IOError:
      print "reading script"
    projections_start=int(projections_start)
    projections_end=int(projections_end)
    projections = projections_end-projections_start+1
    global f
    f=h5py.File(file_name_f + "_0" + str(projections_start) + ".h5")
    print f.filename
    print f['exchange']['images']
    num_channels, num_y, num_x = f['exchange']['images'].shape
    theta = zeros(projections)
    channelname=f['exchange']['images_names']
    data = zeros([num_channels, projections, num_y, num_x-2])
    
    for i in arange(projections):

          ### theta


          ### data
          projection_number = projections_start + i
          file_name = file_name_f + "_0" + str(projection_number)+".h5"
          file_name = os.path.abspath(file_name)
          f = h5py.File(file_name,"r")
          hdfdata = f["/exchange/images"]

          temp=f["MAPS"]["extra_pvs_as_csv"][99]
          theta[i]=temp[temp.rfind(",")+2:]
          
          for j in arange(num_channels):
                data[j,i,:,:]=hdfdata[j,:,:-2]
    # Start reading data.
    #f = h5py.File(file_name, "r")
    #hdfdata = f["/exchange_2/data"]

    f.close()
    
    return data, theta, channelname



def xtomo_writer_f(data, output_file=None, x_start=0,
                 digits=3, axis=0, overwrite=False, 
                 precision=True):
    """ 
    Write 3-D data to a stack of tif files.

    Parameters
    -----------
    output_file : str, optional
        Name of the output file.

    x_start : scalar, optional
        First index of the data on first dimension
        of the array.

    digits : scalar, optional
        Number of digits used for file indexing.
        For example if 4: test_XXXX.tiff
        
    axis : scalar, optional
        Imaages is read along that axis.
        
    overwrite: bool, optional
        if overwrite=True the existing data in the
        reconstruction folder will be overwritten
        
    precision : bool, optional
        Export data type precision. if True it 
        saves 32-bit precision. Otherwise it
        uses 8-bit precision.
    
    Notes
    -----
    If file exists, saves it with a modified name.
    
    If output location is not specified, the data is
    saved inside ``recon`` folder where the input data
    resides. The name of the reconstructed files will
    be initialized with ``recon``
    
    Examples
    --------
    - Save sinogram data:
        
        >>> import tomopy
        >>> 
        >>> # Load data
        >>> myfile = 'demo/data.h5'
        >>> data, white, dark, theta = tomopy.xtomo_reader(myfile)
        >>> 
        >>> # Save data
        >>> output_file='tmp/slice_'
        >>> tomopy.xtomo_writer(data, output_file, axis=1)
        >>> print "Images are succesfully saved at " + output_file + '...'
        
    - Save first 16 projections:
        
        >>> import tomopy
        >>> 
        >>> # Load data
        >>> myfile = 'demo/data.h5'
        >>> data, white, dark, theta = tomopy.xtomo_reader(myfile, projections_start=0, projections_end=16)
        >>> 
        >>> # Save data
        >>> output_file='tmp/projection_'
        >>> tomopy.xtomo_writer(data, output_file, axis=0)
        >>> print "Images are succesfully saved at " + output_file + '...'
        
    - Save reconstructed slices:
        
        >>> import tomopy
        >>> 
        >>> # Load data
        >>> myfile = 'demo/data.h5'
        >>> data, white, dark, theta = tomopy.xtomo_reader(myfile)
        >>> 
        >>> # Perform reconstruction
        >>> d = tomopy.xtomo_dataset(log='error')
        >>> d.dataset(data, white, dark, theta)
        >>> d.center = 661.5
        >>> d.gridrec()
        >>> 
        >>> # Save data
        >>> output_file='tmp/reconstruction_'
        >>> tomopy.xtomo_writer(d.data_recon, output_file, axis=0)
        >>> print "Images are succesfully saved at " + output_file + '...'
    """
    if output_file == None:
        output_file = "tmp/img_" 
    output_file =  os.path.abspath(output_file)
    dir_path = os.path.dirname(output_file)
    
    # Remove TIFF extension if there is.
    if (output_file.endswith('tif') or
        output_file.endswith('tiff')) :
            output_file = output_file.split(".")[-2]
  
    if overwrite:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
            
    # Create new folders.
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    # Select desired x from whole data.
    num_x, num_y, num_z = data.shape
    if axis == 0:
        x_end = x_start+1
    elif axis == 1:
        x_end = x_start+1
    elif axis == 2:
        x_end = x_start+1

    # Write data.
    file_index = ["" for x in range(digits)]
    for m in range(digits):
        file_index[m] = '0' * (digits - m - 1)
    ind = range(x_start, x_end)
    for m in range(len(ind)):
        for n in range(digits):
            if ind[m] < np.power(10, n + 1):
                file_body = output_file + file_index[n] + str(ind[m])
                file_name = file_body + '.tif'
                break
        if precision:
            if axis == 0:
                img = misc.toimage(data[m, :, :], mode='F')
            elif axis == 1:
                img = misc.toimage(data[:, m, :], mode='F')
            elif axis == 2:
                img = misc.toimage(data[:, :, m], mode='F')
        else:
            if axis == 0:
                img = misc.toimage(data[m, :, :])
            elif axis == 1:
                img = misc.toimage(data[:, m, :])
            elif axis == 2:
                img = misc.toimage(data[:, :, m])

        # check if file exists.
        if os.path.isfile(file_name):
            # genarate new file name.
            indq = 1
            FLAG_SAVE = False
            while not FLAG_SAVE:
                new_file_body = file_body + '-' + str(indq)
                new_file_name = new_file_body + '.tif'
                if not os.path.isfile(new_file_name):
                    img.save(new_file_name)
                    FLAG_SAVE = True
                    file_name = new_file_name
                else:
                    indq += 1
        else:
            img.save(file_name)
