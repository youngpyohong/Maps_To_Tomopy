# -*- coding: utf-8 -*-
import tomopy
import string
import h5py
import numpy as np
from scipy import *
from pylab import *
import Maps_To_Tomopy
import Maps_To_Tomopy.readwritetxt
from Maps_To_Tomopy.readwritetxt import readtxt, writetxt
import Maps_To_Tomopy.sinogram
from Maps_To_Tomopy.sinogram import sinogram
from PIL import Image
import os
from os import listdir
from os.path import isfile, join
import Maps_To_Tomopy.xtomo_io_f
from Maps_To_Tomopy.xtomo_io_f import xtomo_reader_f, xtomo_writer_f, xtomo_reader_gui
import sys
import tkFileDialog
from pylab import *
from pylab import *
from PySide import QtGui, QtCore

# Read HDF5 file.
#data, white, dark, theta = tomopy.xtomo_reader('demo/demo3.h5',
#                                               slices_start=0,
#                                               slices_end=1)


class Example(QtGui.QMainWindow):
    
      def __init__(self):
            super(Example, self).__init__()
            self.fn=" "
            self.initUI()
        
      def initUI(self):
            self.cb1 = QtGui.QCheckBox('Element 1', self)
            self.cb1.move(20, 20)
            self.cb1.toggle()
            self.cb1.stateChanged.connect(self.element1)
            dir(self.cb1)

            cb2 = QtGui.QCheckBox('Element 2', self)
            cb2.move(20, 40)
            cb2.toggle()
            cb2.stateChanged.connect(self.element2)

           



            exitAction = QtGui.QAction('Exit', self)
            exitAction.setShortcut('Ctrl+Q')
            exitAction.setStatusTip('Exit application')
            exitAction.triggered.connect(self.close)

            directoryAction = QtGui.QAction('Choose Directory',self)
            directoryAction.triggered.connect(self.directory)

            runAction =QtGui.QAction('Run',self)
            runAction.triggered.connect(self.run)

            self.statusBar()




            toolbar = self.addToolBar('Exit')
            toolbar.addAction(exitAction)
            toolbar.addAction(runAction)
            toolbar.addAction(directoryAction)


            self.setGeometry(300, 300, 350, 250)
            self.setWindowTitle('Main window')    
            self.show()
      def directory(self):

            QtGui.QFileDialog.AnyFile
            directory=QtGui.QFileDialog.getExistingDirectory()
            onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
            onlyfiles = [ f for f in onlyfiles if string.find(f,"h5")!=-1]

            for i in arange(len(onlyfiles)):
                  onlyfiles[i]=directory+"/"+onlyfiles[i]
            self.h5files=onlyfiles

      def run(self):
              
            data_orig,theta,channelname=xtomo_reader_gui(self.h5files)

            writetxt(channelname)

            projections=len(self.h5files)
            diff=array([ -4,   5,   8,   4,   3,  -3, -10, -10,  -7,   3,   1,   6,   6,
                  4,   0,   1,   3,   0,  -9], dtype=int32)

        #diff=diff[int(projections_start)-172:int(projections_end)-172+1]


            element=zeros(1)
            element[0]=7
            count=element.shape[0]

            for numb in arange(count):
                  elem=element[numb]
                  data=data_orig[elem,:,:,:]

            
            for i in arange(projections):
                  olddata=data[i,:,:]

                  data[i,:,:]=np.roll(olddata[:,:],diff[i],axis=1)

            sino=sinogram(data)
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
            d.dataset(data, theta=theta*np.pi/180)
  

            d.center=959
            print reconstructiontype, channelname[elem]
            if reconstructiontype=="mlem":
                  d.mlem()
            if reconstructiontype=="art":
                  d.art()
            if reconstructiontype=="gridrec":
                  d.gridrec()

            xtomo_writer_f(d.data_recon, 'tmp/'+str(channelname[elem])+'/'
                              +reconstructiontype+'/test_'+str(d.center), axis=0)
      def element1(self, state):
            if state == QtCore.Qt.Checked:

                  print dir(self.cb1)
                  self.cb1.setHidden(False)
                  self.cb1.setWindowTitle("rhyme")
                  print self.cb1.windowTitle()
                  self.cb1.move(60,60)
                  self.cb1.setAccessibleName(self,"rhyme")

                  self.cb3=QtGui.QCheckBox("Element 3",self)
                  self.cb3.move(20,20)
                  self.cb3.stateChanged.connect(self.element3)
                  print help(self.cb1.setWindowTitle)
            else:
                  self.setWindowTitle('')
      def element2(self,state):
            if state == QtCore.Qt.Checked:
                  self.setWindowTitle('2')
 
            else:
                  self.setWindowTitle('')
      def element3(self,state):
            print "hip hop"
      
############################################
class Example(QtGui.QMainWindow):
    
      def __init__(self):
            super(Example, self).__init__()
            self.fn=" "
            self.initUI()
        
      def initUI(self):
            self.cb1 = QtGui.QCheckBox('Element 1', self)
            self.cb1.move(20, 20)
            self.cb1.toggle()
            self.cb1.stateChanged.connect(self.element1)
            dir(self.cb1)

            cb2 = QtGui.QCheckBox('Element 2', self)
            cb2.move(20, 40)
            cb2.toggle()
            cb2.stateChanged.connect(self.element2)

           



            exitAction = QtGui.QAction('Exit', self)
            exitAction.setShortcut('Ctrl+Q')
            exitAction.setStatusTip('Exit application')
            exitAction.triggered.connect(self.close)

            directoryAction = QtGui.QAction('Choose Directory',self)
            directoryAction.triggered.connect(self.directory)

            runAction =QtGui.QAction('Run',self)
            runAction.triggered.connect(self.run)

            self.statusBar()




            toolbar = self.addToolBar('Exit')
            toolbar.addAction(exitAction)
            toolbar.addAction(runAction)
            toolbar.addAction(directoryAction)


            self.setGeometry(300, 300, 350, 250)
            self.setWindowTitle('Main window')    
            self.show()
      def directory(self):

            QtGui.QFileDialog.AnyFile
            directory=QtGui.QFileDialog.getExistingDirectory()
            onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
            onlyfiles = [ f for f in onlyfiles if string.find(f,"h5")!=-1]

            for i in arange(len(onlyfiles)):
                  onlyfiles[i]=directory+"/"+onlyfiles[i]
            self.h5files=onlyfiles

      def run(self):
              
            data_orig,theta,channelname=xtomo_reader_gui(self.h5files)

            writetxt(channelname)

            projections=len(self.h5files)
            diff=array([ -4,   5,   8,   4,   3,  -3, -10, -10,  -7,   3,   1,   6,   6,
                  4,   0,   1,   3,   0,  -9], dtype=int32)

        #diff=diff[int(projections_start)-172:int(projections_end)-172+1]


            element=zeros(1)
            element[0]=7
            count=element.shape[0]

            for numb in arange(count):
                  elem=element[numb]
                  data=data_orig[elem,:,:,:]

            
            for i in arange(projections):
                  olddata=data[i,:,:]

                  data[i,:,:]=np.roll(olddata[:,:],diff[i],axis=1)

            sino=sinogram(data)
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
            d.dataset(data, theta=theta*np.pi/180)
  

            d.center=959
            print reconstructiontype, channelname[elem]
            if reconstructiontype=="mlem":
                  d.mlem()
            if reconstructiontype=="art":
                  d.art()
            if reconstructiontype=="gridrec":
                  d.gridrec()

            xtomo_writer_f(d.data_recon, 'tmp/'+str(channelname[elem])+'/'
                              +reconstructiontype+'/test_'+str(d.center), axis=0)
      def element1(self, state):
            if state == QtCore.Qt.Checked:

                  print dir(self.cb1)
                  self.cb1.setHidden(False)
                  self.cb1.setWindowTitle("rhyme")
                  print self.cb1.windowTitle()
                  self.cb1.move(60,60)
                  self.cb1.setObjectName("rhyme")
                  print self.cb1.objectName()
                  self.cb3=QtGui.QCheckBox("Element 3",self)
                  self.cb3.move(20,20)
                  self.cb3.stateChanged.connect(self.element3)
                  print help(self.cb1.setWindowTitle)
            else:
                  self.setWindowTitle('')
      def element2(self,state):
            if state == QtCore.Qt.Checked:
                  self.setWindowTitle('2')
 
            else:
                  self.setWindowTitle('')
      def element3(self,state):
            print "hip hop"
####################################### 
def main():
    
      app = QtGui.QApplication(sys.argv)
      ex = Example()
      sys.exit(app.exec_())

if __name__ == '__main__':
      main()

a,b,c=xtomo_reader_gui(onlyfiles)
print here
####################


