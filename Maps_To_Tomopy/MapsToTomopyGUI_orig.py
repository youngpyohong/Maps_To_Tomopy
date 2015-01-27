#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ZetCode PySide tutorial 

This program creates a skeleton of
a classic GUI application with a menubar,
toolbar, statusbar and a central widget. 

author: Jan Bodnar
website: zetcode.com 
last edited: August 2011
"""

import sys
import tkFileDialog
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy import ndimage,optimize
from PIL import Image
import os
from os.path import isfile, join
import string
import pyqtgraph as pg
from PySide import QtGui, QtCore
import h5py
import tomopy

class Example(QtGui.QMainWindow):

      def __init__(self):
            super(Example, self).__init__()
            self.initUI()
        
      def initUI(self):
              
 
##            textEdit = QtGui.QTextEdit()
##            self.setCentralWidget(textEdit)

            exitAction = QtGui.QAction('Exit', self)
            exitAction.triggered.connect(sys.exit)

            openFileAction =QtGui.QAction('Open File',self)
            openFileAction.triggered.connect(self.openfile)

            openFolderAction = QtGui.QAction('Open Folder', self)
            openFolderAction.triggered.connect(self.openfolder)

            sinogramAction = QtGui.QAction('Sinogram', self)
            sinogramAction.triggered.connect(self.sinogram)

            testAction =QtGui.QAction('test',self)
            testAction.triggered.connect(self.test)

            selectElementAction=QtGui.QAction('Select Element',self)
            selectElementAction.triggered.connect(self.selectElement)

            selectFilesAction=QtGui.QAction('Select Files',self)
            selectFilesAction.triggered.connect(self.selectFiles)

            convertAction = QtGui.QAction('Convert', self)
            convertAction.triggered.connect(self.convert)

            showSinogramAction = QtGui.QAction('Sinogram',self)
            showSinogramAction.triggered.connect(self.showSinogram)

            runReconstructAction = QtGui.QAction("Reconstruction", self)
            runReconstructAction.triggered.connect(self.runReconstruct)

            menubar = self.menuBar()
            fileMenu = menubar.addMenu('&File')
            fileMenu.addAction(openFileAction)
            fileMenu.addAction(openFolderAction)
            fileMenu.addAction(exitAction)

            optionMenu = menubar.addMenu('Convert Option')
            optionMenu.addAction(convertAction)
            
            optionMenu.addAction(selectFilesAction)
            
            

            afterConversionMenu = menubar.addMenu('After Conversion')
            afterConversionMenu.addAction(selectElementAction)
            afterConversionMenu.addAction(showSinogramAction)
            afterConversionMenu.addAction(runReconstructAction)

            toolbar = self.addToolBar('ToolBar')
            toolbar.addAction(exitAction)
            toolbar.addAction(openFileAction)
            toolbar.addAction(openFolderAction)
            toolbar.addAction(runReconstructAction)
            toolbar.addAction(selectElementAction)
            toolbar.addAction(convertAction)
            toolbar.addAction(showSinogramAction)
  
            self.move(400, 150)
            self.setWindowTitle('Maps_To_Tomopy')    
            self.show()
#############################
#! Alignment
      def test(self):
            test


#==========================
      def runReconstruct(self):
            self.recon = QSelect3()
            self.recon.setWindowTitle("Reconstruction Window")
            self.recon.numb=len(self.channelname)
            for j in arange(self.recon.numb):
                  self.recon.combo.addItem(self.channelname[j])
            self.recon.show()
            self.recon.btn.setText("Reconstruction")
            self.recon.btn.clicked.connect(self.reconstruct)
            
      def reconstruct(self):
            self.recon.lbl.setText("Reconstruction is currently running")
            self.d= tomopy.xtomo_dataset(log='debug')
            self.reconelement=self.recon.combo.currentIndex()
            self.d.data=self.data[self.reconelement,:,:,:]
            self.d.data[self.d.data == inf] =1
            self.d.dataset(self.d.data, theta=self.theta*np.pi/180)
            self.d.diagnose_center()
            if self.recon.method.currentIndex()==0:
                  self.d.mlem()
            elif self.recon.method.currentIndex()==1:
                  self.d.gridrec()
            elif self.recon.method.currentIndex()==2:
                  self.d.art()
            
            pg.image(self.d.data_recon)
            self.recon.lbl.setText("Done")
            
      def showSinogram(self):
            self.sino = QSelect2()
            self.sino.setWindowTitle("Sinogram Window")
            self.sino.numb=len(self.channelname)
            for j in arange(self.sino.numb):
                  self.sino.combo.addItem(self.channelname[j])
            self.sino.show()
            self.sino.btn.clicked.connect(self.sinogram)
      def convert(self):
            
            self.convert2array()
            self.x=W2(self)
            self.x.textedit.setText("All the files has been converted")
            self.x.btn.setText("OK")
            self.x.show()
      def selectElement(self):
            self.element=QSelect()
            for i in arange(len(self.channelname)):
                  self.element.button[i].setText(self.channelname[i])
            self.element.setWindowTitle("Select Element")
            self.element.show()

      def selectFiles(self):
            self.filecheck=QSelect()
            for i in arange(len(self.fileNames)):
                  onlyfilename=self.fileNames[i].rfind("/")
                  self.filecheck.button[i].setText(self.fileNames[i][onlyfilename+1:])
                  self.filecheck.button[i].setChecked(True)
            self.filecheck.setWindowTitle("Select files")
            self.filecheck.show()

      def openfile(self):
            try:  
                  fileNametemp = QtGui.QFileDialog.getOpenFileNames(self, "Open File",
                                                                QtCore.QDir.currentPath())
                  global RH2
                  self.fileNames=fileNametemp[0]
                  RH2=self.fileNames
                  self.selectFiles()

            except IndexError:
                  print "no file has been selected"


#================
      def openfolder(self):
            try:
                  folderName = QtGui.QFileDialog.getExistingDirectory(self,"Open Folder",
                                                                   QtCore.QDir.currentPath())
                  global RH,fileName
                  RH=folderName
                  file_name_array = [ f for f in os.listdir(folderName) if isfile(join(folderName,f))]
                  file_name_array = [ f for f in os.listdir(folderName) if string.find(f,"h5")!=-1]
                  file_name_array = [folderName+"/"+f for f in file_name_array]

                  self.directory=1
                  self.fileNames=file_name_array
                  fileName=file_name_array
                  self.selectFiles()
            except IndexError:
                  print "no folder has been selected"
            except OSError:
                  print "no folder has been selected"


      def convert2array(self):
            y=zeros(len(self.fileNames),dtype=bool)
            k=arange(y.shape[0])
            for i in arange(len(self.fileNames)):
                  y[i]=self.filecheck.button[i].isChecked()
      
            self.selectedFiles=[self.fileNames[f] for f in k if y[f]==True]
            
            f=h5py.File(os.path.abspath(self.selectedFiles[0]),"r")
            self.channelname=f["exchange"]["images_names"]

            self.channels,self.y,self.x=f["exchange"]["images"].shape
            self.projections=len(self.selectedFiles)
            self.theta= zeros(self.projections)
            self.data=zeros([self.channels,self.projections,self.y,self.x-2])




            for i in arange(self.projections):
                  file_name = os.path.abspath(self.selectedFiles[i])
                  f = h5py.File(file_name,"r")
                  thetatemp=f["MAPS"]["extra_pvs_as_csv"][99]

                  self.theta[i] = thetatemp[thetatemp.rfind(",")+2:]
                  
                  for j in arange(self.channels):
                        self.data[j,i,:,:]=f["/exchange/images"][j,:,:-2]

            print "worked"
#            yy=self.data
#            yy[isinf(yy)]=1
#            pg.image(yy[37,:,:,:])
      def sinogram(self):
            self.thickness=10
            self.sinoelement=self.sino.combo.currentIndex()
            sinodata=self.data[self.sinoelement,:,:,:]
            sinogram=zeros([sinodata.shape[0]*self.thickness,sinodata.shape[2]],dtype=float32)
            for i in arange(self.projections):
                  for j in arange(self.thickness):
                        sinogram[i*self.thickness+j,:]=sinodata[i,sinodata.shape[1]/2,:]
            sinogram[isinf(sinogram)]=1
            
            pg.image(sinogram)
####==============
class W2(QtGui.QDialog):
      def __init__(self, parent=None):
            super(W2, self).__init__(parent)
            self.textedit=QtGui.QLabel()
            self.textedit.setText("What")
            self.btn = QtGui.QPushButton('Click2')
            self.btn.setText("OK")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.textedit)
            vb.addWidget(self.btn)
            
            self.setLayout(vb)
            self.btn.clicked.connect(self.fireupWindows3)

      def fireupWindows3(self):
            
            self.accept()

####===============

class QSelect(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect, self).__init__()
        
            self.initUI()
        
      def initUI(self):
            names=list()
            for i in arange(45):
                  names.append("")


            self.grid = QtGui.QGridLayout()

            j = 0
            pos = [(0, 0), (1, 0), (2,0), (3, 0),(4,0),(5,0),(6,0),(7,0),(8,0),
                (0, 1), (1, 1), (2,1), (3, 1),(4,1),(5,1),(6,1),(7,1),(8,1),
                (0, 2), (1, 2), (2,2), (3, 2),(4,2),(5,2),(6,2),(7,2),(8,2),
                (0, 3), (1, 3), (2,3), (3, 3),(4,3),(5,3),(6,3),(7,3),(8,3),
                (0, 4), (1, 4), (2,4), (3, 4),(4,4),(5,4),(6,4),(7,4),(8,4)]
            self.button=list()
            for i in names:
                  self.button.append(QtGui.QCheckBox(i))
                  self.grid.addWidget(self.button[j], pos[j][0], pos[j][1])
                  j = j + 1
            self.setLayout(self.grid)   
        
            self.move(300, 370)
            self.setWindowTitle('Calculator')    
            self.show()
####==================
class QSelect2(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect2, self).__init__()
        
            self.initUI()
        
      def initUI(self):
 
                  
            self.combo = QtGui.QComboBox(self)
            self.btn = QtGui.QPushButton('Click2')
            self.btn.setText("Sinogram")
            self.lbl=QtGui.QLabel()
            self.lbl.setText("")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.btn)
            vb.addWidget(self.lbl)
            self.setLayout(vb)
####===================
class QSelect3(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect3, self).__init__()
        
            self.initUI()
        
      def initUI(self):
 
                  
            self.combo = QtGui.QComboBox(self)
            self.method = QtGui.QComboBox(self)
            self.btn = QtGui.QPushButton('Click2')
            self.btn.setText("Sinogram")
            self.lbl=QtGui.QLabel()
            self.lbl.setText("")
            self.methodname=["mlem", "gridrec", "art"]
            for k in arange(3):
                  self.method.addItem(self.methodname[k])
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.method)
            vb.addWidget(self.btn)
            vb.addWidget(self.lbl)
            self.setLayout(vb)


def main():
    
      app = QtGui.QApplication(sys.argv)
      ex = Example()
      sys.exit(app.exec_())

def main2():
      img=np.asarray(Image.open("/Users/youngpyohong/Documents/Work/Python/"
                                +"autocorrection_sinogram/mod_sino.tif"),dtype=np.float32)     
      app = QtGui.QApplication(sys.argv)
      ex = Example2("why")
      sys.exit(app.exec_())

if __name__ == '__main__':
    main()
