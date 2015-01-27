#!/usr/bin/python
# -*- coding: utf-8 -*-


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
            self.ImageTag = "exchange"

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

            selectImageTagAction = QtGui.QAction("Select Image Tag", self)
            selectImageTagAction.triggered.connect(self.selectImageTag)

            xCorAction = QtGui.QAction("Cross Correlation", self)
            xCorAction.triggered.connect(self.CrossCorrelation_test)

            alignFromTextAction = QtGui.QAction("Alignment from Text", self)
            alignFromTextAction.triggered.connect(self.alignFromText)

            restoreAction = QtGui.QAction("Restore", self)
            restoreAction.triggered.connect(self.restore)

            menubar = self.menuBar()
            fileMenu = menubar.addMenu('&File')
            fileMenu.addAction(openFileAction)
            fileMenu.addAction(openFolderAction)
            fileMenu.addAction(exitAction)

            optionMenu = menubar.addMenu('Convert Option')
            optionMenu.addAction(convertAction)
            
            optionMenu.addAction(selectFilesAction)
            optionMenu.addAction(selectImageTagAction)
            
            alignmentMenu = menubar.addMenu("Alignment")
            alignmentMenu.addAction(xCorAction)
            alignmentMenu.addAction(alignFromTextAction)
            alignmentMenu.addAction(restoreAction)
            
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

      def CrossCorrelation_test(self):
            
            self.datacopy=zeros(self.data.shape)
            self.datacopy[...]=self.data[...]
            self.data[np.isnan(self.data)]=1
            self.xcor = QSelect3()
            self.xcor.setWindowTitle("CrossCorrelation Window")
            self.xcor.numb = len(self.channelname)
            for j in arange(self.xcor.numb):
                  self.xcor.combo.addItem(self.channelname[j])
            self.xcor.btn.setText("Cross Correlation")
            print dir(self.xcor.method)
            self.xcor.method.setVisible(False)
            self.xcor.save.setVisible(True)
            self.xcor.btn.clicked.connect(self.xCor)
            self.xcor.save.clicked.connect(self.restore)
            self.xcor.show()
      def restore(self):
            self.data=zeros(self.datacopy.shape)
            self.data[...]=self.datacopy[...]
      def prexCor(self):
            try:
                  self.xcor.savedir=QtGui.QFileDialog.getSaveFileName()
                  if not self.xcor.savedir[0]:
                        raise IndexError
                  self.xcor.savedir=self.xcor.savedir[0]
                  self.xCor()
            except IndexError:
                  print "type the header name"

      
      def xCor(self):
            self.xcor.savedir="texting"
            f=open(self.xcor.savedir+".txt",'w')
            onlyfilename=self.fileNames[0].rfind("/")
            f.write(self.fileNames[0][onlyfilename+1:]+" \n")
            f.write("0 \n")
            self.xcorElement = self.xcor.combo.currentIndex()
            for i in arange(self.projections-1):
                  onlyfilename=self.fileNames[i+1].rfind("/")
                  img1=self.data[self.xcorElement,i,:,:]
                  img2=self.data[self.xcorElement,i+1,:,:]
                  xcorimg=self.xcorrelate(img1,img2)
                  print np.where(xcorimg==np.max(xcorimg))
                  print xcorimg
                  xcorimg[np.isnan(xcorimg)]=1
                  xcorindex=np.where(xcorimg==np.max(xcorimg))[1][0]
                  shift=xcorindex-(xcorimg.shape[1])/2
                  print xcorindex, xcorimg.shape[1]
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],shift,axis=2)
                  f.write(self.fileNames[i+1][onlyfilename+1:]+" \n")
                  f.write(str(shift)+ "\n")
                  print i
            f.close()

      def alignFromText(self):
            try:
                  fileName = QtGui.QFileDialog.getOpenFileName(self, "Open File",
                        QtCore.QDir.currentPath(),"TXT (*.txt)")
                  ##### for future reference "All File (*);;CSV (*.csv *.CSV)"
                  f=open(fileName[0],'r')
                  read = f.readlines()
                  self.datacopy=zeros(self.data.shape)
                  self.datacopy[...]=self.data[...]
                  self.data[np.isnan(self.data)]=1
                  for i in arange(self.projections):
                        onlyfilename=self.fileNames[0].rfind("/")
                        for j in arange(len(read)):
                              if string.find(read[j],self.fileNames[0][onlyfilename+1:])!=-1:
                                    shift=int(read[j+1][:-1])
                                    self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],shift,axis=2)

                  
            except IOError:
                  print "choose file please"
                  
            
      
            


#==========================
      def xcorrelate(self,image1,image2,edgeguass_sigma=4):

            image1_ft=fftshift(fft2(edgegauss(image1)))
            image2_ft=fftshift(fft2(edgegauss(image2)))


            xcor_ft = image1_ft*np.conjugate(image2_ft)
            xcor_image=abs(fftshift(fft(fftshift(xcor_ft))))
            return xcor_image

      
      def edgegauss(self,imagey,sigma=4):
            image=zeros(imagey.shape)
            image[...]=imagey[...]
            nx=image.shape[1]
            ny=image.shape[0]

            n_sigma = -log(10**-6)
            n_roll = max(int(1+sigma*n_sigma),2)
            exparg = float32(arange(n_roll)/float(sigma))
            rolloff = float32(1)-exp(-0.5*exparg*exparg)

            ## Top edge

            xstart = 0
            xstop = nx 
            iy= 0

            for i_roll in arange(n_roll):
                  image[iy,xstart:xstop] = image[iy, xstart:xstop] * rolloff[iy]
                  xstart = min(xstart+1,nx/2-1)
                  xstop = max(xstop - 1,nx/2)
                  iy = min(iy+1,ny-1)


            ## Bottom edge

            xstart = 0
            xstop = nx
            iy = ny-1

            for i_roll in arange(n_roll):
                  image[iy,xstart:xstop] = image[iy,xstart:xstop]*rolloff[ny-1-iy]
                  xstart = min(xstart +1,nx/2-1)
                  xstop = max(xstop-1,nx/2)
                  iy = max(iy-1, 0)

            ## Left edge

            ystart = 1
            ystop = ny - 1
            ix = 0

            for i_roll in arange(n_roll):
                  image[ystart:ystop,ix]=image[ystart:ystop,ix] * rolloff[ix]
                  ystart = min(ystart+1,ny/2-1)
                  ystop = max(ystop-1 , ny/2)
                  ix = min(ix+1,nx-1)

            ## Right edge

            ystart = 1
            ystop = ny-1
            ix = nx-1

            for i_roll in arange(n_roll):
                  image[ystart:ystop, ix] = image[ystart:ystop,ix] * rolloff[nx-1-ix]
                  ystart = min(ystart+1, ny/2-1)
                  ystop = max(ystop-1, ny/2)
                  ix = max(ix - 1,0)
                  
            return image
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
            self.recon.save.clicked.connect(self.saveRecTiff)
            
      def reconstruct(self):
            self.recon.lbl.setText("Reconstruction is currently running")
            self.d= tomopy.xtomo_dataset(log='debug')
            self.reconelement=self.recon.combo.currentIndex()
            self.d.data=self.data[self.reconelement,:,:,:]
            self.d.data[self.d.data == inf] =2
            self.d.data[np.isnan(data)]=2
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
##            print dir(self.recon.save)
            self.recon.save.setHidden(False)

      def saveRecTiff(self):
            try:
                  self.savedir=QtGui.QFileDialog.getSaveFileName()
                  if not self.savedir[0]:
                        raise IndexError
                  self.savedir=self.savedir[0]
                  tomopy.xtomo_writer(self.d.data_recon,self.savedir,axis=0,digits=4)
            except IndexError:
                  print "type the header name"

      def selectImageTag(self):
            self.sit = QSelect3()
            self.sit.setWindowTitle("Seletect Image Tag from h5 file")
            self.sit.data=h5py.File(self.fileNames[0])
            self.sit.firstColumn = self.sit.data.items()
            self.sit.firstColumnNum = len(self.sit.firstColumn)
            for i in arange(self.sit.firstColumnNum):
                  self.sit.combo.addItem(self.sit.firstColumn[i][0])
##            print self.sit.combo.currentIndex()
##            self.sit.secondColumnName = self.sit.firstColumn[self.sit.combo.currentIndex()][0]
##            
##            self.sit.secondColumn = self.sit.data[self.sit.secondColumnName].items()
##            self.sit.secondColumnNum = len(self.sit.secondColumn)
##            for j in arange(self.sit.secondColumnNum):
##                  self.sit.combo.addItem(self.sit.secondColumn[j][0])
            self.sit.btn.setText("Set")
            self.sit.method.setHidden(True)
            self.sit.btn.clicked.connect(self.setImageTag)
            self.sit.show()

      def setImageTag(self):
            self.ImageTag = self.sit.combo.currentText()
            print "Image Tag has been set to \"", self.ImageTag, "\"" 
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
                  if not self.fileNames:
                        raise IndexError
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
            self.channelname=f[self.ImageTag]["images_names"]

            self.channels,self.y,self.x=f[self.ImageTag]["images"].shape
            self.projections=len(self.selectedFiles)
            self.theta= zeros(self.projections)
            self.data=zeros([self.channels,self.projections,self.y,self.x-2])




            for i in arange(self.projections):
                  file_name = os.path.abspath(self.selectedFiles[i])
                  f = h5py.File(file_name,"r")
                  thetatemp=f["MAPS"]["extra_pvs_as_csv"][8]

                  self.theta[i] = thetatemp[thetatemp.rfind(",")+2:]
                  
                  for j in arange(self.channels):
                        self.data[j,i,:,:]=f[self.ImageTag]["images"][j,:,:-2]

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
                  sinogram[i*self.thickness:(i+1)*self.thickness,:]=sinodata[i,sinodata.shape[1]/2-5:sinodata.shape[1]/2+5,:]
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
            for i in arange(100):
                  names.append("")


            self.grid = QtGui.QGridLayout()

            j = 0
##            pos = [(0, 0), (1, 0), (2,0), (3, 0),(4,0),(5,0),(6,0),(7,0),(8,0),
##                (0, 1), (1, 1), (2,1), (3, 1),(4,1),(5,1),(6,1),(7,1),(8,1),
##                (0, 2), (1, 2), (2,2), (3, 2),(4,2),(5,2),(6,2),(7,2),(8,2),
##                (0, 3), (1, 3), (2,3), (3, 3),(4,3),(5,3),(6,3),(7,3),(8,3),
##                (0, 4), (1, 4), (2,4), (3, 4),(4,4),(5,4),(6,4),(7,4),(8,4)]
            pos = list()
            for y in arange(10):
                  for x in arange(10):
                        pos.append((x,y))
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
            self.save = QtGui.QPushButton("Save tiff files")
            self.save.setHidden(True)
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
            vb.addWidget(self.save)
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
