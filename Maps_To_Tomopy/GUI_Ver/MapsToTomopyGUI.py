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
            self.thetaPos = 99


            exitAction = QtGui.QAction('Exit', self)
            exitAction.triggered.connect(sys.exit)

            openFileAction =QtGui.QAction('Open File',self)
            openFileAction.triggered.connect(self.openfile)

            openFolderAction = QtGui.QAction('Open Folder', self)
            openFolderAction.triggered.connect(self.openfolder)

            sinogramAction = QtGui.QAction('Sinogram', self)
            sinogramAction.triggered.connect(self.sinogram)

            saveImageAction = QtGui.QAction('Save Projections', self)
            saveImageAction.triggered.connect(self.saveImage)

            testAction =QtGui.QAction('test',self)
            testAction.triggered.connect(self.test)

            selectElementAction=QtGui.QAction('Select Element',self)
            selectElementAction.triggered.connect(self.selectElement)

            selectFilesAction=QtGui.QAction('Select Files',self)
            selectFilesAction.triggered.connect(self.selectFiles)

            convertAction = QtGui.QAction('Save data in memory', self)
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

            readConfigAction = QtGui.QAction("Read configuration file",self)
            readConfigAction.triggered.connect(self.readConfigFile)



            self.lbl=QtGui.QLabel()
            self.setCentralWidget(self.lbl)
            self.lbl.setText("Starting")

            menubar = self.menuBar()
            fileMenu = menubar.addMenu('&File')
            fileMenu.addAction(openFileAction)
            fileMenu.addAction(openFolderAction)
            fileMenu.addAction(readConfigAction)
            fileMenu.addAction(exitAction)

            self.optionMenu = menubar.addMenu('Convert Option')
            self.optionMenu.addAction(convertAction)
            
            self.optionMenu.addAction(selectFilesAction)
            self.optionMenu.addAction(selectImageTagAction)
            self.optionMenu.setDisabled(True)

            self.alignmentMenu = menubar.addMenu("Alignment")
            self.alignmentMenu.addAction(xCorAction)
            self.alignmentMenu.addAction(alignFromTextAction)
            self.alignmentMenu.addAction(restoreAction)


            self.afterConversionMenu = menubar.addMenu('After saving data in memory')
            self.afterConversionMenu.addAction(saveImageAction)
            self.afterConversionMenu.addAction(selectElementAction)
            self.afterConversionMenu.addAction(showSinogramAction)
            self.afterConversionMenu.addAction(runReconstructAction)
            self.afterConversionMenu.setDisabled(True)

            toolbar = self.addToolBar('ToolBar')
            toolbar.addAction(exitAction)
            toolbar.addAction(openFileAction)
            toolbar.addAction(openFolderAction)
            toolbar.addAction(saveImageAction)
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
            self.d.data[self.d.data==inf]=0.01
            self.d.data[np.isnan(self.d.data)]=0.01

            
            ###TEMP
            #self.d.data = (np.exp(-0.0001*self.d.data)).astype('float32')
            #self.d.center=128
            ###TEMP

            if self.recon.method.currentIndex()==0:
                  self.d.dataset(self.d.data, theta=self.theta*np.pi/180)
                  #self.d.optimize_center()
                  self.d.mlem()
            elif self.recon.method.currentIndex()==1:
                  self.d.dataset(self.d.data, theta=self.theta)
                  #self.d.optimize_center()
                  self.d.gridrec()
            elif self.recon.method.currentIndex()==2:
                  self.d.dataset(self.d.data, theta=self.theta*np.pi/180)
                  #self.d.optimize_center()
                  self.d.art()
            print self.d.center
            pg.image(self.d.data_recon)
            self.recon.lbl.setText("Done")
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
#=============================
                  
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
            self.sit.combo.setCurrentIndex(i)
            self.sit.btn.setText("Set")
            self.sit.method.setHidden(True)
            self.sit.btn.clicked.connect(self.setImageTag)
            self.sit.show()

      def setImageTag(self):
            self.ImageTag = self.sit.combo.currentText()
            self.lbl.setText("Image Tag has been set to \""+self.ImageTag+"\"")
            print "Image Tag has been set to \"", self.ImageTag, "\""
#==============================
      def showSinogram(self):
            self.sino = QSelect2()
            self.sino.setWindowTitle("Sinogram Window")
            self.sino.numb=len(self.channelname)
            for j in arange(self.sino.numb):
                  self.sino.combo.addItem(self.channelname[j])
            self.sino.show()
            self.sino.btn.clicked.connect(self.sinogram)

            
      def convert(self):
            self.afterConversionMenu.setDisabled(False)
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
                  f = h5py.File(os.path.abspath(self.fileNames[i]),"r")
                  thetatemp = f["MAPS"]["extra_pvs_as_csv"][8]
                  thetapos = string.rfind(thetatemp, ",")
                  theta = str(round(float(thetatemp[thetapos+1:])))
                  onlyfilename=self.fileNames[i].rfind("_")
                  self.filecheck.button[i].setText(self.fileNames[i][onlyfilename+1:]+" ("+theta+")")
                  self.filecheck.button[i].setChecked(True)
            self.ImageTag=f.items()[-1][0]
            self.lbl.setText("Image Tag has been set to \""+self.ImageTag+"\"")
            self.filecheck.setWindowTitle("Select files")
            self.filecheck.show()

            self.optionMenu.setEnabled(True)
            self.filecheck.btn.clicked.connect(self.convert)
            self.filecheck.btn2.clicked.connect(self.selectImageTag)

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
            self.data=zeros([self.channels,self.projections,self.y,self.x])




            for i in arange(self.projections):
                  file_name = os.path.abspath(self.selectedFiles[i])
                  f = h5py.File(file_name,"r")
                  thetatemp=f["MAPS"]["extra_pvs_as_csv"][self.thetaPos]

                  self.theta[i] = float(thetatemp[thetatemp.rfind(",")+3:])
                  
                  for j in arange(self.channels):
                        self.data[j,i,:,:]=f[self.ImageTag]["images"][j,:,:]
                  print i, "projection(s) has/have been converted"
            print "worked"
#            yy=self.data
#            yy[isinf(yy)]=1
#            pg.image(yy[37,:,:,:])


#####!!!! just temp need to be fixed
      def saveImage(self):
            for j in arange(self.data.shape[0]):
                  path=os.getcwd()+"/images/"+self.channelname[j]
                  try:
                        os.mkdir(path)
                  except OSError:
                        path=path
                  for i in arange(self.data.shape[1]):
                        temp_img=self.data[j,i,:,:]
                        temp=Image.fromarray(temp_img.astype(np.float32))
                        index=string.rfind(self.fileNames[i],"/")
                        temp.save(path+self.fileNames[i][index:-3]+"_"+self.channelname[j]+".tif")
#####!!!!

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
      def readConfigFile(self):
            try:
                  fileName = QtGui.QFileDialog.getOpenFileName(self, "Select config file",
                        QtCore.QDir.currentPath(),"TXT (*.txt)")
                  ##### for future reference "All File (*);;CSV (*.csv *.CSV)"
                  f=open(fileName[0],'r')
                  l = f.readlines()
                  for i in arange(len(l)):
                          ## image_tag
                        if string.find(l[i],"Image_Tag")!=-1:
                              self.ImageTag=l[i+1][:-1]
                          ## theta postion from extra pvs as csv
                        if string.find(l[i],"Theta position")!=-1:
                              print l[i+1][:-1]
                              self.thetaPos=int(l[i+1][:-1])
                                    
            except IOError:
                  print "choose file please"



            
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
            self.lbl=QtGui.QLabel()
            self.lbl2=QtGui.QLabel()
            self.lbl.setText("closing this window won't affect your selection of the files")
            self.lbl2.setText("You should convert the files in order to generate sinogram or reconstructed data")
            self.btn=QtGui.QPushButton()
            self.btn.setText("Convert")
            self.btn2=QtGui.QPushButton()
            self.btn2.setText("set Image Tag")

            j = 0
            pos=list()
            for y in arange(10):
                  for x in arange(10):
                        pos.append((x,y))

            self.button=list()
            for i in names:
                  self.button.append(QtGui.QCheckBox(i))
                  self.grid.addWidget(self.button[j], pos[j][0], pos[j][1])
                  j = j + 1
            self.setLayout(self.grid)

            self.vb=QtGui.QVBoxLayout()
            self.vb2=QtGui.QVBoxLayout()
            
            self.vb.addWidget(self.lbl,11,0)
            self.vb.addWidget(self.lbl2,12,0)

            self.vb2.addWidget(self.btn2,13,0)
            self.vb2.addWidget(self.btn,14,0)

            self.grid.addLayout(self.vb,11,0,1,10)
            self.grid.addLayout(self.vb2,13,3,1,2)
        
            self.move(100, 100)
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
