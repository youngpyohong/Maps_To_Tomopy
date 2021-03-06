# -*- coding: utf-8 -*-
#!/usr/bin/python


import sys
import tkFileDialog
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy import ndimage,optimize
import scipy.ndimage.interpolation as spni
from skimage.feature import match_template
from PIL import Image
import os
from os.path import isfile, join
import string
import pyqtgraph as pg
from PySide import QtGui, QtCore
from pyqtgraph import QtGui, QtCore
import h5py
import tomopy
import subpixelshift
from subpixelshift import *


class Example(QtGui.QMainWindow):

      def __init__(self):
            super(Example, self).__init__()
            self.initUI()
        
      def initUI(self):
              
 
##            textEdit = QtGui.QTextEdit()
##            self.setCentralWidget(textEdit)
            
            self.ImageTag = "exchange"
            self.thetaPos = 657


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

            saveAlignToTextAction = QtGui.QAction("Save Alignment information to text", self)
            saveAlignToTextAction.triggered.connect(self.saveAlignToText)

            restoreAction = QtGui.QAction("Restore", self)
            restoreAction.triggered.connect(self.restore)

            readConfigAction = QtGui.QAction("Read configuration file",self)
            readConfigAction.triggered.connect(self.readConfigFile)

            runCenterOfMassAction = QtGui.QAction("run center of mass action",self)
            runCenterOfMassAction.triggered.connect(self.centerOfMassWindow)

            alignCenterOfMassAction = QtGui.QAction("Align by fitting center of mass position into sine curve",self)
            alignCenterOfMassAction.triggered.connect(self.alignCenterOfMass)

            matcherAction = QtGui.QAction("match template", self)
            matcherAction.triggered.connect(self.match_window)

            exportDataAction = QtGui.QAction("export data",self)
            exportDataAction.triggered.connect(self.export_data)

            runTransRecAction = QtGui.QAction("Transmission Recon", self)
            runTransRecAction.triggered.connect(self.runTransReconstruct)

            saveHotSpotPosAction = QtGui.QAction("Save Hot Spot Pos", self)
            saveHotSpotPosAction.triggered.connect(self.saveHotSpotPos)

            alignHotSpotPosAction = QtGui.QAction("Align Hot Spot pos", self)
            alignHotSpotPosAction.triggered.connect(self.alignHotSpotPos1)

            ###
            self.frame = QtGui.QFrame()
            self.vl = QtGui.QHBoxLayout()


            
            self.tab_widget = QtGui.QTabWidget()
            self.tab_widget.addTab(self.createMessageWidget(), unicode("Message"))
            self.tab_widget.addTab(self.createProjWidget(), unicode("Projections"))
            self.tab_widget.addTab(self.createSinoWidget(), unicode("Sinogram"))
            self.tab_widget.addTab(self.createReconWidget(), unicode("Reconstruction"))
            #self.tab_widget.addTab(self.sinoGroup, unicode("Sinogram"))
            self.tab_widget.addTab(self.createSaveHotspotWidget(),unicode("Save Hotspot"))

            self.vl.addWidget(self.tab_widget)

##            self.win=pg.ImageView()
##            self.layout1 = QtGui.QVBoxLayout()
##            self.layout1.addWidget(self.win)
##            self.vl.addWidget(self.layout1)
            self.frame.setLayout(self.vl)
            self.setCentralWidget(self.frame)
            self.tab_widget.setDisabled(True)

            ###
            
##            self.lbl=QtGui.QLabel()
##            self.setCentralWidget(self.lbl)
##            self.lbl.setText("Starting")

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
            self.alignmentMenu.addAction(saveAlignToTextAction)
            self.alignmentMenu.addAction(runCenterOfMassAction)
            self.alignmentMenu.addAction(alignCenterOfMassAction)
            self.alignmentMenu.addAction(xCorAction)
            self.alignmentMenu.addAction(matcherAction)
            self.alignmentMenu.addAction(alignFromTextAction)
            self.alignmentMenu.addAction(saveHotSpotPosAction)
            self.alignmentMenu.addAction(alignHotSpotPosAction)
            self.alignmentMenu.addAction(restoreAction)
            self.alignmentMenu.setDisabled(True)

            self.afterConversionMenu = menubar.addMenu('After saving data in memory')
            self.afterConversionMenu.addAction(saveImageAction)
            #self.afterConversionMenu.addAction(selectElementAction)
            self.afterConversionMenu.addAction(showSinogramAction)
            self.afterConversionMenu.addAction(runReconstructAction)
            self.afterConversionMenu.setDisabled(True)

            toolbar = self.addToolBar('ToolBar')
            toolbar.addAction(exitAction)
            toolbar.addAction(openFileAction)
            toolbar.addAction(openFolderAction)
            toolbar.addAction(saveHotSpotPosAction)
            toolbar.addAction(alignHotSpotPosAction)
            toolbar.addAction(exportDataAction)
            toolbar.addAction(runTransRecAction)
            toolbar.addAction(runCenterOfMassAction)
            toolbar.addAction(matcherAction)
            toolbar.addAction(runReconstructAction)
            toolbar.addAction(selectElementAction)
            toolbar.addAction(convertAction)
            toolbar.addAction(showSinogramAction)
            toolbar.setVisible(False)
  
            self.setGeometry(0,0, 1000,500)
            self.setWindowTitle('Maps_To_Tomopy')    
            self.show()
#############################
## creating tab
      def createGridInputWidget(self):
          GridStartVal = '2'
          hBox = QtGui.QHBoxLayout()
          self.GridXIn = QtGui.QLineEdit()
          self.GridYIn = QtGui.QLineEdit()
          self.GridZIn = QtGui.QLineEdit()

          self.GridXIn.setText(GridStartVal)
          self.GridYIn.setText(GridStartVal)
          self.GridZIn.setText(GridStartVal)

          hBox.addWidget(QtGui.QLabel("X"))
          hBox.addWidget(self.GridXIn)
          hBox.addWidget(QtGui.QLabel("Y"))
          hBox.addWidget(self.GridYIn)
          hBox.addWidget(QtGui.QLabel("Z"))
          hBox.addWidget(self.GridZIn)

          gridGroup = QtGui.QGroupBox('Grid Size')
          gridGroup.setLayout(hBox)

          return gridGroup
      def createMessageWidget(self):
            GridStartVal = '2'
            hbox = QtGui.QHBoxLayout()
            self.lbl = QtGui.QLineEdit()
            self.lbl.setText("Starting")

            hbox.addWidget(QtGui.QLabel("Message"))
            hbox.addWidget(self.lbl)

            messageGroup = QtGui.QGroupBox("Message Box")
            messageGroup.setLayout(hbox)

            return messageGroup

      def createSinoWidget(self):
            self.sino = QSelect2()
            self.sinoView = IView()

            sinoBox = QtGui.QHBoxLayout()
            sinoBox.addWidget(self.sino)
            sinoBox.addWidget(self.sinoView,10)
            sinoGroup=QtGui.QGroupBox("Sinogram")
            sinoGroup.setLayout(sinoBox)

            return sinoGroup

      def createReconWidget(self):
            self.recon = QSelect3()
            self.reconView = pg.ImageView()

            reconBox = QtGui.QHBoxLayout()
            reconBox.addWidget(self.recon)
            reconBox.addWidget(self.reconView,10)
            reconGroup=QtGui.QGroupBox("Reconstruction")
            reconGroup.setLayout(reconBox)

            return reconGroup

      def createProjWidget(self):
            self.projection = QSelect2()
            self.projectionView =pg.ImageView()

            projectionBox = QtGui.QHBoxLayout()
            projectionBox.addWidget(self.projection)
            projectionBox.addWidget(self.projectionView,10)
            projectionGroup=QtGui.QGroupBox("Projections")
            projectionGroup.setLayout(projectionBox) 

            return projectionGroup

      def createSaveHotspotWidget(self):
            self.projViewControl = QSelect4()
            self.projView = IView()

            projViewBox = QtGui.QHBoxLayout()
            projViewBox.addWidget(self.projViewControl)
            projViewBox.addWidget(self.projView,10)
            projViewGroup=QtGui.QGroupBox("Save HotSpot")
            projViewGroup.setLayout(projViewBox) 
            return projViewGroup
#############################
#! Alignment
      def test(self):
            test

      def centerOfMassWindow(self):
            self.comer = QSelect3()
            self.comer.setWindowTitle("Center of Mass window")
            self.comer.numb= len(self.channelname)
            for j in arange(self.comer.numb):
                  self.comer.combo.addItem(self.channelname[j])
            self.comer.btn.setText("Center of Mass")
            self.comer.method.setVisible(False)
            self.comer.save.setVisible(True)
            self.comer.save.setText("Restore")
            self.comer.btn.clicked.connect(self.runCenterOfMass)
            self.comer.save.clicked.connect(self.restore)
            self.comer.show()

      def runCenterOfMass(self):
            self.centerOfMass()
            self.fitCenterOfMass()
            self.lbl.setText("Center of Mass: "+str(self.p1[2]))
#### devel
      def runCenterOfMass2(self):
            self.com = zeros(self.projections)
            temp=zeros(self.data.shape[3])
            temp2=zeros(self.data.shape[3])
            self.comelem = self.sino.combo.currentIndex()
            for i in arange(self.projections):
                  temp=sum(self.data[self.comelem,i,self.sino.sld.value()-self.thickness/2:self.sino.sld.value()+self.thickness/2,:]-self.data[self.comelem,i,:10,:10].mean(), axis=0)
                  #temp=sum(self.data[self.comelem,i,:,:]-1, axis=0)
                  numb2=sum(temp)
                  for j in arange(self.data.shape[3]):
                        
                        temp2[j]=temp[j]*j
                  numb=float(sum(temp2))/numb2
                  self.com[i]=numb
            self.fitCenterOfMass()
            self.lbl.setText("Center of Mass: "+str(self.p1[2]))
            self.alignCenterOfMass()
            self.sinogram()
            
#### devel
            
      def centerOfMass(self):
            self.com = zeros(self.projections)
            temp=zeros(self.data.shape[3])
            temp2=zeros(self.data.shape[3])
            self.comelem = self.comer.combo.currentIndex()
            for i in arange(self.projections):
                  temp=sum(self.data[self.comelem,i,:,:]-self.data[self.comelem,i,:10,:10].mean(), axis=0)
                  #temp=sum(self.data[self.comelem,i,:,:]-1, axis=0)
                  numb2=sum(temp)
                  for j in arange(self.data.shape[3]):
                        
                        temp2[j]=temp[j]*j
                  numb=float(sum(temp2))/numb2
                  self.com[i]=numb

            
      def fitCenterOfMass(self):
            x=self.theta
            self.fitfunc = lambda p,x: p[0]*sin(2*pi/360*(x-p[1]))+p[2]
            self.errfunc = lambda p,x,y: self.fitfunc(p,x)-y
            p0=[100,100,100]
            self.p1,success = optimize.leastsq(self.errfunc,p0[:],args=(x,self.com))
            self.centerOfMassDiff=self.fitfunc(self.p1,x)-self.com

      def alignCenterOfMass(self):
            for i in arange(self.projections):
                  self.xshift[i]=int(self.centerOfMassDiff[i])
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], self.xshift[i], axis=2)
            self.lbl.setText("Alignment has been completed")
                  

                  
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
            self.data=zeros(self.oldData.shape)
            self.data[...]=self.oldData[...]
      def prexCor(self):
            try:
                  self.xcor.savedir=QtGui.QFileDialog.getSaveFileName()
                  if not self.xcor.savedir[0]:
                        raise IndexError
                  self.xcor.savedir=self.xcor.savedir[0]
                  self.xCor()
            except IndexError:
                  print "type the header name"
      def match_window(self):
            self.matcher = QSelect3()
            self.matcher.setWindowTitle("Match template window")
            self.matcher.numb= len(self.channelname)
            for j in arange(self.matcher.numb):
                  self.matcher.combo.addItem(self.channelname[j])
            self.matcher.btn.setText("Match Template")
            self.matcher.method.setVisible(False)
            self.matcher.save.setVisible(True)
            self.matcher.save.setText("Restore")
            self.matcher.btn.clicked.connect(self.match)
            self.matcher.save.clicked.connect(self.restore)
            self.matcher.show()
      def match(self):
            self.matchElem=self.matcher.combo.currentIndex()
            for i in arange(self.projections-1):
                  img=self.data[self.matchElem,i,:,:]
                  img1=ones([img.shape[0]*2,img.shape[1]*2])*self.data[self.matchElem,i,:10,:10].mean()
                  img1[img.shape[0]/2:img.shape[0]*3/2,img.shape[1]/2:img.shape[1]*3/2]=img
                  img2=self.data[self.matchElem,i+1,:,:]
                  result=match_template(img1,img2)
                  result=np.where(result==np.max(result))
                  self.yshift[i+1]=result[0][0]-img.shape[0]/2
                  self.xshift[i+1]=result[1][0]-img.shape[1]/2
                  print self.xshift[i+1], self.yshift[i+1]
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], self.xshift[i+1], axis=2)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], self.yshift[i+1], axis=1)

      def xCor(self):
##            self.xcor.savedir="texting"
##            f=open(self.xcor.savedir+".txt",'w')
##            onlyfilename=self.fileNames[0].rfind("/")
##            f.write(self.fileNames[0][onlyfilename+1:]+" \n")
##            f.write("0 \n")
            self.xcorElement = self.xcor.combo.currentIndex()
            for i in arange(self.projections-1):
                  onlyfilename=self.fileNames[i+1].rfind("/")
                  img1=self.data[self.xcorElement,i,:,:]
                  img2=self.data[self.xcorElement,i+1,:,:]
                  
                  self.t0,self.t1=self.xcorrelate(img1,img2)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],self.t0,axis=1)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],self.t1,axis=2)

                  
##                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],shift,axis=2)
##                  f.write(self.fileNames[i+1][onlyfilename+1:]+" \n")
##                  f.write(str(shift)+ "\n")
##                  print i
##            f.close()

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
                        onlyfilename=self.selectedFiles[i].rfind("/")
                        for j in arange(len(read)):
                              if string.find(read[j],self.selectedFiles[i][onlyfilename+1:])!=-1:
                                    secondcol=read[j].rfind(",")
                                    firstcol=read[j][:secondcol].rfind(",")
                                    self.yshift[i]=int(float(read[j][secondcol+1:-1]))
                                    self.xshift[i]=int(float(read[j][firstcol+1:secondcol]))
                                    self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],self.xshift[i],axis=2)
                                    self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],self.yshift[i],axis=1)
            
                  f.close()

                  self.lbl.setText("Alignment using values from Text has been completed")
            except IOError:
                  print "choose file please"

      def saveAlignToText(self):
            try:
                  self.alignFileName = QtGui.QFileDialog.getSaveFileName()
                  f=open(self.alignFileName[0],"w")
                  for i in arange(self.projections):
                        onlyfilename=self.selectedFiles[i].rfind("/")
                        print self.selectedFiles[i]
                        f.writelines(self.selectedFiles[i][onlyfilename+1:]+", "+str(self.xshift[i])+", "+str(self.yshift[i])+ "\n")
                  f.close()
            except IOError:
                  print "choose file please"
                        

#==========================
      def xcorrelate(self,image1,image2,edgeguass_sigma=4):
            t0,t1=tomopy.algorithms.preprocess.align_projections.cross_correlate(image1,image2)
            return t0,t1

      def phasecorrelate(self,image1,image2):
            t0,t1=tomopy.algorithms.preprocess.align_projections.phase_correlate(image1,image2)
            return t0,t1
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

      def showSaveHotSpotPos(self):
            self.projViewControl.numb=len(self.channelname)
            for j in arange(self.projViewControl.numb):
                  self.projViewControl.combo.addItem(self.channelname[j])
            self.projViewControl.combo.currentIndexChanged.connect(self.saveHotSpotPos)
            self.projViewControl.btn.clicked.connect(self.alignHotSpotPos3)

      def saveHotSpotPos(self):
            self.projViewElement = self.projViewControl.combo.currentIndex()
            self.projView.data=self.data[self.projViewElement,:,:,:]
            self.projView.posMat=zeros([self.data.shape[1],2])
            self.projView.setImage(self.data[self.projViewElement,0,:,:])
            
      def alignHotSpotPos1(self):
            
            for i in arange(self.projections-1):
                  print "shifted"
                  
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], int(round(self.projView.posMat[0,0]-self.projView.posMat[i+1,0])), axis=2)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], int(round(self.projView.posMat[0,1]-self.projView.posMat[i+1,1])), axis=1)            

      def alignHotSpotPos2(self):
            for i in arange(self.projections-1):
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], int(round(self.projView.posMat[0,1]-self.projView.posMat[i+1,1])), axis=1)     
            x=self.theta
            self.fitfunc = lambda p,x: p[0]*sin(2*pi/360*(x-p[1]))+p[2]
            self.errfunc = lambda p,x,y: self.fitfunc(p,x)-y
            p0=[100,100,100]
            self.p1,success = optimize.leastsq(self.errfunc,p0[:],args=(x,self.projView.posMat[:,0]))
            self.hotSpotPosDiff=self.fitfunc(self.p1,x)-self.projView.posMat[:,0]

            for i in arange(self.projections):
                  self.xshift[i]=int(self.hotSpotPosDiff[i])
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], self.xshift[i], axis=2)
            self.lbl.setText("Alignment has been completed")
            
      def alignHotSpotPos3(self):
            #self.projView.data2=self.data[7,:,:,:]
            self.boxSize=40
            self.boxSize2=self.boxSize/2
            self.xPos=zeros(self.projections)
            self.yPos=zeros(self.projections)
            self.boxPos=zeros([self.projections,self.boxSize,self.boxSize])
            for i in arange(self.projections):
                  
                  self.xPos[i]=int(round(self.projView.posMat[i,0]))
                  self.yPos[i]=int(round(self.projView.posMat[i,1]))
                  #self.boxPos[i,:,:]=self.projView.data[i,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2]
                  self.boxPos[i,:,:]=self.projView.data[i,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2]
            print self.boxPos.shape
            pg.image(self.boxPos)
            self.alignHotSpotPos3_2()

      def alignHotSpotPos3_1(self):
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[0,:,:]=self.boxPos[0,:,:]
            self.boxShift=zeros([2,self.projections])
            for i in arange(self.projections-1):
                  #t0,t1=self.phasecorrelate(self.boxPos[i,:,:],self.boxPos[i+1,:,:])
                  t0,t1=self.xcorrelate(self.boxPos[i,:,:],self.boxPos[i+1,:,:])
                  print t0,t1
                  self.newBoxPos[i+1,:,:]=np.roll(np.roll(self.boxPos[i+1,:,:],t1,axis=1),t0,axis=0)
                  self.boxShift[0,i+1]=t0
                  self.boxShift[1,i+1]=t1
            self.testtest=pg.ImageView()
            self.testtest.setWindowTitle("rhyme")
            self.testtest.setImage(self.newBoxPos)
            self.testtest.show()

            for j in arange(self.projections-1):
                  self.data[:,j,:,:]=np.roll(np.roll(self.data[:,j,:,:],int(round(self.boxShift[1,j+1]-self.xPos[j+1]+self.xPos[0])),axis=2),int(round(self.boxShift[0,j+1]-self.yPos[j+1]+self.yPos[0])),axis=1)
            self.test2=pg.ImageView()
            self.test2.setWindowTitle("test2")
            self.test2.setImage(self.data[self.projViewElement,:,:,:])
            self.test2.show()

      def alignHotSpotPos3_2(self):
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[0,:,:]=self.boxPos[0,:,:]
            self.boxShift=zeros([2,self.projections])
            for i in arange(self.projections-1):
                  
                  img1=self.boxPos[i,:,:]
                  img2=self.boxPos[i+1,:,:]
                  
                  corr=scipy.signal.correlate2d(img1,img2,fillvalue=np.average(img1[0:2,0:2]))
                  t0, t1 = np.unravel_index(np.argmax(corr), corr.shape)
                  self.newBoxPos[i+1,:,:]=np.roll(self.boxPos[i+1,:,:],40-t0,axis=0)
                  self.newBoxPos[i+1,:,:]=np.roll(self.boxPos[i+1,:,:],40-t1,axis=1)
                  print t0, t1
            self.testtest=pg.ImageView()
            self.testtest.setWindowTitle("rhyme")
            self.testtest.setImage(self.newBoxPos)
            self.testtest.show()

      def alignHotSpotPos3_3(self):
            print "test"

#==========================
      def export_data(self):
            a=h5py.File("export_data.h5")
            print np.where(self.data==inf)
            
            a.create_dataset("data",data=self.data,compression="gzip")
            a.close()


#==========================

      def runTransReconstruct(self):

            self.recon.numb=len(self.channelname)
            for j in arange(self.recon.numb):
                  self.recon.combo.addItem(self.channelname[j])
            self.recon.show()
            self.recon.btn.setText("Reconstruction")
            self.recon.btn.clicked.connect(self.reconstruct)
            self.recon.save.clicked.connect(self.saveRecTiff)
            self.recon.reconvalue=1


      
      def runReconstruct(self):

            self.recon.numb=len(self.channelname)
            for j in arange(self.recon.numb):
                  self.recon.combo.addItem(self.channelname[j])
            self.recon.show()
            self.recon.btn.setText("Reconstruction")
            self.recon.btn.clicked.connect(self.reconstruct)
            self.recon.save.clicked.connect(self.saveRecTiff)
            self.recon.reconvalue=0
            
      def reconstruct(self):
            self.recon.lbl.setText("Reconstruction is currently running")
            self.d= tomopy.xtomo_dataset(log='debug')
            self.reconelement=self.recon.combo.currentIndex()
            self.d.data=self.data[self.reconelement,:,:,:]
            self.d.data[self.d.data==inf]=0.01
            self.d.data[np.isnan(self.d.data)]=0.01

            ###TEMP
            if self.recon.reconvalue==0:
                  self.d.data = (np.exp(-0.0001*self.d.data)).astype('float32')
            else:
                  print "transmission"
            #self.d.center=128
            ###TEMP
            self.d.center=self.p1[2]
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
            elif self.recon.method.currentIndex()==3:
                  self.d.dataset(self.d.data, theta=self.theta*np.pi/180)
                  #self.d.optimize_center()
                  self.d.pml()
                  
            print self.d.center

            self.reconView.setImage(self.d.data_recon)
            self.reconView.setWindowTitle("Slices of reconstructed model")
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
            print self.sit.combo.currentIndex()

            self.sit.combo.setCurrentIndex(i)
            for i in arange(30):
                  self.sit.method.removeItem(0)
            self.sit.secondColumnName = self.sit.firstColumn[self.sit.combo.currentIndex()][0]
            
            self.sit.secondColumn = self.sit.data[self.sit.secondColumnName].items()
            self.sit.secondColumnNum = len(self.sit.secondColumn)
            for j in arange(self.sit.secondColumnNum):
                  self.sit.method.addItem(self.sit.secondColumn[j][0])

            self.sit.method.setHidden(False)
            
            self.sit.combo.currentIndexChanged.connect(self.selectImageTag_image)
            self.sit.btn.setText("Set")

            self.sit.btn.clicked.connect(self.setImageTag)
            self.sit.show()

      def selectImageTag_image(self):
            for i in arange(30):
                  self.sit.method.removeItem(0)
            self.sit.secondColumnName = self.sit.firstColumn[self.sit.combo.currentIndex()][0]
            
            self.sit.secondColumn = self.sit.data[self.sit.secondColumnName].items()
            self.sit.secondColumnNum = len(self.sit.secondColumn)
            for j in arange(self.sit.secondColumnNum):
                  self.sit.method.addItem(self.sit.secondColumn[j][0])

            self.sit.method.setHidden(False)

      def setImageTag(self):
            self.ImageTag = self.sit.combo.currentText()
            self.lbl.setText("Image Tag has been set to \""+self.ImageTag+"\"")
            print "Image Tag has been set to \"", self.ImageTag, "\""
#==============================


            
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
            degree_sign= u'\N{DEGREE SIGN}'
            for i in arange(len(self.fileNames)):
                  self.fileNames[i]=str(self.fileNames[i])
                  f = h5py.File(os.path.abspath(self.fileNames[i]),"r")
                  thetatemp = f["MAPS"]["extra_pvs_as_csv"][self.thetaPos]
                  thetapos = string.rfind(thetatemp, ",")
                  theta = str(round(float(thetatemp[thetapos+1:])))
                  onlyfilename=self.fileNames[i].rfind("_")
                  self.filecheck.button[i].setText(self.fileNames[i][onlyfilename+1:]+" ("+theta+degree_sign+")")
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
                  self.fileNames=str(fileNametemp.join("\n")).split("\n")
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
                  folderName=str(folderName)
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
            self.channelname=f[self.ImageTag]["channel_names"]

            self.channels,self.y,self.x=f[self.ImageTag]["data"].shape
            self.projections=len(self.selectedFiles)
            self.theta= zeros(self.projections)
            self.data=zeros([self.channels,self.projections,self.y,self.x])

            self.xshift=zeros(self.projections,int)
            self.yshift=zeros(self.projections,int)

            
            

            for i in arange(self.projections):
                  file_name = os.path.abspath(self.selectedFiles[i])
                  f = h5py.File(file_name,"r")
                  thetatemp=f["MAPS"]["extra_pvs_as_csv"][self.thetaPos]

                  self.theta[i] = float(thetatemp[thetatemp.rfind(",")+3:])
                  
                  for j in arange(self.channels):
                        self.data[j,i,:,:-2]=f[self.ImageTag]["data"][j,:,:-2]
                  print i+1, "projection(s) has/have been converted"
            print "worked"

            self.data[isnan(self.data)]=0.0001
            self.data[self.data==inf]=0.0001
            self.p1=[100,100,self.data.shape[3]/2]

            self.alignmentMenu.setEnabled(True)
            self.tab_widget.setEnabled(True)
            
            self.oldData=zeros(self.data.shape)
            self.oldData[...]=self.data[...]

            self.projView.hotSpotNumb=0
            self.showProjections()
            self.showSinogram()
            self.sinogram()
            self.viewProjections()
            self.runReconstruct()
            self.showSaveHotSpotPos()
            
#            yy=self.data
#            yy[isinf(yy)]=1
#            pg.image(yy[37,:,:,:])


#####!!!! just temp need to be fixed
      def saveImage(self):
            self.saveImageDir=QtGui.QFileDialog.getExistingDirectory()
            print self.saveImageDir
            for j in arange(self.data.shape[0]):
                  path=self.saveImageDir+"/"+self.channelname[j]
                  try:
                        os.makedirs(path)
                  except OSError:
                        path=path
                  for i in arange(self.data.shape[1]):
                        temp_img=self.data[j,i,:,:]
                        temp=Image.fromarray(temp_img.astype(np.float32))
                        index=string.rfind(self.fileNames[i],"/")
                        temp.save(path+self.fileNames[i][index:-3]+"_"+self.channelname[j]+".tif")

      def showProjections(self):
            self.projection.numb=len(self.channelname)
            for j in arange(self.projection.numb):
                  self.projection.combo.addItem(self.channelname[j])

            self.projection.btn.setVisible(False)
            self.projection.btn2.setVisible(False)
            self.projection.sld.setVisible(False)
            self.projection.lcd.setVisible(False)
            self.projection.combo.currentIndexChanged.connect(self.viewProjections)

      def viewProjections(self):
            self.projelement=self.projection.combo.currentIndex()
            projdata=self.data[self.projelement,:,:,:]
            self.projectionView.setImage(projdata)




#####!!!!

##      def sinogram(self):
##            self.thickness=10
##            self.sinoelement=self.sino.combo.currentIndex()
##            sinodata=self.data[self.sinoelement,:,:,:]
##            sinogram=zeros([sinodata.shape[0]*self.thickness,sinodata.shape[2]],dtype=float32)
##            for i in arange(self.projections):
##                  for j in arange(self.thickness):
##                        sinogram[i*self.thickness+j,:]=sinodata[i,sinodata.shape[1]/2,:]
##            sinogram[isinf(sinogram)]=1
##
##            iView=IView()
##            iView.show()
##            iView.setImage(sinogram)
##      def sinogram(self):
##            self.thickness=10
##            self.sinoelement=self.sino.combo.currentIndex()
##            sinodata=self.data[self.sinoelement,:,:,:]
##            self.sinoz=self.data.shape[2]-9
##            sinogram=zeros([self.sinoz,sinodata.shape[0]*self.thickness,sinodata.shape[2]],dtype=float32)
##            for j in arange(self.sinoz):
##                  sinoZ = j + 5
##                  for i in arange(self.projections):
##                        sinogram[j,i*self.thickness:(i+1)*self.thickness,:]=sinodata[i,sinoZ-self.thickness/2:sinoZ+self.thickness/2,:]
##
##            sinogram[isinf(sinogram)]=0.001
##
##            self.iView=IView()
##            self.iView.show()
##            self.iView.setImage(sinogram)
##            self.iView.getShape()

      def showSinogram(self):
            #self.sino = QSelect2()
            #self.sino.setWindowTitle("Sinogram Window")


            self.sino.numb=len(self.channelname)
            for j in arange(self.sino.numb):
                  self.sino.combo.addItem(self.channelname[j])
            self.sino.show()
            self.sino.btn.clicked.connect(self.runCenterOfMass2)
            self.sino.btn.setText("center of mass")
            self.sino.btn2.clicked.connect(self.sinoShift)
            self.sino.sld.setRange(5,self.y-5)
            self.sino.lcd.display(5)
            self.sino.sld.valueChanged.connect(self.sino.lcd.display)
            self.sino.sld.valueChanged.connect(self.sinogram)
            
            self.sinoView.show()

      def sinogram(self):
            self.thickness=10
            self.sinoelement=self.sino.combo.currentIndex()
            sinodata=self.data[self.sinoelement,:,:,:]
            self.sinogramData=zeros([sinodata.shape[0]*self.thickness,sinodata.shape[2]],dtype=float32)

            for i in arange(self.projections):
                  self.sinogramData[i*self.thickness:(i+1)*self.thickness,:]=sinodata[i,self.sino.sld.value()-self.thickness/2:self.sino.sld.value()+self.thickness/2,:]

            self.sinogramData[isinf(self.sinogramData)]=0.001
            self.sinoView.setImage(self.sinogramData)
            self.sinoView.getShape()
            self.sinoView.setWindowTitle("Sinogram "+self.sino.combo.itemText(self.sinoelement)+" "+str(self.sino.sld.value()))



      def sinoShift(self):
            for i in arange(self.projections):
                  
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], self.iView.regShift[i], axis=2)
            
            

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
                  self.selectFiles()                                   
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
            self.btn.setText("Save Data in Memory")
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
            
            self.vb.addWidget(self.lbl,11)
            self.vb.addWidget(self.lbl2,12)

            self.vb2.addWidget(self.btn2,13)
            self.vb2.addWidget(self.btn,14)

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
 
            self.sld=QtGui.QSlider(QtCore.Qt.Horizontal, self)
            self.lcd=QtGui.QLCDNumber(self)
            self.combo = QtGui.QComboBox(self)
            self.btn = QtGui.QPushButton('Click2')
            self.btn.setText("Sinogram")
            self.btn2 = QtGui.QPushButton("shift data")
            self.lbl=QtGui.QLabel()
            self.lbl.setText("")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.btn)
            vb.addWidget(self.btn2)
            vb.addWidget(self.lcd)
            vb.addWidget(self.sld)
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
            self.methodname=["mlem", "gridrec", "art","pml"]
            for k in arange(4):
                  self.method.addItem(self.methodname[k])
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.method)
            vb.addWidget(self.btn)
            vb.addWidget(self.save)
            vb.addWidget(self.lbl)
            self.setLayout(vb)

class QSelect4(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect4, self).__init__()
        
            self.initUI()
        
      def initUI(self):
 
            self.sld=QtGui.QSlider(QtCore.Qt.Horizontal, self)
            self.lcd=QtGui.QLCDNumber(self)
            hb1=QtGui.QVBoxLayout()
            hb2=QtGui.QVBoxLayout()
            self.lbl11=QtGui.QLabel()
            self.lbl11.setText("Current Position")
            self.lbl12=QtGui.QLabel()
            self.lbl21=QtGui.QLabel()
            self.lbl21.setText("Total Position")
            self.lbl22=QtGui.QLabel()
            hb1.addWidget(self.lbl11)
            hb1.addWidget(self.lbl12)
            hb2.addWidget(self.lbl21)
            hb2.addWidget(self.lbl22)
            self.combo = QtGui.QComboBox(self)
            self.btn=QtGui.QPushButton("Draw boxes")
##            self.btn = QtGui.QPushButton('Click2')
##            self.btn.setText("Sinogram")
##            self.btn2 = QtGui.QPushButton("shift data")
##            self.lbl=QtGui.QLabel()
##            self.lbl.setText("")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
##            vb.addWidget(self.btn)
##            vb.addWidget(self.btn2)
            vb.addWidget(self.lcd)
            vb.addWidget(self.sld)
            vb.addLayout(hb1)
            vb.addLayout(hb2)
            vb.addWidget(self.btn)
##            vb.addWidget(self.lbl)
            self.setLayout(vb)

class IView(pg.ImageView):
    
      def __init__(self):
            super(IView, self).__init__()
        
            self.initUI()

            self.hotSpotNumb=0
        
      def initUI(self):
            self.show()

      def keyPressEvent(self, ev):
            if ev.key() == QtCore.Qt.Key_M:
                  self.getMousePos()
                  self.shiftnumb = int(self.imageItem.difY)
                  self.shift()
                  self.setImage(self.copy)
                  self.regShift[self.numb2]+=self.shiftnumb
                  
            if ev.key() == QtCore.Qt.Key_Down:

                  self.getMousePos()
                  self.shiftnumb=1
                  self.shift()
                  self.setImage(self.copy)
                  self.regShift[self.numb2]+=self.shiftnumb

            if ev.key() == QtCore.Qt.Key_Up:

                  self.getMousePos()
                  self.shiftnumb=-1
                  self.shift()
                  self.setImage(self.copy)
                  self.regShift[self.numb2]+=self.shiftnumb

            if ev.key() == QtCore.Qt.Key_N:
                  print "n"
                  print self.data.shape[0], self.hotSpotNumb
                  print self.imageItem.iniX, self.imageItem.iniY
                  self.posMat[self.hotSpotNumb,0]=self.imageItem.iniY
                  self.posMat[self.hotSpotNumb,1]=self.imageItem.iniX
                  if self.hotSpotNumb<self.data.shape[0]-1:
                        self.hotSpotNumb+=1
                        self.setImage(self.data[self.hotSpotNumb,:,:])

            if ev.key() ==QtCore.Qt.Key_A:
                  
                  self.posMat[self.hotSpotNumb,0]=self.imageItem.iniY
                  self.posMat[self.hotSpotNumb,1]=self.imageItem.iniX
                  if self.hotSpotNumb<self.data.shape[0]-1:
                        self.hotSpotNumb+=1
                        self.setImage(self.data[self.hotSpotNumb,:,:])
                  
      def getMousePos(self):
            numb=self.imageItem.iniX
            self.numb2=int(numb/10)


      def shift(self):
            self.copy=self.getProcessedImage()
            self.copy[self.numb2*10:self.numb2*10+10,:]=np.roll(self.copy[self.numb2*10:self.numb2*10+10,:],self.shiftnumb,axis=1)
      def getShape(self):
            self.regShift=zeros(self.getProcessedImage().shape[0],dtype=int)
            
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
