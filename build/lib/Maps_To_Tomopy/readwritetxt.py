import string
from scipy import *

def readtxt(textfile=None):
    if textfile==None:
        textfile="configuration.txt"
    f=open(textfile)
    l=f.readlines()
    for i in arange(len(l)):
        ## file name
        if string.find(l[i],"file")!=-1:
            filen=l[i+1][:-1]
        ## folder name
        if string.find(l[i],"folder")!=-1:
            if l[i+1][-2]!='/':
                foldern=l[i+1][:-1]+'/'
            else:
                foldern=l[i+1][:-1]
        ## projection starts
        if string.find(l[i],"start")!=-1:
            projections_start=l[i+1][:-1]
        ## projections end
        if string.find(l[i],"end")!=-1:
            projections_end=l[i+1][:-1]
        ## excluding numbers
        if string.find(l[i],"exclude")!=-1:
            count=string.count(l[i+1][:-1],',')+1
            exclude_numbers=zeros(count,int)
            temp=l[i+1][:-1]
            for j in arange(count-1):
                exclude_numbers[j]=int(temp[:string.find(temp,',')])
                temp=temp[string.find(temp,',')+1:]
            try:
                exclude_numbers[-1]=int(temp)
            except ValueError:
                exclude_numbers=array([])
        
        ## element
        if string.find(l[i],"element")!=-1:
            if string.find(l[i],"number")!=-1:
                count=string.count(l[i+1][:-1],',')+1
                element=zeros(count,int)
                temp=l[i+1][:-1]
                for i in arange(count-1):
                    element[i]=int(temp[:string.find(temp,',')])
                    temp=temp[string.find(temp,',')+1:]
                element[-1]=int(temp)
        ## reconstruct method
        if string.find(l[i],"reconstruct")!=-1:
            rtype=l[i+1][:-1]
        ## save sinogram
        if string.find(l[i],"sinogram")!=-1:
            if string.find(l[i+1],"yes")!=-1:
                sinogramsave=True
            else:
                sinogramsave=False
        ## cross correlation
        if string.find(l[i],"correlation")!=-1:
            count=string.count(l[i+2][:-1],',')+1
            cor_number=zeros(2,int)
            temp=l[i+2][:-1]
            for j in arange(count-1):
                cor_number[j]=int(temp[:string.find(temp,',')])
                temp=temp[string.find(temp,',')+1:]
            try:
                cor_number[-1]=int(temp)
            except ValueError:
                cor_number=array([])
    filename=foldern+filen
    return filename,projections_start,projections_end,exclude_numbers,element,rtype,sinogramsave,cor_number

def writetxt(elementlist):
    f=open("./elementlist.txt","w")
    y=elementlist[...]
    for i in arange(y.shape[0]):
        f.write(str(i)+"       "+y[i]+"\n")
    f.close()
