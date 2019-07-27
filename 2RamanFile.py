# NOTE: This first cell has code to take 2 files, pull the relevant peaks from each, and find a 2D/G Raman ratio.
# File names ought to be the same, except for the ending charatcters which end in 2D and G.
# Future changes to the code can modify it this as needed.
# The other cells further below are used for files from "extended scans" with both peaks in 1 file.

###This is the code for sorted files (C1 folder, C2 folder...)
import numpy as np
#import simplejson
import pandas as pd
import matplotlib.pyplot as plt
import random
import glob
import itertools
import os

import plotly.plotly as py
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF

import scipy
from scipy.signal import chirp, find_peaks, peak_widths
import peakutils
from peakutils.plot import plot as pplot
from matplotlib import pyplot



# Note, there are 2 things we can do. 1: Count all txt files with glob and stores that value into an txtcount variable
# Use txtcount to run through a for loop and run the below code on each data file
# OR 2: Parse each txt file in folder, run slice/plot code on each file as it is found.
# This code uses method one. I'd like to change it to method 2 if possible.

###Change myPath
myPath = '/Users/DerekChang/Documents/Ragan Lab/Test/063019-BijelBoxC2-C3 copy 2/C2'
SamplePath = '/Users/DerekChang/Documents/Ragan Lab/Test/063019-BijelBoxC2-C3 copy 2/C2'
CombinedPath = '/Users/DerekChang/Documents/Ragan Lab/Test/063019-BijelBoxC2-C3 copy 2/C2/Combined'
os.chdir(myPath)

FileNames = glob.glob1(myPath, '*.txt') #list(?)
DataNames = glob.glob1(CombinedPath, '*.txt')

numberofplots = len(FileNames)


###Declared variables
list = {}
list2D = {}
listG = {}
ratio = {}
a = {}
b = {}
VisN = 0
Number = 0
Name2D = {}
NameG = 0
avgRatio = 0
lengthR = 0
sumR = 0
wavenumber = {}
wavenumber1 = {}

NewFile = []
intensity = {}
intensity1 = {}

filenames = {}
fileNamesTwoD = {}
fileNamesG = {}
gPeakArray = {}
twoDPeakArray = {}
# initialize G and 2D counter iterator for for loops
gCounter = 0
twoDCounter = 0


fig, axs = plt.subplots(numberofplots, 1, figsize=(10, 60), constrained_layout=False)
fig.tight_layout()
plt.rc('axes', labelsize='20')
plt.rc('axes', titlesize='20')

for counter, datafile in enumerate(FileNames):  # enumerate splits a list up into a counter and its value (in this case, filenames.txt) into a 2 column list
    path = myPath + '/' + datafile
    # print(counter, datafile) #Should rename to filename later.
    #   print('counter is %i' % counter)
    # Save each filename in a list
    filenames[counter] = datafile

    data_file = np.loadtxt(path, delimiter='\t')  # put skiprows=1 into loadtxt for new sets of data with headers still in.
    wavenumber[counter] = data_file[:, 0]
    intensity[counter] = data_file[:, 1]
    # Strip txt from datafile and change to uppercase.
    strippedDatafile = datafile.rstrip('.txt')
    strippedDatafile = strippedDatafile.replace('-', '')
    strippedDatafile = strippedDatafile.upper()

    x = wavenumber[counter]

    # Subtract baseline
    y2 = intensity[counter] + np.polyval([0.002, -0.08, 5], wavenumber[counter])
    base = peakutils.baseline(y2, 2)
    y3 = y2 - base

    # Find peaks of baseline-subtracted curve
    indexes3 = peakutils.indexes(y3, thres=0.12, min_dist=90)
    print(strippedDatafile)
    print(len(indexes3))

    # Code added 032019 to select intensities from known X-range in case of multiple peaks
    for iterator, peakpos in enumerate(x[indexes3]):
        print(peakpos)
        if peakpos > 1540 and peakpos < 1600:
            fileNamesG[gCounter] = np.array(strippedDatafile[:-1])
            #   print('filename G is %s' % fileNamesG[counter])
            #listG = np.array([gCounter] = np.array(data_file)
            listG[gCounter] = np.array(data_file)
            gpeak = y3[indexes3[iterator]]
            gPeakArray[gCounter] = y3[indexes3[iterator]]
            a[gCounter] = np.array(gPeakArray[gCounter])

            gCounter = gCounter + 1
            #   print('gCounter is %i' % gCounter)
            # x[indexes3[iterator]] #This is the x position of the G peak.
            # y3[indexes3[iterator]] #This is the intensity.
            print('G indexes are %f' % y3[indexes3[iterator]])
            # a = x[indexes3[iterator]]
        if peakpos > 2660 and peakpos < 2730:
            fileNamesTwoD[twoDCounter] = np.array(strippedDatafile[:-2])
            # print('saved filename 2D is %s' % fileNamesTwoD[counter])
            list2D[twoDCounter] = np.array(data_file)
            twodpeak = y3[indexes3[iterator]]
            twoDPeakArray[twoDCounter] = y3[indexes3[iterator]]
            b[twoDCounter] = np.array(twoDPeakArray[twoDCounter])
            twoDCounter = twoDCounter + 1
            # b = x[indexes3[iterator]]
            print('2D indexes are %f' % y3[indexes3[iterator]])  # this is 2D intensity
        # print(b/a)

        #print(list2D[i])

        #print('this is new file', NewFile)
        #NewFile[i] = listG[i].append(list2D[i])
    #f = open("newfiletesting", "w")

    #    f.write(listG[0])
    #f.close()

    #with open('newfiletesting.txt', 'w') as f:
    #    for listG in range(len(listG)):
    #        f.write("%s\n" % listG)

    #   print('twoDCounter is %i' % twoDCounter)
    # This needs to be executed on the last run of this for loop
    # intensity_ratio_2dg = twoDPeakArray[twoDCounter]/gPeakArray[gCounter]
    #   print("gpeak is %f"% gpeak)a

    # Graphing
    # axs[counter].plot(wavenumber[counter], y3)
    # axs[counter].set_xlabel(r'$Wavenumber\ (cm^{-1})$')
    # axs[counter].set_ylabel('$Intensity\ (a.u.)$')
    # if len(FileNames) - 1 == counter:
    #     plt.show()
    #     plt.savefig('Graphene-on-quartz.png')

### Matching
for i in range(len(fileNamesTwoD)):
    for j in range(len(fileNamesG)):
        if np.array_equal(fileNamesTwoD[i] ,fileNamesG[j]):
            ratio[i] = np.array(b[i] / a[j])
            def createFolder(directory):         ##creates folder if there isn't one exist already
                try:
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                except OSError:
                    print('Error: Creating directory. ' + directory)
            createFolder('./Combined/')          ##folder name
## try to append the text files here listG[j] with list2D[i]
            os.chdir(CombinedPath)
            with open("2DG_Ratio_Combined.txt", "wb") as CR:
            # np.savetxt(listG[j].append(list2D[i]))
                np.savetxt(CR, list2D[i], delimiter='\t')
                np.savetxt(CR, listG[i], delimiter='\t')
            # print(listG[j])
            # np.savetxt('2DG_Ratio_Combined.txt',listG[j], delimiter = '\t')
                VisN = len(fileNamesTwoD)-1
            # print('this is VisN', VisN)
                NameG = fileNamesTwoD[j]
                os.rename('2DG_Ratio_Combined.txt', '%s' %NameG)
                os.chdir(SamplePath)
                if fileNamesTwoD[i] != fileNamesG[j]:
                    j = j + 1
                if fileNamesTwoD[i] == len(fileNamesG)-1:
                    break

###Graphing
for counter2, dataNames in enumerate(DataNames):
    os.chdir(CombinedPath)
    CombinedPath = CombinedPath + '/' + dataNames
    dataNames = np.loadtxt(CombinedPath,delimiter='\t')  
    wavenumber1[counter2] = dataNames[:, 0]
    intensity1[counter2] = dataNames[:, 1]
    print('this is date name',dataNames)
    axs[counter2].plot(wavenumber1[counter2], intensity1[counter2])
    axs[counter2].set_xlabel(r'$Wavenumber\ (cm^{-1})$')
    axs[counter2].set_ylabel('$Intensity\ (a.u.)$')
    if len(DataNames) - 1 == counter2:
        plt.show()
        os.chdir(SamplePath)
        plt.savefig('Graphene-on-quartz.png')

# f = open('Test/output.txt', 'w')
# for i in range(len(listG)):
#    print('this is list', list2D)
# for i in range(len(listG)):
#    f.write("%s ", list)
# f.close()

#with open('Test/output.txt', 'w') as f:
#    for item,name in listG:
#        f.write("%s \n" item)

#print('this is a = ', a)
#print('\tthis is g = ', a[5])
#print('this is b = ', b)
#print('this is ratio = ', b / a)

# print('this is listg = ', listG)
# print('file name',fileNamesTwoD)


for w in range(len(ratio)):    ##this calculates the average ratio
     sumR = sumR + ratio[w]
     if w > len(ratio):
         break
avgRatio = sumR/len(ratio)

print('this is average', avgRatio)
# print('this is names ', NameG)
