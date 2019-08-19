# NOTE: This first cell has code to take 2 files, pull the relevant peaks from each, and find a 2D/G Raman ratio.
# File names ought to be the same, except for the ending charatcters which end in 2D and G.
# Future changes to the code can modify it this as needed.
# The other cells further below are used for files from "extended scans" with both peaks in 1 file.

###This is the code for sorted files (C1 folder, C2 folder...)
import numpy as np
from array import *
#import simplejson
import pandas as pd
import array
import matplotlib.pyplot as plt
import random
import glob
import itertools
import os
import math
import plotly.plotly as py
import plotly.graph_objs as go
import statistics as stat
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
myPath = '//Users/DerekChang/Documents/Ragan Lab/Characterization/08142019-Bx7 C1-C5/C1'
CombinedPath = '/Users/DerekChang/Documents/Ragan Lab/Characterization/08142019-Bx7 C1-C5/C1/Combined'
os.chdir(myPath)

# SamplePath = '/Users/DerekChang/Documents/Ragan Lab/Test/072919-NiFilmsBox6_D2-D6/D3'

FileNames = glob.glob1(myPath, '*.txt') #list(?)

numberofplots = len(FileNames)


#[] array
#{} dictionary
###Declared variables
list = []
list2D = []
listG = []
ratio = []
a = []
b = []
timmer = 0
VisN = 0
Number = 0
Name2D = []
NameG = 0
avgRatio = 0
stdv = 0
lengthR = 0
sumR = 0
sum1 = 0
extra = 0
N = 0
wavenumber = []
wavenumber1 = []
ratio_list = []
NewFile = []
intensity = []
intensity1 = []
filenames = []
filenames2 = []
fileNamesTwoD = []
fileNamesG = []
gPeakArray = []
twoDPeakArray = []
gCounter = 0
twoDCounter = 0


fig, axs = plt.subplots(numberofplots, 1, figsize=(30, 25), constrained_layout=False)  ###nice figures
# fig, axs = plt.subplots(numberofplots, 1, figsize=(5, 10), constrained_layout=False)  ###Computation
fig.tight_layout()
plt.rc('axes', labelsize='10')
plt.rc('axes', titlesize='10')
print('')

for counter, datafile in enumerate(FileNames):  # enumerate splits a list up into a counter and its value (in this case, filenames.txt) into a 2 column list
    path = myPath + '/' + datafile
    # print(counter, datafile) #Should rename to filename later.
    #   print('counter is %i' % counter)
    # Save each filename in a list
    filenames.append(datafile)

    data_file = np.loadtxt(path, delimiter='\t', dtype="float")  # put skiprows=1 into loadtxt for new sets of data with headers still in.
    wavenumber.append(data_file[:, 0])
    intensity.append(data_file[:, 1])
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
    # print(strippedDatafile)
    # print(len(indexes3))

    # Code added 032019 to select intensities from known X-range in case of multiple peaks
    for iterator, peakpos in enumerate(x[indexes3]):
        # print(peakpos)
        if peakpos > 1540 and peakpos < 1600:
            fileNamesG.append(np.array(strippedDatafile[:-1]))
            #   print('filename G is %s' % fileNamesG[counter])
            #listG = np.array([gCounter] = np.array(data_file)
            listG.append(np.array(data_file))
            gpeak = y3[indexes3[iterator]]
            gPeakArray.append(y3[indexes3[iterator]])
            a.append(np.array(gPeakArray[gCounter]))

            gCounter = gCounter + 1
            #   print('gCounter is %i' % gCounter)
            # x[indexes3[iterator]] #This is the x position of the G peak.
            # y3[indexes3[iterator]] #This is the intensity.
            # print('G indexes are %f' % y3[indexes3[iterator]])
            # a = x[indexes3[iterator]]
        if peakpos > 2660 and peakpos < 2730:
            fileNamesTwoD.append(np.array(strippedDatafile[:-2]))
            # print('saved filename 2D is %s' % fileNamesTwoD[counter])
            list2D.append(np.array(data_file))
            twodpeak = y3[indexes3[iterator]]
            twoDPeakArray.append(y3[indexes3[iterator]])
            b.append(np.array(twoDPeakArray[twoDCounter]))
            twoDCounter = twoDCounter + 1
            # b = x[indexes3[iterator]]
            # print('2D indexes are %f' % y3[indexes3[iterator]])  # this is 2D intensity
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

##delete extra plots
for w in range(len(fileNamesTwoD)):
    extra = len(fileNamesTwoD)+w
    fig.delaxes(axs[extra])
### Matching
for i in range(len(fileNamesTwoD)):
    timmer = len(fileNamesTwoD) - i
    print('Ready in %s' % timmer + "...")
    for j in range(len(fileNamesG)):
        if np.array_equal(fileNamesTwoD[i] ,fileNamesG[j]):
            # ratio[i] = np.array(b[i] / a[j], dtype="float")
            ratio.append(np.array(b[i] / a[j], dtype="float"))
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
                np.savetxt(CR, listG[j], delimiter='\t')
            # print(listG[j])
            # np.savetxt('2DG_Ratio_Combined.txt',listG[j], delimiter = '\t')
            #     VisN = len(fileNamesTwoD)-1
            # print('this is VisN', VisN)
                NameG = fileNamesTwoD[i]
                os.rename('2DG_Ratio_Combined.txt', '%s' %NameG+".txt")
                # print('Name:', NameG)
                DataNames = glob.glob1(CombinedPath, '*.txt')
                for counter2, dataNames in enumerate(DataNames):
                    os.chdir(CombinedPath)
                    Combined_Path = CombinedPath + '/' + dataNames
                    filenames2.append(dataNames)
                    data_Names = np.loadtxt(Combined_Path, delimiter='\t', dtype="float")
                    wavenumber1.append(data_Names[:, 0])
                    intensity1.append(data_Names[:, 1])
                    axs[counter2].plot(wavenumber1[counter2], intensity1[counter2])
                    axs[counter2].set_xlabel(r'$Wavenumber\ (cm^{-1})$')
                    axs[counter2].set_ylabel('$Intensity\ (a.u.)$')
                    if len(DataNames) - 1 == counter2:
                        os.chdir(myPath)
                        plt.savefig('Graphene-on-quartz.png')
                os.chdir(myPath)
            if fileNamesTwoD[i] != fileNamesG[j]:
                j = j + 1
                break
        if fileNamesTwoD[i] == len(fileNamesG)-1:
            break

print('........done!')
# plt.show()


###Graphing


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
# print(ratio)
# print(len(ratio))

# for w in range(len(ratio)):    ##this calculates the average ratio
#      sumR = sumR + ratio[w]
#      if w > len(ratio):
#          break\

from array import *


# print('this is array', ratio)
# print('this is array list', ratio_list)
# print('this is list ', ratio_list)
#this is wrong, counting the index but not data
# print('this is ratio ', ratio)

# avgRatio = sumR/len(ratio)
ratio_list = np.array(ratio).tolist()
avgRatio = stat.mean(ratio_list)
stdv = stat.stdev(ratio_list)

print('')
print('Average = ', avgRatio)
print('Standard Deviation = ', stdv)
