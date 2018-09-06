#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 12:54:48 2018

@author: lucynewnham
"""


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import interpolation as ipt
import scipy



def usage():
    print('barlengthvshi expects 3 params:')
    print('\tpython barlengthvshi_fixedy.py [input file (mom0 map)] [rotation - dont include -ve signs.] [pixel scale]\n\n')

    
def barlengthvshi_fixedy(fname, rotation, pix):
    hifits= fits.open(fname)
    hi= hifits[0].data

    #(data, value to rotate by in degrees)
    hirotate = ipt.rotate(hi,rotation)
    shape = hirotate.shape
    
    row=int(shape[0]/2)
    row1=row+1
    row2=row-1
    
    #this is adding the rows that we want to look at together
    m=hirotate[row]+hirotate[row1]+hirotate[row2]
    
    fig=plt.figure()
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.5)    
    #fig,ax=plt.subplots(1,2,sharey=False,figsize=(10,2))
    #print(ax)
    #ax1, ax2 = ax[0]
    #fig.subplots_adjust(wspace=0,hspace=0)
    
    #plt.figure(0)
    
    #Use this to see if the image has been rotated
    im = ax1.imshow(hirotate, origin="lower")
    ax1.plot(np.zeros(hirotate.shape[1])+row, color="steelblue")
    ax1.plot(np.zeros(hirotate.shape[1])+row1, color="steelblue")
    ax1.plot(np.zeros(hirotate.shape[1])+row2, color="steelblue")
    
    diff=row/3
    
    lim1=row-diff
    lim2=row+diff
    ax1.set_xlim(lim1,lim2)
    ax1.set_ylim(lim1,lim2)
    
    ticks = ax1.get_xticks()
    x=(ticks-row)*pix
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(x.astype(int))

    ticksy = ax1.get_yticks()
    y=(ticksy-row)*pix
    ax1.set_yticks(ticksy)
    ax1.set_yticklabels(y.astype(int))
    ax1.set_xlabel('Arcsecs from centre')
    ax1.set_ylabel('Arcsecs from centre')
    
    cbaxes = fig.add_axes([0.5, 0.1, 0.03, 0.8]) 
    fig.colorbar(im,cbaxes)
    
    #plt.show()
    #number is the channel on the y axis you want to look at.
    
    #plt.figure(1)
    
    #plots the one pixel and the surrounding rows according to m
    ax2.plot(m)
    ax2.set_xlim(lim1,lim2)
    ticks2 = ax2.get_xticks()
    x2=(ticks2-row)*pix
    ax2.set_xticks(ticks2)
    ax2.set_xticklabels(x2.astype(int))
    ax2.set_xlabel('Arcsecs from centre')

    

    #plots only the one pixel row
    #plt.plot(hirotate[y])
    
    plt.show()    
    
if __name__ == "__main__":
#    if len(sys.argv) != 5 or not os.path.isfile(sys.argv[1]) or not os.path.isfile(sys.argv[3]):
#        usage()
#    else:
#        print(spectraoverlay(*sys.argv[1:]))
#we use the above if all the inputs are strings and the below if they are not all strings
#this allows us to covert them all to strings to make it work.
    if len(sys.argv)<4:
        usage()
    else:
        print(barlengthvshi_fixedy(*[i(j) for i,j in zip([str, float, float], sys.argv[1:])]))




























