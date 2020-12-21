# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 10:45:13 2020

@author: Rhys Shaw

This program calculates the Morphology metric as discussed in my Final year project. It also estimates it error using a bootstrapping technique. Included in this folder is a
.fits file what will need to be downloaded to test and use this code.
"""

from astropy.io import fits
from astropy.table import Table
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import time

class morph_statistic:
    
    def extract_fits(fitsfile,N):
        hdu_list = fits.open(fitsfile) #make sure you look at directry
        
        hdu_list.info()                #this gives you ta types in the file you will need the primary 

        evt_data = Table(hdu_list[1].data)
        x = evt_data['x']
        y = evt_data['y']
        NBINS = (N,N)

        A = np.histogram2d(x,y,NBINS)[0]
        #print(A)
        plt.imshow(A,norm=LogNorm())
        return A

    def cluster_symetry(A,N):
        n0 = 0
        Atot = 0
        for i in range(0,len(A)-1):
            for j in range(0,len(A)-1):
                if A[i,j] != 0:
                    n0 = 1 + n0
                    Atot =Atot + A[i,j]
    
        cpx = cpy = N/2
        cp = np.array([[cpx],[cpy]])
        res = 0
        for i in range(0,N-1):
            for j in range(0,N-1):
                if i == cpx:
                    if j ==cpy:
                        res = res + 0
                else:
                    point = np.array([[i],[j]])
                    x = np.subtract(point,cp)
                    y = A[i,j]
                    res = np.add(res,x*y)  
    
        meanA = Atot/n0
        vec = res/(meanA*n0)
    #print('Mean number of counts: ',meanA)     
        c = np.linalg.norm(res)
        sym = c/(meanA*n0)
    
    #print('Normalised vector:     ',c)
    #print('Symetry value:         ',c/(meanA*n0))

        return sym, c, vec


    def sym_error(A,N,Ni):
        
        disym = []
        for n in range(0,Ni):
        
            #print(n)
            An = np.zeros((len(A),len(A)))
            for i in range(0,len(A)-1):
                for j in range(0,len(A)-1):
                    if A[i,j] != 0:
                        #print(A[i,j])
                        An[i,j] = np.random.poisson(A[i,j])
            syme = morph_statistic.cluster_symetry(An,N)[0]
            disym.append(syme)
        #return di
        
        Mean = np.mean(disym)
        sd = np.std(disym)
        #fig, ax = plt.subplots()
        #n, bins, patches = ax.hist(disym, 100, density=1)
        #ax.set_xlabel('Morphology Metric')
        #ax.set_ylabel('Frequecy')
        #plt.show()
        
        return Mean, sd
        

fitsfile = np.array([['13517_remorph.fits',62]])
N = int(fitsfile[0,1])*2
A = morph_statistic.extract_fits(fitsfile[0,0], N)
morph_stat, norm ,vec = morph_statistic.cluster_symetry(A,N)
Mean, sd = morph_statistic.sym_error(A,N, 100) #iteration number is made smaller so that computational time is slower
print('Morph stat: ',morph_stat)
print('Mean :',Mean)
print('Sd :',sd)

