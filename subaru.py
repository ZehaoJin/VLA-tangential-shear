# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 00:51:01 2018

@author: zehaojin
"""

import numpy as np
import matplotlib.pyplot as plt


cluster='0416'

subarucatalog=np.loadtxt('hlsp_clash_subaru_suprimecam_macs%s_photoz-cat.txt' %cluster,usecols=(0,12,11,13,14,17,19,10))  
                                          #(0_id,12_zbest,11_uncertainty in z,13_zbmin,14_zbmax,17_zml,19_Poorness of BPZ fit,10_z_Subaru band magnitude)    
bins=100
plt.figure('J'+cluster+' histogram '+str(bins)+'bins')    
plt.hist(np.abs(subarucatalog[:,3]-subarucatalog[:,4]),bins,color='red',label='zmin')
#plt.axvline(x=z,color='g')
#plt.xlim(0,5)
#plt.xlabel('zmin')