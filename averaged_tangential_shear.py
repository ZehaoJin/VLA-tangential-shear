# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 21:42:03 2018

@author: zehaojin
"""

import pyfits

import numpy as np
import matplotlib.pyplot as plt


cluster=raw_input('which cluster?\n 1. J0416 (default)\n 2. J0717\n 3. J1149\n')
if cluster=="1" or cluster=="":
	cluster='0416'
if cluster=="2":
	cluster='0717'
if cluster=="3":
	cluster='1149'

'''
http://iopscience.iop.org/article/10.1088/0004-637X/795/2/163/pdf     page 13-shear to compare with
http://adsabs.harvard.edu/abs/2014MNRAS.442.1507G                     page 1530-shear for J0416
'''


 

###How many rings you want?
ring=[4]
#ring=[3,4]


###show image?
show=True
###shuffle time
shuffle_time=1000


##J0416 z=0.42
##RA:04h 16m 08.38s
##DEC:−24° 04′ 20.80″   from Zitrin et al (http://adsabs.harvard.edu/abs/2013ApJ...762L..30Z),http://iopscience.iop.org/article/10.1088/2041-8205/762/2/L30/pdf
if cluster=='0416':
    ra_center=4*15+16*(1.0/4)+8.38*(1.0/240)
    dec_center=-(24+4*(1.0/60)+20.80*(1.0/3600))
    print 'center of cluster: (',ra_center,',',dec_center,')'
    
##J0717 z=0.548
##RA 07:17:32.63
##DEC 37:44:59.7      from http://iopscience.iop.org/article/10.1088/0004-637X/795/2/163/pdf
if cluster=='0717':
    ra_center=7*15+17*(1.0/4)+32.63*(1.0/240)
    dec_center=37+44*(1.0/60)+59.7*(1.0/3600)
    print 'center of cluster: (',ra_center,',',dec_center,')'

##J1149 z=0.544
##RA 11:49:35.69
##DEC 22:23:54.6      from http://iopscience.iop.org/article/10.1088/0004-637X/795/2/163/pdf
if cluster=='1149':
    ra_center=11*15+49*(1.0/4)+35.69*(1.0/240)
    dec_center=22+23*(1.0/60)+54.6*(1.0/3600)
    print 'center of cluster: (',ra_center,',',dec_center,')'



 


###read tangential shear fits file
fname='fits/J%s_ra[0]_dec[1]_tangential_shear[2]_randomized[3].fits' %cluster
hdulist = pyfits.open(fname)
data=hdulist[0].data
###data(ra[0] dec[1] tangential_shear[3],[i])
hdulist.close()

entries=data[0].size

newrow=np.zeros((1,entries))
data=np.vstack([data,newrow])
###data(ra[0] dec[1] tangential_shear[2] randomized_shear[3:-2] radius^2[-1],[i])
data[-1]=((data[0]-ra_center)*np.cos(dec_center*np.pi/180))**2+(data[1]-dec_center)**2
##ra_position=(catalog[i][1]-ra_center)*np.cos(dec_center*np.pi/180)

##sort according to data[3]
data=data[:,np.argsort(data[-1])]


result=np.zeros((len(ring),max(ring)))
result2=np.zeros((shuffle_time,len(ring),max(ring)))
error_mean=np.zeros((len(ring),max(ring)))
error_stdev=np.zeros((len(ring),max(ring)))

count=-1
for ring_quantity in ring:
    xticks=np.array([])
    count+=1
    plt.figure(count*2)
    plt.title('tangential shear scatter with '+str(ring_quantity)+' rings')
    plt.gca().invert_xaxis()
    plt.scatter(data[0],data[1],c=data[2],cmap='gray')
    plt.plot(ra_center,dec_center,marker='+', markersize=1000, color="red")
    plt.colorbar()
    plt.xlabel('RA(deg)')
    plt.ylabel('DEC(deg)')
    point_in_each_ring=entries//ring_quantity
    for i in range(ring_quantity-1):
        for j in range(point_in_each_ring):
            index=i*point_in_each_ring+j
            result[count][i]+=data[2][index]
            #print index
        result[count][i]/=point_in_each_ring
        radius=(np.sqrt(data[-1][index])+np.sqrt(data[-1][index+1]))/2
        xticks=np.append(xticks,radius)
        circle=plt.Circle((ra_center,dec_center),radius,color='r',fill=False)
        ax=plt.gca()
        ax.add_patch(circle)
        
    for j in range(point_in_each_ring*(ring_quantity-1),entries):
        result[count][ring_quantity-1]+=data[2][j]
        #print j
    result[count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))


    ###
    for time in range(shuffle_time):
        for i in range(ring_quantity-1):
            for j in range(point_in_each_ring):
                index=i*point_in_each_ring+j
                result2[time][count][i]+=data[3+time][index]
            result2[time][count][i]/=point_in_each_ring
        
        for j in range(point_in_each_ring*(ring_quantity-1),entries):
            result2[time][count][ring_quantity-1]+=data[3+time][j]
        result2[time][count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))
    
    
    ###calculate mean and stdev
    '''
    for time in range(shuffle_time):
        error_mean[count]+=result2[time][count]
    error_mean[count]/=shuffle_time
    '''
    error_mean[count]=np.average(result2[:,count,:],0)


    '''
    for time in range(shuffle_time):
        for i in range(ring_quantity-1):
            for j in range(point_in_each_ring):
                index=i*point_in_each_ring+j
                error_stdev[count][i]+=(data[3+time][index]-error_mean[count][i])**2
                #error_stdev[count][i]+=(data[3+time][index]-result2[time][count][i])**2
        for j in range(point_in_each_ring*(ring_quantity-1),entries):
            error_stdev[count][ring_quantity-1]+=(data[3+time][j]-error_mean[count][ring_quantity-1])**2
            #error_stdev[count][ring_quantity-1]+=(data[3+time][j]-result2[time][count][ring_quantity-1])**2
    error_stdev[count][0:ring_quantity-1]=np.sqrt(error_stdev[count][0:ring_quantity-1]/(point_in_each_ring*shuffle_time))
    error_stdev[count][ring_quantity-1]=np.sqrt(error_stdev[count][ring_quantity-1]/((entries-point_in_each_ring*3)*shuffle_time))
    '''
    
    
    error_stdev[count]=np.std(result2[:,count,:],0)



    xticks=(xticks+np.insert(xticks,0,0)[0:ring_quantity-1])/2
    xticks=np.append(xticks,xticks[-1]+xticks[0])##add last tick
    xticks=xticks*60 ##deg to arcmin

    plt.figure(count*2+1)
    plt.xlabel('radius(arcmin)')
    plt.ylabel('averaged tangential shear')
    #plt.scatter(np.array(range(ring_quantity))+1,result[count][:ring_quantity], color="red",marker='s',label="actual data")
    plt.errorbar(xticks,result[count][:ring_quantity],error_stdev[count][:ring_quantity],ecolor="green",fmt='ro',label="actual data")

    plt.scatter(xticks,error_mean[count][:ring_quantity], color="blue",marker='s',label="shuffled data")
    plt.xticks(xticks)
    plt.title('$\overline{\gamma}$  for'+str(ring_quantity)+' rings')
    #plt.legend(loc='center left', bbox_to_anchor=(0.8, 1))
    plt.legend(loc=0)


    print 'averged tangential shear:\n',result
    print 'mean of randomized tangential shear:\n',error_mean
    print 'standard deviation(1/2length of whole error bar):\n',error_stdev
    print 'ratio of std: ',error_stdev/np.min(error_stdev)


plt.show()

if show==False:
    plt.close('all')
   
      
#randomize for error bars