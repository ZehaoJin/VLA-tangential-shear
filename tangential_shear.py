# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 15:44:40 2018

@author: zehaojin
"""


import pyfits

import numpy as np
import matplotlib.pyplot as plt
import random


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
Kaiser Squares Broodhurt KSB
1.put cut b>0.2
2.ssds vla   (a^2+b^2)psf/(a^2+b^2)im
''' 



###show image?
show=True
###shuffle time
shuffle_time=1000




##load catalog
catalog=np.loadtxt('VLA-HFF_%s_compact_optical_rasort.txt' %cluster,dtype={'names': ('name', 'RA', 'DEC','IM_MAJ','IM_MIN','M_STAR','IM_PA','SIGMA_IM_MAJ','SIGMA_IM_MIN','ZBEST'),'formats': ('|S26',np.float,np.float,np.float,np.float,'|S26',np.float,np.float,np.float,'|S26')},usecols=(0,1,2,21,23,43,25,22,24,40))
#catalog=np.loadtxt('VLA-HFF_%s_compact_optical_rasort.txt' %cluster,dtype={'names': ('name', 'RA', 'DEC','IM_MAJ','IM_MIN','M_STAR','IM_PA','SIGMA_IM_MAJ','SIGMA_IM_MIN','ZBEST'),'formats': ('|S26',np.float,np.float,np.float,np.float,'|S26',np.float,np.float,np.float,'|S26')},usecols=(0,1,2,27,29,43,31,22,24,40))
ra,dec=np.loadtxt('VLA-HFF_%s_compact_optical_rasort.txt' %cluster,usecols=(1,2),unpack=True)
###ra,dec are for making plots
entries=catalog.size
if show==True:
    for i in range(entries):
        print catalog[i][3],catalog[i][7],catalog[i][4],catalog[i][8],catalog[i][9],catalog[i][7]/catalog[i][3],catalog[i][8]/catalog[i][4]
    print 'IM_MAJ,SIGMA_IM_MAJ,IM_MIN,SIGMA_IM_MIN,ZBEST,SIGMA/MAJ,SIGMA/MIN'
##catalog[i,(name_0,ra_1,dec_2,maj_3,min_4,mass_5,pa_6)]



########construct an array for results
result=np.zeros((5+shuffle_time,entries))
###[(0-phi,1-ellipticity,2-tangential shear,3-real[complex_shear],4-img[complex_shear],5~last-tangential_shear calculated after shuffle),i]



########Determine center of Mass
'''
total_mass=0.
ra_summation=0.
dec_summation=0.
for i in range(entries):
    if catalog[i][5]!='-':
        catalog[i][5]=np.float(catalog[i][5])
        total_mass+=float(catalog[i][5])
        ra_summation+=(float(catalog[i][5])*catalog[i][1])
        dec_summation+=(float(catalog[i][5])*catalog[i][2])

ra_center=ra_summation/total_mass
dec_center=dec_summation/total_mass

print 'center of mass: (',ra_center,',',dec_center,')' 
'''


##J0416 z=0.42
##RA:04h 16m 08.38s
##DEC:−24° 04′ 20.80″   from Zitrin et al (http://adsabs.harvard.edu/abs/2013ApJ...762L..30Z)
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
    
    
    


########Determine phi
for i in range(entries):
    ##catalog[i,(name_0,ra_1,dec_2,maj_3,min_4,mass_5,pa_6)]
    ra_position=catalog[i][1]-ra_center
    dec_position=catalog[i][2]-dec_center
    if ra_position==0 and dec_position==0:
        result[0][i]=0
    elif ra_position==0 and dec_position>0:
        result[0][i]=0
    elif ra_position==0 and dec_position<0:
        result[0][i]=np.pi
    elif ra_position>0 and dec_position==0:
        result[0][i]=np.pi/2
    elif ra_position<0 and dec_position==0:
        result[0][i]=np.pi*3/2
    elif ra_position>0 and dec_position>0:
        result[0][i]=np.arctan(float(ra_position)/dec_position)
    elif ra_position>0 and dec_position<0:
        result[0][i]=np.arctan(float(dec_position)/ra_position*-1)+(np.pi/2)
    elif ra_position<0 and dec_position<0:
        result[0][i]=np.arctan(float(ra_position)/dec_position)+(np.pi)
    elif ra_position<0 and dec_position>0:
        result[0][i]=np.arctan(float(dec_position)/ra_position*-1)+(np.pi*3/2)
    
    
    result[0][i]=(-result[0][i])+2*np.pi
    result[0][i]+=np.pi/2
    if result[0][i]>2*np.pi:
        result[0][i]-=2*np.pi
    


plt.figure('phi J%s' %cluster)
plt.gca().invert_xaxis()
plt.scatter(ra,dec,c=result[0],cmap='gray')
plt.plot(ra_center,dec_center,marker='+', markersize=1000, color="red")
plt.colorbar()
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')



########Determine ellipticity
for i in range(entries):
    ##catalog[i,(name_0,ra_1,dec_2,maj_3,min_4,mass_5,pa_6)]
    r=catalog[i][4]/catalog[i][3]
    result[1][i]=(1-r*r)/(1+r*r)

result[1]=np.nan_to_num(result[1])
    
plt.figure('ellipticity J%s' %cluster)
plt.gca().invert_xaxis()
plt.scatter(ra,dec,c=result[1],cmap='gray')
plt.plot(ra_center,dec_center,marker='+', markersize=1000, color="red")
plt.colorbar()
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')

    
    
########Determine tangential shear  https://www.astro.umd.edu/~richard/ASTR680/Schneider_weak_lensing.pdf
real_part_sum=0
imag_part_sum=0 
for i in range(entries):
    ##catalog[i,(name_0,ra_1,dec_2,maj_3,min_4,mass_5,pa_6)]   
    complex_shear=result[1][i]*np.exp(2j*catalog[i][6]*np.pi/180)
    ##result[(0-phi,1-ellipticity,2-tangential shear,3-real[complex_shear],4-imag[complex_shear]),i]    
    result[3][i]=np.real(complex_shear)
    result[4][i]=np.imag(complex_shear)
    
    real_part_sum+=result[3][i]
    imag_part_sum+=result[4][i]
    
real_noise=real_part_sum/entries
imag_noise=imag_part_sum/entries
print 'noise: ',real_noise,'+',imag_noise,'i'


'''

##for test file, might need to ## these sentences below
result[3]=(result[3]-real_noise)*entries/(entries-1)
result[4]=(result[4]-imag_noise)*entries/(entries-1)

'''

###randomize stuff and get tangential shear(randomized)
random_array=np.zeros((shuffle_time,entries))
for i in range(shuffle_time):
    random_list=list(range(entries))
    for j in range(10):
            random.shuffle(random_list)
    random_array[i]=random_list

random_array=random_array.astype(int)

for j in range(shuffle_time):
    for i in range(entries):
###[(0-phi,1-ellipticity,2-tangential shear,3-real[complex_shear],4-img[complex_shear],5~end-tangential_shear calculated after shuffle),i]
        complex_shear=result[3][random_array[j][i]]+1j*result[4][random_array[j][i]]
        result[5+j][i]=-np.real(complex_shear*np.exp(-2j*result[0][i]))



###calculate tangential shear
for i in range(entries):
    complex_shear=result[3][i]+1j*result[4][i]    
    result[2][i]=-np.real(complex_shear*np.exp(-2j*result[0][i]))
   

plt.figure('Tangential Shear J%s' %cluster)
plt.gca().invert_xaxis()
plt.scatter(ra,dec,c=result[2],cmap='gray')
plt.plot(ra_center,dec_center,marker='+', markersize=1000, color="red")
plt.colorbar()
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')

plt.show()
if show==False:
    plt.close('all')


##save data to fits file
data=np.zeros((3+shuffle_time,entries))
data[0]=ra
data[1]=dec
data[2]=result[2]
data[3:3+shuffle_time]=result[5:5+shuffle_time]
tangential_shear_hdu=pyfits.PrimaryHDU(data)
tangential_shear_hdu.writeto('fits/J%s_ra[0]_dec[1]_tangential_shear[2]_randomized[3].fits' %cluster,clobber=True)




##5/20 6 last2mon 7 oneweek














'''    
plt.figure('Shear')
plt.gca().invert_xaxis()
plt.scatter(ra,dec,c=result[1],cmap='gray')
plt.plot(ra_center,dec_center,marker='+', markersize=1000, color="red")
plt.colorbar()
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')
'''    

'''
########Determine tangential shear
for i in range(entries):
    
    complex_shear=result[1][i]*np.exp(2j*result[0][i])
    shear_1=np.real(complex_shear)
    shear_2=np.imag(complex_shear)
    result[2][i]=-((shear_1*np.cos(2*result[0][i]))+(shear_2*np.sin(2*result[0][i])))


    result[2][i]=-np.real(result[1][i]*np.exp(-2j*result[0][i]))

    
plt.figure('Tangential Shear')
plt.gca().invert_xaxis()
plt.scatter(ra,dec,c=result[2],cmap='gray')
plt.plot(ra_center,dec_center,marker='+', markersize=1000, color="red")
plt.colorbar()
plt.xlabel('RA(deg)')
plt.ylabel('DEC(deg)')
'''