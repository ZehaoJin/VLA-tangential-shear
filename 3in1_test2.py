# -*- coding: utf-8 -*-
"""
Created on Thu May 24 22:50:17 2018

@author: zehaojin
"""
import functions as fun
import pyfits

import numpy as np
import matplotlib.pyplot as plt


'''
ned
clash paper 1000hkpc~5arcmin

typical redshift distribution for faint radio sources
VLA primary beam(PB)+dN/dS

z>zcluster
z<zcluster
z~zcluster(0.1)
'''


''' 
constants settings
'''



cluster_list=('0416','0717','1149')


##NVSS

'''original
####DC_MIN>DC_cut are accepted
DC_cut=0.2
####R<R_cut are accepted
R_cut=0.7
'''

'''one extreme
####DC_MIN>DC_cut are accepted
DC_cut=0.1
####R<R_cut are accepted
R_cut=0.95
'''

'''another extreme
####DC_MIN>DC_cut are accepted
DC_cut=0.5
####R<R_cut are accepted
R_cut=0.45
'''

####DC_MIN>DC_cut sources are accepted
DC_cut=0.
####R<R_cut sources are accepted
R_cut=2.0/3
####Distance between sources(arcsec)>Jet_cut are accepted
Jet_cut=(7.121,8.1,15) ##"

#ellipticity_cutter=5
ellipticity_cutter_list=[10000,10,5,4,3,2,1.5]
ellipticity_cutter_list=[9,8,7,6,5,4,3,2]
ellipticity_cutter_list=[7,6,5.5,5,4.5,4,3]
#DC_a/DC_b


shuffle_time=1000




ring=[5]

seperate_or_not=False
#seperate_or_not=True

show_scatter=False

seperate_DC_IM=False

behind1_infrontof2_within3=1



ellipticity_move_forward=(0.5/len(ellipticity_cutter_list))/1.25
for ellipticity_counter,ellipticity_cutter in enumerate(ellipticity_cutter_list):
###tangential shear scatter
    for cluster in cluster_list:
    
        #print 'calculating J',cluster    
    
        ##get basic facts
        ra_center,dec_center,z,B_MAJ,B_MIN,B_PA = fun.facts(cluster)
        
        ##load catalog
        catalog = fun.loadcatalog(cluster) 
        ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
        entries_full = catalog.size
        
        
        catalog,subaru_zb,subaru_zml = fun.z_cut(cluster,catalog,z,behind1_infrontof2_within3)
        
        entries_z = catalog.size
    
    
    
        ##DC cut
        catalog = fun.DC_cutter(catalog,DC_cut)
        entries_DC = catalog.size
        
        catalog = fun.ellipticity_cut(catalog,ellipticity_cutter)
        entries_ellipticity = catalog.size
    
        ##R cut
        catalog,R = fun.R_cutter(catalog,R_cut,B_MAJ,B_MIN)
        entries_R = catalog.size
   
    
        #'''
        ###get grid
        grid=fun.get_phi_and_grid(catalog,ra_center,dec_center)
        ###grid[(phi,x,y,dc_a,dc_b),i]
    
    
    
        #################DC calculatrions############################
        complex_shear_DC,real_avg,imag_avg=fun.get_complex_shear(catalog,True)
        ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
        #print 'DC mean complex shear=',real_avg,'+',imag_avg,'i'    
    
    
        tangential_shear_DC,random_array=fun.shuffle_complex_shear_get_tangential_shear(complex_shear_DC,grid,shuffle_time,True)
        ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_b,6-dc_PA,7~end-shuffled_tangential_shear),i]
        
        
        ###################IM calculations#############################
        complex_shear_IM=fun.get_complex_shear(catalog,False)
        ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
    
        complex_shear_IM,real_avg,imag_avg=fun.apply_Bernstein_Jarvis_correction(complex_shear_IM,R,B_MAJ,B_MIN,B_PA)
        ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
        #print 'IM mean complex shear=',real_avg,'+',imag_avg,'i' 
    
        tangential_shear_IM,random_array=fun.shuffle_complex_shear_get_tangential_shear(complex_shear_IM,grid,shuffle_time,False,random_array)
        ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_b,6-dc_PA,7~end-shuffled_tangential_shear),i]
    
    
    
    
        ##save data to fits file
        tangential_shear_DC_hdu=pyfits.PrimaryHDU(tangential_shear_DC)
        tangential_shear_IM_hdu=pyfits.PrimaryHDU(tangential_shear_IM)
        
        tangential_shear_DC_hdu.writeto('fits/J%s_x[0]_y[1]_tangential_shear[2]_randomized[3]_DC.fits' %cluster,clobber=True)
        tangential_shear_IM_hdu.writeto('fits/J%s_x[0]_y[1]_tangential_shear[2]_randomized[3]_IM.fits' %cluster,clobber=True)
        
        
        ###averaged tangential shear

    
    tangential_shear_DC_total=np.zeros((7+shuffle_time,0))
    tangential_shear_IM_total=np.zeros((7+shuffle_time,0))
    for cluster in cluster_list:
            
        ###read tangential shear fits file
        fname_DC='fits/J%s_x[0]_y[1]_tangential_shear[2]_randomized[3]_DC.fits' %cluster
        hdulist = pyfits.open(fname_DC)
        tangential_shear_DC=hdulist[0].data
        hdulist.close()
        
        fname_IM='fits/J%s_x[0]_y[1]_tangential_shear[2]_randomized[3]_IM.fits' %cluster
        hdulist = pyfits.open(fname_IM)
        tangential_shear_IM=hdulist[0].data
        hdulist.close()
        ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_b,6-dc_PA,7~end-shuffled_tangential_shear),i]
        
        tangential_shear_DC_total=np.append(tangential_shear_DC_total,tangential_shear_DC,axis=1)
        tangential_shear_IM_total=np.append(tangential_shear_IM_total,tangential_shear_IM,axis=1)
            
    print ellipticity_counter,'. a/b_cut:',ellipticity_cutter   
    fun.fixed_bins_averaged_tangential_shear_MORE_FUN('0416+0717+1149',tangential_shear_DC_total,shuffle_time,ellipticity_counter,ellipticity_move_forward)
    

print 'ellipticity_cutter_list:',ellipticity_cutter_list
plt.title('ellipticity_cutter_list:'+str(ellipticity_cutter_list))
plt.show()




        