# -*- coding: utf-8 -*-
"""
Created on Thu May 24 22:50:17 2018

@author: zehaojin
"""
import functions_display as fun
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
DC_cut=0
####R<R_cut sources are accepted
R_cut=0.67
####Distance between sources(arcsec)>Jet_cut are accepted
Jet_cut=(7.121,8.1,15) ##"


shuffle_time=1000




ring=[4]

seperate_or_not=False
seperate_or_not=True

show_scatter=False

seperate_DC_IM=False

behind1_infrontof2_within3=1





###tangential shear scatter
for cluster in cluster_list:
    
    print 'J',cluster    
    
    ##get basic facts
    ra_center,dec_center,z,B_MAJ,B_MIN,B_PA = fun.facts(cluster)
    
    ##load catalog
    catalog = fun.loadcatalog(cluster) 
    ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
    entries_full = catalog.size
    
    
    catalog,subaru_zb,subaru_zml = fun.z_cut(cluster,catalog,z,behind1_infrontof2_within3)
    
    entries_z = catalog.size
    
    ##Jet cut
    #catalog,cut_list = fun.Jet_cutter(catalog,Jet_cut,ra_center,dec_center,cluster,subaru_zb,subaru_zml)
    entries_Jet = catalog.size
    
    
    
    
    
    #######find Jet cut
    #r_list=fun.Look_for_Jet_cut(catalog,ra_center,dec_center,cluster,Jet_cut)    
    
    #'''
    #position cut
    #catalog = fun.position_cut(catalog,cluster)
    entries_position = catalog.size
    #'''
    
    ##DC cut
    #catalog = fun.DC_cutter(catalog,DC_cut)
    entries_DC = catalog.size
    ##R cut
    catalog,R = fun.R_cutter(catalog,R_cut,B_MAJ,B_MIN)
    entries_R = catalog.size
    
    #find Jet cut
    #r_list=fun.Look_for_Jet_cut(catalog,ra_center,dec_center,cluster)
    
    
    #print 'sources selected:',entries_full,'->',entries_DC,'->',entries_DC_R,'sources'
    
    
    
    

    
    #print 'sources selected:',entries_full,'->',entries_z,'->',entries_Jet,'->',entries_position,'->',entries_DC,'->',entries_R,'sources'
    print 'sources selected:',entries_full,'->',entries_z,'->',entries_DC,'->',entries_R,'sources'
    #print 'Jet cut size=',cut_list.size
    
    #'''
    #position cut
    #catalog = fun.position_cut(catalog,cluster)
    #entries_position = catalog.size
    #print "sources after position cut",entries_position
    #'''    
    
    
    
    
    #'''
    ###get grid
    grid=fun.get_phi_and_grid(catalog,ra_center,dec_center)
    ###grid[(phi,x,y,dc_a,dc_b),i]
    
    
    
    #################DC calculatrions############################
    complex_shear_DC,real_avg,imag_avg=fun.get_complex_shear(catalog,True)
    ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
    print 'DC mean complex shear=',real_avg,'+',imag_avg,'i'    
    
    
    tangential_shear_DC,random_array=fun.shuffle_complex_shear_get_tangential_shear(complex_shear_DC,grid,shuffle_time,True)
    ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_b,6-dc_PA,7~end-shuffled_tangential_shear),i]
    
  
    ###################IM calculations#############################
    complex_shear_IM=fun.get_complex_shear(catalog,False)
    ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
    
    complex_shear_IM,real_avg,imag_avg=fun.apply_Bernstein_Jarvis_correction(complex_shear_IM,R,B_MAJ,B_MIN,B_PA)
    ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
    print 'IM mean complex shear=',real_avg,'+',imag_avg,'i' 
    
    tangential_shear_IM,random_array=fun.shuffle_complex_shear_get_tangential_shear(complex_shear_IM,grid,shuffle_time,False,random_array)
    ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_b,6-dc_PA,7~end-shuffled_tangential_shear),i]
    
    
    
    
    ##save data to fits file
    tangential_shear_DC_hdu=pyfits.PrimaryHDU(tangential_shear_DC)
    tangential_shear_IM_hdu=pyfits.PrimaryHDU(tangential_shear_IM)
    
    tangential_shear_DC_hdu.writeto('fits/J%s_x[0]_y[1]_tangential_shear[2]_randomized[3]_DC.fits' %cluster,clobber=True)
    tangential_shear_IM_hdu.writeto('fits/J%s_x[0]_y[1]_tangential_shear[2]_randomized[3]_IM.fits' %cluster,clobber=True)
    
    
    
    #fun.area_tangential_shear(catalog,tangential_shear_IM,cluster)

'''
####four quardants plot
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
    
    result_DC_1,result_DC_2,result_DC_3,result_DC_4=fun.seperate_quardants_plot(cluster,tangential_shear_DC,tangential_shear_IM,ring,shuffle_time)
    
#'''











    




#   '''
###averaged tangential shear
if seperate_or_not==True:
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
        
        ###averaged_tangential_shear(cluster,data_DC,data_IM,seperate_DC_IM,ring,shuffle_time,show_scatter)
        fun.averaged_tangential_shear(cluster,tangential_shear_DC,tangential_shear_IM,seperate_DC_IM,ring,shuffle_time,show_scatter,R_cut)
      
if seperate_or_not==False:
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
        
    fun.averaged_tangential_shear('0416+0717+1149',tangential_shear_DC_total,tangential_shear_IM_total,seperate_DC_IM,ring,shuffle_time,show_scatter,R_cut)
    '''
    ###comparison to http://iopscience.iop.org/article/10.1088/0004-637X/795/2/163/pdf page 16
    #plt.figure('J0416+0717+1149averaged tangential shear for 4 rings')    
    avg_z=np.average((0.396,0.548,0.544),weights=(144,208,116))###=0.500
    ##http://arcsec2parsec.joseonorbe.com/index.html  gives 1 arcmin=378.9166 kiloparsecs
    #1kpc=1/378.9166 arcmin
    kpc2arcmin=1/378.9166
    x_clustercentric_R=np.array([220.,300.,390.,500.,650.,840.,1050.,1350.,1800.,2300.,3050.])
    x_clustercentric_R*=kpc2arcmin
    y1_optical_tangential_shear=np.array([1.3e-1,1.5e-1,8.1e-2,9.1e-2,6.8e-2,5e-2,4.2e-2,3.1e-2,2.5e-2,1.5e-2,1.3e-2])
    y2_optical_cross_shear=np.array([0.035,-0.019,0.001,-0.002,0,0,-0.002,-0.006,0,-0.002,0])
    plt.scatter(x_clustercentric_R,y1_optical_tangential_shear,s=80,marker='*',color='black',label='optical tangential shear')
    plt.scatter(x_clustercentric_R,y2_optical_cross_shear,s=80,marker='x',color='green',label='optical cross shear')
    plt.legend(loc=0)
    plt.show()
#'''    
        


        