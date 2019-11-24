# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:19:47 2018

@author: zehaojin
"""

import pyfits

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random






def facts(cluster):
    ##J0416 z=0.42
    ##RA:04h 16m 08.38s
    ##DEC:−24° 04′ 20.80″   from Zitrin et al (http://adsabs.harvard.edu/abs/2013ApJ...762L..30Z),http://iopscience.iop.org/article/10.1088/2041-8205/762/2/L30/pdf
    if cluster=='0416':
        ra_center=4*15+16*(1.0/4)+8.38*(1.0/240)
        dec_center=-(24+4*(1.0/60)+20.80*(1.0/3600))
        
        z=0.396
        
        B_MAJ    = 0.000261401928209888*3600                                                 
        B_MIN    = 0.000140935584386377*3600                                                 
        B_PA     =     1.92392664856767
    
    ##J0717 z=0.548
    ##RA 07:17:32.63
    ##DEC 37:44:59.7      from http://iopscience.iop.org/article/10.1088/0004-637X/795/2/163/pdf
    if cluster=='0717':
        ra_center=7*15+17*(1.0/4)+32.63*(1.0/240)
        dec_center=37+44*(1.0/60)+59.7*(1.0/3600)
        
        z=0.548
        
        B_MAJ    = 0.000201486997608722*3600                                                 
        B_MIN    = 0.000170685105066699*3600                                                
        B_PA     =     93.5441097813981


    ##J1149 z=0.544
    ##RA 11:49:35.69
    ##DEC 22:23:54.6      from http://iopscience.iop.org/article/10.1088/0004-637X/795/2/163/pdf
    if cluster=='1149':
        ra_center=11*15+49*(1.0/4)+35.69*(1.0/240)
        dec_center=22+23*(1.0/60)+54.6*(1.0/3600)
        
        z=0.544
        
        B_MAJ    = 0.000142663419847312*3600                                                 
        B_MIN    = 0.000134522499313709*3600                                                 
        B_PA     =     35.8609177934992
        
    return ra_center,dec_center,z,B_MAJ,B_MIN,B_PA







def loadcatalog(cluster):
    ##load catalog
    
    catalog=np.loadtxt('VLA-HFF_%s_compact_optical_rasort.txt' %cluster,dtype={'names': ('name', 'RA', 'DEC','IM_MAJ','IM_MIN','IM_PA','DC_MAJ','DC_MIN','DC_PA','SUBARU_ID','SUBARU_SEP','ZBPZ','ZBPZ_LOWER'),'formats': ('|S26',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,'|S26','|S26','|S26','|S26')},usecols=(0,1,2,21,23,25,27,29,31,50,53,56,57))    
      
    ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8,subaru_id_9,subaru_sep_10,zbpz_11,zbpz_lower_12)]
    '''  
    catalog=np.loadtxt('VLA-HFF_%s_compact_rasort.txt' %cluster,dtype={'names': ('name', 'RA', 'DEC','IM_MAJ','IM_MIN','IM_PA','DC_MAJ','DC_MIN','DC_PA'),'formats': ('|S26',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float)},usecols=(0,1,2,21,23,25,27,29,31))
    '''
    return catalog




def z_cut(cluster,catalog,z,behind1_infrontof2_within3):
    
    z_within_sigma=0.3   
    
    ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8,subaru_id_9)]
    if cluster=='0416':
        subarucatalog=np.loadtxt('hlsp_clash_subaru_suprimecam_macs%s_photoz-cat.txt' %cluster,usecols=(0,12,13,17)) #(0_id,1_zb,2_zmin,3_zml)
    if cluster=='0717':
        subarucatalog=np.loadtxt('hlsp_clash_subaru_suprimecam_macs%s_photoz-cat.txt' %cluster,usecols=(0,20,21,25))
    if cluster=='1149':
        subarucatalog=np.loadtxt('hlsp_clash_subaru_suprimecam_macs%s_photoz-cat.txt' %cluster,usecols=(0,16,17,21))
    #subarucols=(0,12,11,13,14,17,19,10)
    #(0_id,12_zbest,11_uncertainty in z,13_zbmin,14_zbmax,17_zml,19_Poorness of BPZ fit,10_z_Subaru band magnitude)
    entries=catalog.size
    cut_list=np.array([])
    subaru_zb=np.array([])
    subaru_zmin=np.array([])
    subaru_zml=np.array([])


    
    for i in range(entries):
        if catalog[i][9]=='-':
            subaru_zb=np.append(subaru_zb,0)
            subaru_zmin=np.append(subaru_zmin,0)
            subaru_zml=np.append(subaru_zml,0)
        else:
            subaru_index=np.int(catalog[i][9])-1
            subaru_zb=np.append(subaru_zb,subarucatalog[subaru_index][1])
            subaru_zmin=np.append(subaru_zmin,subarucatalog[subaru_index][2])
            subaru_zml=np.append(subaru_zml,subarucatalog[subaru_index][3])
        
        #if subaru_zmin[i]<z and subaru_zb[i]<z and subaru_zml[i]<z:
        #if subaru_zml[i]<z:
        #if behind1_infrontof2_within3==1 and subaru_zb[i]<=(z+z_within_sigma):
        if behind1_infrontof2_within3==1 and subaru_zmin[i]<=(z):
        #if subaru_zb[i]<z or subaru_zml[i]<z:
        #if subaru_zb[i]<z and subaru_zml[i]<z:
        #if (subaru_zmin[i])<z:
            cut_list=np.append(cut_list,i)
        if behind1_infrontof2_within3==2 and (subaru_zb[i])>=(z-z_within_sigma):
        #if behind1_infrontof2_within3==2 and (subaru_zb[i])>=(z):
            cut_list=np.append(cut_list,i)
        if behind1_infrontof2_within3==3 and (subaru_zb[i]>(z+z_within_sigma) or subaru_zb[i]<(z-z_within_sigma)):
            cut_list=np.append(cut_list,i)


        
    catalog=np.delete(catalog,cut_list,0)
    #R=np.delete(R,cut_list,0)
    '''
    bins=100
    plt.figure('J'+cluster+' zmin histogram '+str(bins)+'bins')    
    plt.hist(subaru_zmin,bins,color='red',label='zmin')
    plt.axvline(x=z,color='g')
    plt.xlim(0,5)
    plt.xlabel('zmin')
    
    plt.figure('J'+cluster+' zb histogram '+str(bins)+'bins')    
    plt.hist(subaru_zb,bins,color='red',label='zb')
    plt.axvline(x=z,color='g')
    plt.xlim(0,5)
    plt.xlabel('zb')
    
    plt.figure('J'+cluster+' zml histogram '+str(bins)+'bins')
    plt.hist(subaru_zml,bins,color='blue',label='zml')
    plt.axvline(x=z,color='g')
    plt.xlim(0,5)
    plt.xlabel('zml')
    
    plt.show()
    '''
    
    '''
    for i in range(entries):
        print cluster,subaru_zb[i],subaru_zmin[i],subaru_zml[i]
    '''
        
    return catalog,subaru_zb,subaru_zml















def position_cut(catalog,cluster):
    ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
    '''
    J0717: https://arxiv.org/pdf/0905.3650.pdf page3  
                ra1:7:17:35,  7:17:36.5    dec1:37:45:45,37:46:15   
                ra2:7:17:30.7,7:17:32.2    dec2:37:44:45,37:45:00
    J1149: https://arxiv.org/pdf/1608.01329.pdf page2  (main)
           http://www.mergingclustercollaboration.org/macs-j114952223.html
                ra1:11:49:43  ,11:49:47     dec1:22:21:00,22:22:08
                ra2:11:49:21.5,11:49:24.5   dec2:22:22:54,22:24:06
    J0416: no cuts are done
    '''
        
        
    J0717_ra_cut1=np.array([7*15+17*(1.0/4)+35*(1.0/240),7*15+17*(1.0/4)+36.5*(1.0/240)])
    J0717_dec_cut1=np.array([37+45*(1.0/60)+45*(1.0/3600),37+46*(1.0/60)+15*(1.0/3600)])
    J0717_ra_cut2=np.array([7*15+17*(1.0/4)+30.7*(1.0/240),7*15+17*(1.0/4)+32.2*(1.0/240)])
    J0717_dec_cut2=np.array([37+44*(1.0/60)+45*(1.0/3600),37+45*(1.0/60)+0*(1.0/3600)])
    
    J1149_ra_cut1=np.array([11*15+49*(1.0/4)+43*(1.0/240),11*15+49*(1.0/4)+47*(1.0/240)])
    J1149_dec_cut1=np.array([22+21*(1.0/60)+0*(1.0/3600),22+22*(1.0/60)+8*(1.0/3600)])
    J1149_ra_cut2=np.array([11*15+49*(1.0/4)+21.5*(1.0/240),11*15+49*(1.0/4)+24.5*(1.0/240)])
    J1149_dec_cut2=np.array([22+22*(1.0/60)+54*(1.0/3600),22+24*(1.0/60)+6*(1.0/3600)])
    

    
    entries=catalog.size
    cut_list=np.array([])
    if cluster=='0717':
        for i in range(entries):
            if (J0717_ra_cut1[0] <= catalog[i][1] <= J0717_ra_cut1[1]) and (J0717_dec_cut1[0] <= catalog[i][2] <= J0717_dec_cut1[1]):
                cut_list=np.append(cut_list,i)
            if (J0717_ra_cut2[0] <= catalog[i][1] <= J0717_ra_cut2[1]) and (J0717_dec_cut2[0] <= catalog[i][2] <= J0717_dec_cut2[1]):
                cut_list=np.append(cut_list,i)
        
        catalog=np.delete(catalog,cut_list,0)
        #print 'position cut',cut_list

        
        
    if cluster=='1149':
        for i in range(entries):
            if (J1149_ra_cut1[0] <= catalog[i][1] <= J1149_ra_cut1[1]) and (J1149_dec_cut1[0] <= catalog[i][2] <= J1149_dec_cut1[1]):
                cut_list=np.append(cut_list,i)
            if (J1149_ra_cut2[0] <= catalog[i][1] <= J1149_ra_cut2[1]) and (J1149_dec_cut2[0] <= catalog[i][2] <= J1149_dec_cut2[1]):
                cut_list=np.append(cut_list,i)
        
        catalog=np.delete(catalog,cut_list,0)
        #print 'position cut',cut_list
        
        
        
    return catalog




def ellipticity_cut(catalog,ellipticity_cutter):
    #catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
    entries=catalog.size
    cut_list=np.array([])
    for i in range(entries):
        if (catalog[i][6]/catalog[i][7])>=ellipticity_cutter:
            #print catalog[i][7],catalog[i][6]
            cut_list=np.append(cut_list,i)
    
    catalog=np.delete(catalog,cut_list,0)
    
    return catalog









    
def DC_cutter(catalog,DC_cut):
    '''
    definition of DC and IM datas: http://sundog.stsci.edu/first/catalogs/readme.html
    '''
    ####DC_MIN>DC_cut are accepted
    catalog1=catalog
    c=0
    entries=catalog.size
    for i in range(entries):
        if catalog[i][7]<=DC_cut:
            catalog1=np.delete(catalog1,i-c,0)
            c+=1
    return catalog1
    
    

def R_cutter(catalog,R_cut,B_MAJ,B_MIN):
    '''
    see functions.apply_Bernstein_Jarvis_correction for reference
    '''
    entries=catalog.size
    R=np.zeros(entries)
    for i in range(entries):
        ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
        R[i]=(B_MAJ**2+B_MIN**2)/(catalog[i][3]**2+catalog[i][4]**2)
        
    R1=R
    c=0
    for i in range(entries):
        if R[i]>=R_cut:
            catalog=np.delete(catalog,i-c,0)
            R1=np.delete(R1,i-c,0)
            c+=1

    return catalog,R1
    
    
def Look_for_Jet_cut(catalog,ra_center,dec_center,cluster,Jet_cut):
    entries=catalog.size
    cluster_list=('0416','0717','1149')
    if cluster==cluster_list[0]:
        Jet_cut=Jet_cut[0]
    if cluster==cluster_list[1]:
        Jet_cut=Jet_cut[1]
    if cluster==cluster_list[2]:
        Jet_cut=Jet_cut[2]
    
    
    bins=100
    
    
    r_list=np.array([])
    factor=np.cos(dec_center*np.pi/180)
    for i in range(entries-1):
        ra1=catalog[i][1]
        dec1=catalog[i][2]
        for j in range(i,entries):
            if i!=j:
                ra2=catalog[j][1]
                dec2=catalog[j][2]
                r=np.sqrt(((ra1-ra2)*factor)**2+(dec1-dec2)**2)
                r_list=np.append(r_list,r)
                
    r_list*=3600
    plt.figure('J'+cluster+' r histogram '+str(bins)+'bins')    
    plt.hist(r_list,bins,(0,300))
    plt.xlabel('arcsec')
    plt.axvline(x=Jet_cut,color='r')



    return r_list









    
def Jet_cutter(catalog,Jet_cut,ra_center,dec_center,cluster,subaru_zb,subaru_zml):
    '''
    The intention would be to cut both sources when they are too close in position. 
    The idea is that adjacent sources may in fact be a single radio source with a compact core and jets coming out. 
    The jet alignment would be something strongly anisotropic, 
    in a way that is completely unrelated to the gravitational effect of the galaxy cluster.
    '''
    cluster_list=('0416','0717','1149')
    
    z_limit=1    
    
    entries=catalog.size
    if cluster==cluster_list[0]:
        Jet_cut=Jet_cut[0]
    if cluster==cluster_list[1]:
        Jet_cut=Jet_cut[1]
    if cluster==cluster_list[2]:
        Jet_cut=Jet_cut[2]
        
        
    Jet_cut=Jet_cut/3600.0 ##arcsec to deg
    Jet_cut=Jet_cut**2     ##deg^2, r^2
    
    
    cut_list=np.array([])
    factor=np.cos(dec_center*np.pi/180)
    for i in range(entries-1):
        ra1=catalog[i][1]
        dec1=catalog[i][2]
        for j in range(i,entries):
            if i!=j:
                ra2=catalog[j][1]
                dec2=catalog[j][2]
                r=((ra1-ra2)*factor)**2+(dec1-dec2)**2
                if r<Jet_cut:
                    if np.abs(subaru_zb[i]-subaru_zb[j])<z_limit or np.abs(subaru_zml[i]-subaru_zml[j])<z_limit:
                        cut_list=np.append(cut_list,i)
                        cut_list=np.append(cut_list,j)
    
    '''            
    cut_list=np.array([])
    factor=np.cos(dec_center*np.pi/180)
    for i in range(entries-1):
        r_list=np.array([])
        j_list=np.array([])
        ra1=catalog[i][1]
        dec1=catalog[i][2]
        for j in range(i,entries):
            if i!=j:
                ra2=catalog[j][1]
                dec2=catalog[j][2]
                r=((ra1-ra2)*factor)**2+(dec1-dec2)**2
                r_list=np.append(r_list,r)
                j_list=np.append(j_list,j)
        if np.min(r_list)<Jet_cut:
            cut_list=np.append(cut_list,i)
            cut_list=np.append(cut_list,j_list[np.argmin(r_list)])
    '''
                    
    cut_list=np.unique(cut_list)
    cut_list=np.sort(cut_list)
    
    c=0
    for i in cut_list:        
        catalog=np.delete(catalog,i-c,0)
        c+=1
                    
            
            
        
    
    
    
    
    return catalog,cut_list
    
    
    
def get_phi_and_grid(catalog,ra_center,dec_center):
    '''
    change in x=change in RA cos(dec or dec of center)
    change in y=change in DEC
    definition of phi: https://www.astro.umd.edu/~richard/ASTR680/Schneider_weak_lensing.pdf
    '''
    entries=catalog.size
    result=np.zeros((6,entries))
            
    
    
    
    
    ########Determine phi
    for i in range(entries):
        ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
    
        ra_position=(catalog[i][1]-ra_center)*np.cos(dec_center*np.pi/180)
        dec_position=catalog[i][2]-dec_center
        
        result[1][i]=ra_position
        result[2][i]=dec_position
        result[3][i]=catalog[i][6]
        result[4][i]=catalog[i][7]
        result[5][i]=catalog[i][8]
        
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
        
    return result   ###grid[(phi,x,y,dc_a,dc_b,dc_PA),i]
    

def get_complex_shear(catalog,DC_or_not):
    '''
    https://www.astro.umd.edu/~richard/ASTR680/Schneider_weak_lensing.pdf
    '''
    entries=catalog.size
    ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
    if DC_or_not==False:###IM
        a_index=3
        b_index=4
        pa_index=5
    if DC_or_not==True:###DC
        a_index=6
        b_index=7
        pa_index=8
    
    result=np.zeros((3,entries))
    ###result[(0-ellipcity,e1-1,e2-2),i]
    ########Determine ellipticity
    for i in range(entries):
        r=catalog[i][b_index]/catalog[i][a_index]
        result[0][i]=(1-r*r)/(1+r*r)

    result[0]=np.nan_to_num(result[0])
    
    real_part_sum=0
    imag_part_sum=0 
    for i in range(entries):
        ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]  
        complex_shear=result[0][i]*np.exp(2j*np.deg2rad(90-catalog[i][pa_index]))
        complex_shear/=2
        result[1][i]=np.real(complex_shear)
        result[2][i]=np.imag(complex_shear)
    
        real_part_sum+=result[1][i]
        imag_part_sum+=result[2][i]
    
    real_avg=real_part_sum/entries
    imag_avg=imag_part_sum/entries
    
    if DC_or_not==False:###IM
        return result
    if DC_or_not==True:###DC
        return result,real_avg,imag_avg
    ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
        




def shuffle_complex_shear_get_tangential_shear(complex_shear,grid,shuffle_time,DC_or_not,random_array=0):
    '''
    formulaism for tangential shear: https://www.astro.umd.edu/~richard/ASTR680/Schneider_weak_lensing.pdf
    note the thing calculated is X. for shear, shear=X/2
    '''
    ###randomize stuff and get tangential shear(randomized)
    entries=complex_shear[0].size
    tangential_shear=np.zeros((7+shuffle_time,entries))
    ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_a,6-dc_PA,7~end-shuffled_tangential_shear),i]
    ###complex_shear[(0-ellipcity,e1-1,e2-2),i]
    ###grid[(phi,x,y,dc_a,dc_b,dc_PA),i]
    
    tangential_shear[0]=grid[1]
    tangential_shear[1]=grid[2]
    tangential_shear[4]=grid[3]
    tangential_shear[5]=grid[4]
    tangential_shear[6]=grid[5]
    
    ##tangential shear
    for i in range(entries):
        complex_shear_complex=complex_shear[1][i]+1j*complex_shear[2][i]    
        tangential_shear[2][i]=-np.real(complex_shear_complex*np.exp(-2j*grid[0][i]))
        tangential_shear[3][i]=np.imag(complex_shear_complex*np.exp(-2j*grid[0][i]))
    
    
    
    ##randomized tangential shear
    if DC_or_not==True:
        random_array=np.zeros((shuffle_time,entries))
        for i in range(shuffle_time):
            random_list=list(range(entries))
            for j in range(10):
                    random.shuffle(random_list)
                    random_array[i]=random_list
        random_array=random_array.astype(int)

    for j in range(shuffle_time):
        for i in range(entries):
            complex_shear_complex=complex_shear[1][random_array[j][i]]+1j*complex_shear[2][random_array[j][i]]
            tangential_shear[7+j][i]=-np.real(complex_shear_complex*np.exp(-2j*grid[0][i]))
    return tangential_shear,random_array
    
def apply_Bernstein_Jarvis_correction(complex_shear,R,B_MAJ,B_MIN,B_PA):
    '''
    Bernstein & Jarvis (2002)     https://arxiv.org/pdf/1507.05977.pdf    http://iopscience.iop.org/article/10.1086/338085/pdf
    e1 = (e1meas - R*e1psf)/(1-R)
    e2 = (e2meas - R*e2psf)/(1-R)
    R = mrr_cc_psf/mrr_cc*(4/mcr4_psf-1)/(4/mcr4-1)=(a^2+b^2)psf/(a^2+b^2)im
    '''
    #psf
    r=B_MIN/B_MAJ
    complex_shear_psf=((1-r*r)/(1+r*r))*np.exp(2j*np.deg2rad(90-B_PA))
    complex_shear_psf/=2
    e1_psf=np.real(complex_shear_psf)
    e2_psf=np.imag(complex_shear_psf)
    
    #print 'psf_shear: ',e1_psf,'+',e2_psf,'i'
    
    ###complex_shear[(0-ellipcity,e1-1,e2-2),i]    
    complex_shear[1]=(complex_shear[1]-R*e1_psf)/(1-R)
    complex_shear[2]=(complex_shear[2]-R*e2_psf)/(1-R)
    
    real_avg=np.average(complex_shear[1])
    imag_avg=np.average(complex_shear[2])
    
   
    return complex_shear,real_avg,imag_avg



def averaged_tangential_shear(cluster,data_DC,data_IM,seperate_DC_IM,ring,shuffle_time,show_scatter):
    '''
    https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot
    '''
    J0717_ra_cut1=(7*15+17*(1.0/4)+35*(1.0/240),7*15+17*(1.0/4)+36.5*(1.0/240))
    J0717_dec_cut1=(37+45*(1.0/60)+45*(1.0/3600),37+46*(1.0/60)+15*(1.0/3600))
    J0717_ra_cut2=(7*15+17*(1.0/4)+30.7*(1.0/240),7*15+17*(1.0/4)+32.2*(1.0/240))
    J0717_dec_cut2=(37+44*(1.0/60)+45*(1.0/3600),37+45*(1.0/60)+0*(1.0/3600))
    
    J1149_ra_cut1=(11*15+49*(1.0/4)+43*(1.0/240),11*15+49*(1.0/4)+47*(1.0/240))
    J1149_dec_cut1=(22+21*(1.0/60)+0*(1.0/3600),22+22*(1.0/60)+8*(1.0/3600))
    J1149_ra_cut2=(11*15+49*(1.0/4)+21.5*(1.0/240),11*15+49*(1.0/4)+24.5*(1.0/240))
    J1149_dec_cut2=(22+22*(1.0/60)+54*(1.0/3600),22+24*(1.0/60)+6*(1.0/3600))
    
    

    
    data_fmt_DC='ro'  ##red circles
    error_bar_color_DC='g' ##green
    shuffled_data_color_DC='c' ##cyan

    data_fmt_IM='yo' ##yellow circles
    error_bar_color_IM='b' ##blue
    shuffled_data_color_IM='m'##magenta    
    
    entries=data_DC[0].size
    
    ###tangential_shear[(0-x,1-y,2-tangential_shear,3-cross_shear,4-dc_a,5-dc_b,6-dc_PA,7~end-shuffled_tangential_shear),i]
    ellipses=[]
    for i in range(entries):
        ellipses.append(patches.Ellipse((data_DC[0][i],data_DC[1][i]),data_DC[5][i]/3600*30,data_DC[4][i]/3600*30,-data_DC[6][i],edgecolor='black',fill=False))
    
        
    
    newrow=np.zeros((1,entries))
    data_DC=np.vstack([data_DC,newrow])
    data_IM=np.vstack([data_IM,newrow])
    ###old data(ra[0] dec[1] tangential_shear[2] randomized_shear[3:-2] radius^2[-1],[i])
    ###new data=tangential_shear(x[0] y[1] tangential_shear[2] cross_shear[3] randomized_shear[4:-2] radius^2[-1],[i])
    data_DC[-1]=data_DC[0]**2+data_DC[1]**2
    data_IM[-1]=data_IM[0]**2+data_IM[1]**2
    ##sort according to data[-1] (radius^2)
    data_DC=data_DC[:,np.argsort(data_DC[-1])]
    data_IM=data_IM[:,np.argsort(data_IM[-1])]


    result_DC=np.zeros((len(ring),max(ring)))
    result2_DC=np.zeros((shuffle_time,len(ring),max(ring)))
    result3_DC=np.zeros((len(ring),max(ring)))
    error_mean_DC=np.zeros((len(ring),max(ring)))
    error_stdev_DC=np.zeros((len(ring),max(ring)))
    
    result_IM=np.zeros((len(ring),max(ring)))
    result2_IM=np.zeros((shuffle_time,len(ring),max(ring)))
    result3_IM=np.zeros((len(ring),max(ring)))
    error_mean_IM=np.zeros((len(ring),max(ring)))
    error_stdev_IM=np.zeros((len(ring),max(ring)))
    
    count=-1
    for ring_quantity in ring:
        
        count+=1
        
        ###################DC
        xticks=np.array([])
        if show_scatter==True:
            plt.figure('J'+cluster+' tangential shear scatter with '+str(ring_quantity)+' rings(DC)')
            plt.title('J'+cluster+'$\gamma $ with'+str(ring_quantity)+' rings(DC)')
            plt.xlim(-0.20,0.25)
            plt.ylim(-0.20,0.20)
            plt.gca().invert_xaxis()
            ax=plt.gca()
            plt.scatter(data_DC[0],data_DC[1],c=data_DC[2],cmap='rainbow',edgecolors='none')
            scatter_xticks=np.linspace(-0.20,0.25,20)
            scatter_xticks=np.append(scatter_xticks,0)
            scatter_xticks=np.append(scatter_xticks,0)
            scatter_yticks=np.linspace(-0.20,0.20,20)
            plt.xticks(scatter_xticks,np.round(scatter_xticks*60,2))
            plt.yticks(scatter_yticks,np.round(scatter_yticks*60,2))
            for ellipse in ellipses:
                ax.add_patch(ellipse)
                
            plt.plot(0,0,marker='+', markersize=1000, color="red")
            plt.colorbar()
            plt.xlabel("$\Delta x(arcmin)$")
            plt.ylabel("$\Delta y(arcmin)$")
            

   
        point_in_each_ring=entries//ring_quantity
        for i in range(ring_quantity-1):
            for j in range(point_in_each_ring):
                index=i*point_in_each_ring+j
                result_DC[count][i]+=data_DC[2][index]
                result3_DC[count][i]+=data_DC[3][index]
                #print index
            result_DC[count][i]/=point_in_each_ring
            result3_DC[count][i]/=point_in_each_ring
            radius=(np.sqrt(data_DC[-1][index])+np.sqrt(data_DC[-1][index+1]))/2
            print cluster,'bin radius(arcmin) ',radius*60
            xticks=np.append(xticks,radius)
            if show_scatter==True:
                circle=plt.Circle((0,0),radius,color='r',fill=False)
                ax.add_patch(circle)

                
                if cluster=='0717':
                    ra_center=7*15+17*(1.0/4)+32.63*(1.0/240)
                    dec_center=37+44*(1.0/60)+59.7*(1.0/3600)
                    width1=(J0717_ra_cut1[1]-J0717_ra_cut1[0])*np.cos(np.deg2rad(dec_center))
                    height1=J0717_dec_cut1[1]-J0717_dec_cut1[0]
                    width2=(J0717_ra_cut2[1]-J0717_ra_cut2[0])*np.cos(np.deg2rad(dec_center))
                    height2=J0717_dec_cut2[1]-J0717_dec_cut2[0]
                    xy1=((J0717_ra_cut1[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J0717_dec_cut1[0]-dec_center)
                    xy2=((J0717_ra_cut2[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J0717_dec_cut2[0]-dec_center)
                    rect1=patches.Rectangle(xy1,width1,height1,fill=False,edgecolor='b')
                    ax.add_patch(rect1)
                    rect2=patches.Rectangle(xy2,width2,height2,fill=False,edgecolor='b')
                    ax.add_patch(rect2)

        
                if cluster=='1149':
                    ra_center=11*15+49*(1.0/4)+35.69*(1.0/240)
                    dec_center=22+23*(1.0/60)+54.6*(1.0/3600)
                    width1=(J1149_ra_cut1[1]-J1149_ra_cut1[0])*np.cos(np.deg2rad(dec_center))
                    height1=J1149_dec_cut1[1]-J1149_dec_cut1[0]
                    width2=(J1149_ra_cut2[1]-J1149_ra_cut2[0])*np.cos(np.deg2rad(dec_center))
                    height2=J1149_dec_cut2[1]-J1149_dec_cut2[0]
                    xy1=((J1149_ra_cut1[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J1149_dec_cut1[0]-dec_center)
                    xy2=((J1149_ra_cut2[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J1149_dec_cut2[0]-dec_center)
                    rect1=patches.Rectangle(xy1,width1,height1,fill=False,edgecolor='b')
                    ax.add_patch(rect1)
                    rect2=patches.Rectangle(xy2,width2,height2,fill=False,edgecolor='b')
                    ax.add_patch(rect2)
                    
                    
                    
        
        for j in range(point_in_each_ring*(ring_quantity-1),entries):
            result_DC[count][ring_quantity-1]+=data_DC[2][j]
            result3_DC[count][ring_quantity-1]+=data_DC[3][j]
            #print j
        result_DC[count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))
        result3_DC[count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))


        ###
        for time in range(shuffle_time):
            for i in range(ring_quantity-1):
                for j in range(point_in_each_ring):
                    index=i*point_in_each_ring+j
                    result2_DC[time][count][i]+=data_DC[7+time][index]
                result2_DC[time][count][i]/=point_in_each_ring
                    
            for j in range(point_in_each_ring*(ring_quantity-1),entries):
                result2_DC[time][count][ring_quantity-1]+=data_DC[7+time][j]
            result2_DC[time][count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))
    
    
        ###calculate mean and stdev
    
        error_mean_DC[count]=np.average(result2_DC[:,count,:],0)

    
        error_stdev_DC[count]=np.std(result2_DC[:,count,:],0)
        
        ###################IM
        xticks=np.array([])
        if show_scatter==True:
            plt.figure('J'+cluster+' tangential shear scatter with '+str(ring_quantity)+' rings(IM)')
            plt.title('J'+cluster+'$\gamma $ with'+str(ring_quantity)+' rings(IM)')
            plt.gca().invert_xaxis()
            plt.scatter(data_IM[0],data_IM[1],c=data_IM[2],cmap='rainbow')
            plt.plot(0,0,marker='+', markersize=1000, color="red")
            plt.colorbar()
            plt.xlabel('$\Delta x^\circ$')
            plt.ylabel('$\Delta y^\circ$')
        
        point_in_each_ring=entries//ring_quantity
        for i in range(ring_quantity-1):
            for j in range(point_in_each_ring):
                index=i*point_in_each_ring+j
                result_IM[count][i]+=data_IM[2][index]
                result3_IM[count][i]+=data_IM[3][index]
                #print index
            result_IM[count][i]/=point_in_each_ring
            result3_IM[count][i]/=point_in_each_ring
            radius=(np.sqrt(data_IM[-1][index])+np.sqrt(data_IM[-1][index+1]))/2
            xticks=np.append(xticks,radius)
            if show_scatter==True:
                circle=plt.Circle((0,0),radius,color='r',fill=False)
                ax=plt.gca()
                ax.add_patch(circle)
                
                if cluster=='0717':
                    ra_center=7*15+17*(1.0/4)+32.63*(1.0/240)
                    dec_center=37+44*(1.0/60)+59.7*(1.0/3600)
                    width1=(J0717_ra_cut1[1]-J0717_ra_cut1[0])*np.cos(np.deg2rad(dec_center))
                    height1=J0717_dec_cut1[1]-J0717_dec_cut1[0]
                    width2=(J0717_ra_cut2[1]-J0717_ra_cut2[0])*np.cos(np.deg2rad(dec_center))
                    height2=J0717_dec_cut2[1]-J0717_dec_cut2[0]
                    xy1=((J0717_ra_cut1[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J0717_dec_cut1[0]-dec_center)
                    xy2=((J0717_ra_cut2[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J0717_dec_cut2[0]-dec_center)
                    rect1=patches.Rectangle(xy1,width1,height1,fill=False,edgecolor='b')
                    ax.add_patch(rect1)
                    rect2=patches.Rectangle(xy2,width2,height2,fill=False,edgecolor='b')
                    ax.add_patch(rect2)

        
                if cluster=='1149':
                    ra_center=11*15+49*(1.0/4)+35.69*(1.0/240)
                    dec_center=22+23*(1.0/60)+54.6*(1.0/3600)
                    width1=(J1149_ra_cut1[1]-J1149_ra_cut1[0])*np.cos(np.deg2rad(dec_center))
                    height1=J1149_dec_cut1[1]-J1149_dec_cut1[0]
                    width2=(J1149_ra_cut2[1]-J1149_ra_cut2[0])*np.cos(np.deg2rad(dec_center))
                    height2=J1149_dec_cut2[1]-J1149_dec_cut2[0]
                    xy1=((J1149_ra_cut1[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J1149_dec_cut1[0]-dec_center)
                    xy2=((J1149_ra_cut2[0]-ra_center)*np.cos(np.deg2rad(dec_center)),J1149_dec_cut2[0]-dec_center)
                    rect1=patches.Rectangle(xy1,width1,height1,fill=False,edgecolor='b')
                    ax.add_patch(rect1)
                    rect2=patches.Rectangle(xy2,width2,height2,fill=False,edgecolor='b')
                    ax.add_patch(rect2)
        
        for j in range(point_in_each_ring*(ring_quantity-1),entries):
            result_IM[count][ring_quantity-1]+=data_IM[2][j]
            result3_IM[count][ring_quantity-1]+=data_IM[3][j]
            #print j
        result_IM[count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))
        result3_IM[count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))


        ###
        for time in range(shuffle_time):
            for i in range(ring_quantity-1):
                for j in range(point_in_each_ring):
                    index=i*point_in_each_ring+j
                    result2_IM[time][count][i]+=data_IM[7+time][index]
                result2_IM[time][count][i]/=point_in_each_ring
                    
            for j in range(point_in_each_ring*(ring_quantity-1),entries):
                result2_IM[time][count][ring_quantity-1]+=data_IM[7+time][j]
            result2_IM[time][count][ring_quantity-1]/=len(range(point_in_each_ring*(ring_quantity-1),entries))
    
    
        ###calculate mean and stdev
    
        error_mean_IM[count]=np.average(result2_IM[:,count,:],0)

    
        error_stdev_IM[count]=np.std(result2_IM[:,count,:],0)
        
        




        ########back to heaven, or it's stil hell?
        xticks=(xticks+np.insert(xticks,0,0)[0:ring_quantity-1])/2
        xticks=np.append(xticks,xticks[-1]+(xticks[-1]-xticks[-2]))##add last tick
        xticks=xticks*60 ##deg to arcmin
        
        if seperate_DC_IM==False:
            
             plt.figure('J'+cluster+'averaged tangential shear for '+str(ring_quantity)+' rings')
             plt.xlabel('$r$  (arcmin)')
             plt.ylabel('$\overline{\gamma}$')
             #plt.ylim(-0.20,0.20)
             #plt.scatter(np.array(range(ring_quantity))+1,result[count][:ring_quantity], color="red",marker='s',label="actual data")
             plt.scatter(xticks,result3_DC[count][:ring_quantity],s=100, color="r",marker='x',label="DC cross shear")
             plt.scatter(xticks,result3_IM[count][:ring_quantity],s=100, color="b",marker='x',label="IM cross shear")
             plt.errorbar(xticks,result_DC[count][:ring_quantity],error_stdev_DC[count][:ring_quantity],ecolor=error_bar_color_DC,fmt=data_fmt_DC,label="DC actual data")
             plt.errorbar(xticks,result_IM[count][:ring_quantity],error_stdev_IM[count][:ring_quantity],ecolor=error_bar_color_IM,fmt=data_fmt_IM,label="IM actual data")
             
             plt.scatter(xticks,error_mean_DC[count][:ring_quantity], color=shuffled_data_color_DC,marker='s',label="DC shuffled data")
             plt.scatter(xticks,error_mean_IM[count][:ring_quantity], color=shuffled_data_color_IM,marker='s',label="IM shuffled data")
             
             plt.xticks(xticks)
             plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings')
             #plt.legend(loc='center left', bbox_to_anchor=(0.8, 1))
             plt.legend(loc=0)
             
        
        if seperate_DC_IM==True:
            
             plt.figure('J'+cluster+'averaged tangential shear for '+str(ring_quantity)+' rings(DC)')
             plt.xlabel('$r$  (arcmin)')
             plt.ylabel('$\overline{\gamma}$')
             plt.errorbar(xticks,result_DC[count][:ring_quantity],error_stdev_DC[count][:ring_quantity],ecolor=error_bar_color_DC,fmt=data_fmt_DC,label="DC actual data")
             plt.scatter(xticks,error_mean_DC[count][:ring_quantity], color=shuffled_data_color_DC,marker='s',label="DC shuffled data")
             plt.scatter(xticks,result3_DC[count][:ring_quantity],s=100,color="r",marker='x',label="DC cross shear")

             
             plt.xticks(xticks)
             plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings(DC)')
             #plt.legend(loc='center left', bbox_to_anchor=(0.8, 1))
             plt.legend(loc=0)            
            
             plt.figure('J'+cluster+'averaged tangential shear for '+str(ring_quantity)+' rings(IM)')
             plt.xlabel('$r$  (arcmin)')
             plt.ylabel('$\overline{\gamma}$')
             plt.errorbar(xticks,result_IM[count][:ring_quantity],error_stdev_IM[count][:ring_quantity],ecolor=error_bar_color_IM,fmt=data_fmt_IM,label="IM actual data")
             plt.scatter(xticks,error_mean_IM[count][:ring_quantity], color=shuffled_data_color_IM,marker='s',label="IM shuffled data")
             plt.scatter(xticks,result3_IM[count][:ring_quantity],s=100,color="b",marker='x',label="IM cross shear")
             
             plt.xticks(xticks)
             plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings(IM)')
             #plt.legend(loc='center left', bbox_to_anchor=(0.8, 1))
             plt.legend(loc=0)
        
        print 'J',cluster
        print 'DC averged tangential shear:\n',result_DC
        print 'DC mean of randomized tangential shear:\n',error_mean_DC
        print 'DC standard deviation(1/2length of whole error bar):\n',error_stdev_DC
        print 'DC ratio of std: ',error_stdev_DC/np.min(error_stdev_DC)
             
        print 'IM averged tangential shear:\n',result_IM
        print 'IM mean of randomized tangential shear:\n',error_mean_IM
        print 'IM standard deviation(1/2length of whole error bar):\n',error_stdev_IM
        print 'IM ratio of std: ',error_stdev_IM/np.min(error_stdev_IM)
        
    plt.show()
    return





def seperate_quardants_plot(cluster,data_DC,data_IM,ring,shuffle_time):
    
    data_fmt_DC='ro'  ##red circles
    error_bar_color_DC='g' ##green
    shuffled_data_color_DC='c' ##cyan

    data_fmt_IM='yo' ##yellow circles
    error_bar_color_IM='b' ##blue
    shuffled_data_color_IM='m'##magenta    
    
    entries=data_DC[0].size
    
        
    
    newrow=np.zeros((1,entries))
    data_DC=np.vstack([data_DC,newrow])
    data_IM=np.vstack([data_IM,newrow])
    ###old data(ra[0] dec[1] tangential_shear[2] randomized_shear[3:-2] radius^2[-1],[i])
    ###new data=tangential_shear(x[0] y[1] tangential_shear[2] randomized_shear[3:-2] radius^2[-1],[i])
    data_DC[-1]=data_DC[0]**2+data_DC[1]**2
    data_IM[-1]=data_IM[0]**2+data_IM[1]**2
    ##sort according to data[-1] (radius^2)
    data_DC=data_DC[:,np.argsort(data_DC[-1])]
    data_IM=data_IM[:,np.argsort(data_IM[-1])]


    result_DC_1=np.zeros((len(ring),max(ring)))
    result2_DC_1=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_DC_1=np.zeros((len(ring),max(ring)))
    error_stdev_DC_1=np.zeros((len(ring),max(ring)))
    
    result_DC_2=np.zeros((len(ring),max(ring)))
    result2_DC_2=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_DC_2=np.zeros((len(ring),max(ring)))
    error_stdev_DC_2=np.zeros((len(ring),max(ring)))
    
    result_DC_3=np.zeros((len(ring),max(ring)))
    result2_DC_3=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_DC_3=np.zeros((len(ring),max(ring)))
    error_stdev_DC_3=np.zeros((len(ring),max(ring)))
    
    result_DC_4=np.zeros((len(ring),max(ring)))
    result2_DC_4=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_DC_4=np.zeros((len(ring),max(ring)))
    error_stdev_DC_4=np.zeros((len(ring),max(ring)))
    
    
    
    
    result_IM_1=np.zeros((len(ring),max(ring)))
    result2_IM_1=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_IM_1=np.zeros((len(ring),max(ring)))
    error_stdev_IM_1=np.zeros((len(ring),max(ring)))
    
    result_IM_2=np.zeros((len(ring),max(ring)))
    result2_IM_2=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_IM_2=np.zeros((len(ring),max(ring)))
    error_stdev_IM_2=np.zeros((len(ring),max(ring)))
    
    result_IM_3=np.zeros((len(ring),max(ring)))
    result2_IM_3=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_IM_3=np.zeros((len(ring),max(ring)))
    error_stdev_IM_3=np.zeros((len(ring),max(ring)))
    
    result_IM_4=np.zeros((len(ring),max(ring)))
    result2_IM_4=np.zeros((shuffle_time,len(ring),max(ring)))
    error_mean_IM_4=np.zeros((len(ring),max(ring)))
    error_stdev_IM_4=np.zeros((len(ring),max(ring)))
    
    
    
    
    
    count=-1
    for ring_quantity in ring:
        
        count+=1
        
        ###################DC
        xticks=np.array([])

   
        point_in_each_ring=entries//ring_quantity
        for i in range(ring_quantity-1):
            first_quardant_count=0
            second_quardant_count=0
            third_quardant_count=0
            fourth_quardant_count=0
            for j in range(point_in_each_ring):
                ###new data=tangential_shear(x[0] y[1] tangential_shear[2] randomized_shear[3:-2] radius^2[-1],[i])
                index=i*point_in_each_ring+j
                
                if data_DC[0][index]>0 and data_DC[1][index]>=0:
                    result_DC_1[count][i]+=data_DC[2][index]
                    first_quardant_count+=1
                if data_DC[0][index]<=0 and data_DC[1][index]>0:
                    result_DC_2[count][i]+=data_DC[2][index]
                    second_quardant_count+=1
                if data_DC[0][index]<0 and data_DC[1][index]<=0:
                    result_DC_3[count][i]+=data_DC[2][index]
                    third_quardant_count+=1
                if data_DC[0][index]>=0 and data_DC[1][index]<0:
                    result_DC_4[count][i]+=data_DC[2][index]
                    fourth_quardant_count+=1
                
                
            result_DC_1[count][i]/=first_quardant_count
            result_DC_2[count][i]/=second_quardant_count
            result_DC_3[count][i]/=third_quardant_count
            result_DC_4[count][i]/=fourth_quardant_count
            
            #if first_quardant_count+second_quardant_count+third_quardant_count+fourth_quardant_count==point_in_each_ring:
                #print 'sources in 1-4 quardants in each bin:',first_quardant_count,second_quardant_count,third_quardant_count,fourth_quardant_count
            
            radius=(np.sqrt(data_DC[-1][index])+np.sqrt(data_DC[-1][index+1]))/2
            xticks=np.append(xticks,radius)



        first_quardant_count=0
        second_quardant_count=0
        third_quardant_count=0
        fourth_quardant_count=0
        for j in range(point_in_each_ring*(ring_quantity-1),entries):
            if data_DC[0][j]>0 and data_DC[1][j]>=0:
                result_DC_1[count][ring_quantity-1]+=data_DC[2][j]
                first_quardant_count+=1
            if data_DC[0][j]<=0 and data_DC[1][j]>0:
                result_DC_2[count][ring_quantity-1]+=data_DC[2][j]
                second_quardant_count+=1
            if data_DC[0][j]<0 and data_DC[1][j]<=0:
                result_DC_3[count][ring_quantity-1]+=data_DC[2][j]
                third_quardant_count+=1
            if data_DC[0][j]>=0 and data_DC[1][j]<0:
                result_DC_4[count][ring_quantity-1]+=data_DC[2][j]
                fourth_quardant_count+=1


        result_DC_1[count][ring_quantity-1]/=first_quardant_count
        result_DC_2[count][ring_quantity-1]/=second_quardant_count
        result_DC_3[count][ring_quantity-1]/=third_quardant_count
        result_DC_4[count][ring_quantity-1]/=fourth_quardant_count
        
        #if first_quardant_count+second_quardant_count+third_quardant_count+fourth_quardant_count==len(range(point_in_each_ring*(ring_quantity-1),entries)):
            #print 'sources in 1-4 quardants in each bin:',first_quardant_count,second_quardant_count,third_quardant_count,fourth_quardant_count


        ###
        for time in range(shuffle_time):
            for i in range(ring_quantity-1):
                first_quardant_count=0
                second_quardant_count=0
                third_quardant_count=0
                fourth_quardant_count=0
                for j in range(point_in_each_ring):
                    index=i*point_in_each_ring+j
                    if data_DC[0][index]>0 and data_DC[1][index]>=0:
                        result2_DC_1[time][count][i]+=data_DC[7+time][index]
                        first_quardant_count+=1
                    if data_DC[0][index]<=0 and data_DC[1][index]>0:
                        result2_DC_2[time][count][i]+=data_DC[7+time][index]
                        second_quardant_count+=1
                    if data_DC[0][index]<0 and data_DC[1][index]<=0:
                        result2_DC_3[time][count][i]+=data_DC[7+time][index]
                        third_quardant_count+=1
                    if data_DC[0][index]>=0 and data_DC[1][index]<0:
                        result2_DC_4[time][count][i]+=data_DC[7+time][index]
                        fourth_quardant_count+=1
                result2_DC_1[time][count][i]/=first_quardant_count
                result2_DC_2[time][count][i]/=second_quardant_count
                result2_DC_3[time][count][i]/=third_quardant_count
                result2_DC_4[time][count][i]/=fourth_quardant_count
                    
                    

            first_quardant_count=0
            second_quardant_count=0
            third_quardant_count=0
            fourth_quardant_count=0        
            for j in range(point_in_each_ring*(ring_quantity-1),entries):
                if data_DC[0][j]>0 and data_DC[1][j]>=0:
                    result2_DC_1[time][count][ring_quantity-1]+=data_DC[7+time][j]
                    first_quardant_count+=1
                if data_DC[0][j]<=0 and data_DC[1][j]>0:
                    result2_DC_2[time][count][ring_quantity-1]+=data_DC[7+time][j]
                    second_quardant_count+=1
                if data_DC[0][j]<0 and data_DC[1][j]<=0:
                    result2_DC_3[time][count][ring_quantity-1]+=data_DC[7+time][j]
                    third_quardant_count+=1
                if data_DC[0][j]>=0 and data_DC[1][j]<0:
                    result2_DC_4[time][count][ring_quantity-1]+=data_DC[7+time][j]
                    fourth_quardant_count+=1
                
            result2_DC_1[time][count][ring_quantity-1]/=first_quardant_count
            result2_DC_2[time][count][ring_quantity-1]/=second_quardant_count
            result2_DC_3[time][count][ring_quantity-1]/=third_quardant_count
            result2_DC_4[time][count][ring_quantity-1]/=fourth_quardant_count 
                

    
    
        ###calculate mean and stdev
    
        error_mean_DC_1[count]=np.average(result2_DC_1[:,count,:],0)
        error_mean_DC_2[count]=np.average(result2_DC_2[:,count,:],0)
        error_mean_DC_3[count]=np.average(result2_DC_3[:,count,:],0)
        error_mean_DC_4[count]=np.average(result2_DC_4[:,count,:],0)

    
        error_stdev_DC_1[count]=np.std(result2_DC_1[:,count,:],0)
        error_stdev_DC_2[count]=np.std(result2_DC_2[:,count,:],0)
        error_stdev_DC_3[count]=np.std(result2_DC_3[:,count,:],0)
        error_stdev_DC_4[count]=np.std(result2_DC_4[:,count,:],0)
        
        
        ###################IM
        xticks=np.array([])

   
        point_in_each_ring=entries//ring_quantity
        for i in range(ring_quantity-1):
            first_quardant_count=0
            second_quardant_count=0
            third_quardant_count=0
            fourth_quardant_count=0
            for j in range(point_in_each_ring):
                ###new data=tangential_shear(x[0] y[1] tangential_shear[2] randomized_shear[3:-2] radius^2[-1],[i])
                index=i*point_in_each_ring+j
                
                if data_IM[0][index]>0 and data_IM[1][index]>=0:
                    result_IM_1[count][i]+=data_IM[2][index]
                    first_quardant_count+=1
                if data_IM[0][index]<=0 and data_IM[1][index]>0:
                    result_IM_2[count][i]+=data_IM[2][index]
                    second_quardant_count+=1
                if data_IM[0][index]<0 and data_IM[1][index]<=0:
                    result_IM_3[count][i]+=data_IM[2][index]
                    third_quardant_count+=1
                if data_IM[0][index]>=0 and data_IM[1][index]<0:
                    result_IM_4[count][i]+=data_IM[2][index]
                    fourth_quardant_count+=1
                
                
            result_IM_1[count][i]/=first_quardant_count
            result_IM_2[count][i]/=second_quardant_count
            result_IM_3[count][i]/=third_quardant_count
            result_IM_4[count][i]/=fourth_quardant_count
            

            
            radius=(np.sqrt(data_IM[-1][index])+np.sqrt(data_IM[-1][index+1]))/2
            xticks=np.append(xticks,radius)



        first_quardant_count=0
        second_quardant_count=0
        third_quardant_count=0
        fourth_quardant_count=0
        for j in range(point_in_each_ring*(ring_quantity-1),entries):
            if data_IM[0][j]>0 and data_IM[1][j]>=0:
                result_IM_1[count][ring_quantity-1]+=data_IM[2][j]
                first_quardant_count+=1
            if data_IM[0][j]<=0 and data_IM[1][j]>0:
                result_IM_2[count][ring_quantity-1]+=data_IM[2][j]
                second_quardant_count+=1
            if data_IM[0][j]<0 and data_IM[1][j]<=0:
                result_IM_3[count][ring_quantity-1]+=data_IM[2][j]
                third_quardant_count+=1
            if data_IM[0][j]>=0 and data_IM[1][j]<0:
                result_IM_4[count][ring_quantity-1]+=data_IM[2][j]
                fourth_quardant_count+=1


        result_IM_1[count][ring_quantity-1]/=first_quardant_count
        result_IM_2[count][ring_quantity-1]/=second_quardant_count
        result_IM_3[count][ring_quantity-1]/=third_quardant_count
        result_IM_4[count][ring_quantity-1]/=fourth_quardant_count
        



        ###
        for time in range(shuffle_time):
            for i in range(ring_quantity-1):
                first_quardant_count=0
                second_quardant_count=0
                third_quardant_count=0
                fourth_quardant_count=0
                for j in range(point_in_each_ring):
                    index=i*point_in_each_ring+j
                    if data_IM[0][index]>0 and data_IM[1][index]>=0:
                        result2_IM_1[time][count][i]+=data_IM[7+time][index]
                        first_quardant_count+=1
                    if data_IM[0][index]<=0 and data_IM[1][index]>0:
                        result2_IM_2[time][count][i]+=data_IM[7+time][index]
                        second_quardant_count+=1
                    if data_IM[0][index]<0 and data_IM[1][index]<=0:
                        result2_IM_3[time][count][i]+=data_IM[7+time][index]
                        third_quardant_count+=1
                    if data_IM[0][index]>=0 and data_IM[1][index]<0:
                        result2_IM_4[time][count][i]+=data_IM[7+time][index]
                        fourth_quardant_count+=1
                result2_IM_1[time][count][i]/=first_quardant_count
                result2_IM_2[time][count][i]/=second_quardant_count
                result2_IM_3[time][count][i]/=third_quardant_count
                result2_IM_4[time][count][i]/=fourth_quardant_count
                    
                    

            first_quardant_count=0
            second_quardant_count=0
            third_quardant_count=0
            fourth_quardant_count=0        
            for j in range(point_in_each_ring*(ring_quantity-1),entries):
                if data_IM[0][j]>0 and data_IM[1][j]>=0:
                    result2_IM_1[time][count][ring_quantity-1]+=data_IM[7+time][j]
                    first_quardant_count+=1
                if data_IM[0][j]<=0 and data_IM[1][j]>0:
                    result2_IM_2[time][count][ring_quantity-1]+=data_IM[7+time][j]
                    second_quardant_count+=1
                if data_IM[0][j]<0 and data_IM[1][j]<=0:
                    result2_IM_3[time][count][ring_quantity-1]+=data_IM[7+time][j]
                    third_quardant_count+=1
                if data_IM[0][j]>=0 and data_IM[1][j]<0:
                    result2_IM_4[time][count][ring_quantity-1]+=data_IM[7+time][j]
                    fourth_quardant_count+=1
                
            result2_IM_1[time][count][ring_quantity-1]/=first_quardant_count
            result2_IM_2[time][count][ring_quantity-1]/=second_quardant_count
            result2_IM_3[time][count][ring_quantity-1]/=third_quardant_count
            result2_IM_4[time][count][ring_quantity-1]/=fourth_quardant_count
            
            
                ###calculate mean and stdev
    
        error_mean_IM_1[count]=np.average(result2_IM_1[:,count,:],0)
        error_mean_IM_2[count]=np.average(result2_IM_2[:,count,:],0)
        error_mean_IM_3[count]=np.average(result2_IM_3[:,count,:],0)
        error_mean_IM_4[count]=np.average(result2_IM_4[:,count,:],0)

    
        error_stdev_IM_1[count]=np.std(result2_IM_1[:,count,:],0)
        error_stdev_IM_2[count]=np.std(result2_IM_2[:,count,:],0)
        error_stdev_IM_3[count]=np.std(result2_IM_3[:,count,:],0)
        error_stdev_IM_4[count]=np.std(result2_IM_4[:,count,:],0)



        ########back to heaven, or it's stil hell?
        xticks=(xticks+np.insert(xticks,0,0)[0:ring_quantity-1])/2
        xticks=np.append(xticks,xticks[-1]+(xticks[-1]-xticks[-2]))##add last tick
        xticks=xticks*60 ##deg to arcmin
        

        '''    
        plt.figure('J'+cluster+'averaged tangential shear for '+str(ring_quantity)+' rings for four quardants (subplot position matches scatter plot position)')
        plt.subplot(2,2,1)
        plt.xlabel('$r$  (arcmin)')
        plt.ylabel('$\overline{\gamma}$')
        plt.ylim((-0.25,0.4))
        plt.errorbar(xticks,result_DC_1[count][:ring_quantity],error_stdev_DC_1[count][:ring_quantity],ecolor=error_bar_color_DC,fmt=data_fmt_DC,label="DC actual data")
        plt.errorbar(xticks,result_IM_1[count][:ring_quantity],error_stdev_IM_1[count][:ring_quantity],ecolor=error_bar_color_IM,fmt=data_fmt_IM,label="IM actual data")             
        plt.scatter(xticks,error_mean_DC_1[count][:ring_quantity], color=shuffled_data_color_DC,marker='s',label="DC shuffled data")
        plt.scatter(xticks,error_mean_IM_1[count][:ring_quantity], color=shuffled_data_color_IM,marker='s',label="IM shuffled data")
        plt.xticks(xticks)
        plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings 1st quardant')
        plt.legend(loc=0)
        
        plt.subplot(2,2,2)
        plt.xlabel('$r$  (arcmin)')
        plt.ylabel('$\overline{\gamma}$')
        plt.ylim((-0.25,0.4))
        plt.errorbar(xticks,result_DC_2[count][:ring_quantity],error_stdev_DC_2[count][:ring_quantity],ecolor=error_bar_color_DC,fmt=data_fmt_DC,label="DC actual data")
        plt.errorbar(xticks,result_IM_2[count][:ring_quantity],error_stdev_IM_2[count][:ring_quantity],ecolor=error_bar_color_IM,fmt=data_fmt_IM,label="IM actual data")             
        plt.scatter(xticks,error_mean_DC_2[count][:ring_quantity], color=shuffled_data_color_DC,marker='s',label="DC shuffled data")
        plt.scatter(xticks,error_mean_IM_2[count][:ring_quantity], color=shuffled_data_color_IM,marker='s',label="IM shuffled data")
        plt.xticks(xticks)
        plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings 2nd quardant')
        plt.legend(loc=0)
        
        plt.subplot(2,2,4)
        plt.xlabel('$r$  (arcmin)')
        plt.ylabel('$\overline{\gamma}$')
        plt.ylim((-0.25,0.4))
        plt.errorbar(xticks,result_DC_3[count][:ring_quantity],error_stdev_DC_3[count][:ring_quantity],ecolor=error_bar_color_DC,fmt=data_fmt_DC,label="DC actual data")
        plt.errorbar(xticks,result_IM_3[count][:ring_quantity],error_stdev_IM_3[count][:ring_quantity],ecolor=error_bar_color_IM,fmt=data_fmt_IM,label="IM actual data")             
        plt.scatter(xticks,error_mean_DC_3[count][:ring_quantity], color=shuffled_data_color_DC,marker='s',label="DC shuffled data")
        plt.scatter(xticks,error_mean_IM_3[count][:ring_quantity], color=shuffled_data_color_IM,marker='s',label="IM shuffled data")
        plt.xticks(xticks)
        plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings 3rd quardant')
        plt.legend(loc=0)
        
        plt.subplot(2,2,3)
        plt.xlabel('$r$  (arcmin)')
        plt.ylabel('$\overline{\gamma}$')
        plt.ylim((-0.25,0.4))
        plt.errorbar(xticks,result_DC_4[count][:ring_quantity],error_stdev_DC_4[count][:ring_quantity],ecolor=error_bar_color_DC,fmt=data_fmt_DC,label="DC actual data")
        plt.errorbar(xticks,result_IM_4[count][:ring_quantity],error_stdev_IM_4[count][:ring_quantity],ecolor=error_bar_color_IM,fmt=data_fmt_IM,label="IM actual data")             
        plt.scatter(xticks,error_mean_DC_4[count][:ring_quantity], color=shuffled_data_color_DC,marker='s',label="DC shuffled data")
        plt.scatter(xticks,error_mean_IM_4[count][:ring_quantity], color=shuffled_data_color_IM,marker='s',label="IM shuffled data")
        plt.xticks(xticks)
        plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings 4th quardant')
        plt.legend(loc=0)
        '''
        
        plt.figure('J'+cluster+'averaged tangential shear for '+str(ring_quantity)+' rings for four quardants')
        plt.xlabel('$r$  (arcmin)')
        plt.ylabel('$\overline{\gamma}$')
        
        plt.errorbar(xticks,result_IM_1[count][:ring_quantity],error_stdev_IM_1[count][:ring_quantity],ecolor='r',fmt='ro',label="1st quadrant")             
        plt.errorbar(xticks,result_IM_2[count][:ring_quantity],error_stdev_IM_2[count][:ring_quantity],ecolor='g',fmt='go',label="2nd quadrant")             
        plt.errorbar(xticks,result_IM_3[count][:ring_quantity],error_stdev_IM_3[count][:ring_quantity],ecolor='b',fmt='bo',label="3rd quadrant")             
        plt.errorbar(xticks,result_IM_4[count][:ring_quantity],error_stdev_IM_4[count][:ring_quantity],ecolor='c',fmt='co',label="4th quadrant")             

        plt.xticks(xticks)
        plt.title('J'+cluster+'   $\overline{\gamma}$   for'+str(ring_quantity)+' rings 1-4 quardrants')
        plt.legend(loc=0)        
        
        
        
    plt.show()
    return  result_DC_1,result_DC_2,result_DC_3,result_DC_4  
    
    
def area_tangential_shear(catalog,tangential_shear,cluster):
        ##catalog[i,(name_0,ra_1,dec_2,im_maj_3,im_min_4,im_pa_5,dc_maj_6,dc_min_7,dc_pa_8)]
    ###tangential_shear[(0-x,1-y,2-tangential_shear,3~end-shuffled_tangential_shear),i]
    entries=catalog.size
    x=np.zeros(entries)
    y=np.zeros(entries)
    for i in range(entries):
        y[i]=catalog[i][3]**2+catalog[i][4]**2
        x[i]=tangential_shear[2][i]
        
    plt.figure('tangential_shear_vs_a^2+b^2'+'J'+cluster)
    plt.scatter(x,y)
    plt.title('tangential_shear_vs_a^2+b^2')
    plt.xlabel('tangential_shear')
    plt.ylabel('a^2+b^2')
    plt.show()
    return
        
