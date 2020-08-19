#!/usr/local/bin/python
#-*- coding: utf-8 -*-


from __future__ import division

import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

def read_in_1d_fits(path):
    """ Read in 1d fits file
        F: flux, E: error, W: wavelength (Angstrom), hd: header
        Returns wavelength, flux, error, header
    """
    data_arr = pf.open(path)
    hdf = data_arr[0].header
    hde = data_arr[1].header
    hdm = data_arr[2].header

    F = data_arr[0].data
    E = data_arr[1].data
    M = data_arr[2].data

    # Wavelength in nanometers or Angstrom
    if not hdf.comments['CRVAL1']:# nm
        W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])*10
    elif hdf.comments['CRVAL1'] == 'AA': # AA
        W = (hdf['CRVAL1'] + (hdf['CRPIX1'] - 1 + np.arange(hdf['NAXIS1']))*hdf['CDELT1'])

    return W, F, E, M, hdf, hde, hdm


def plot_abs_em_lines(w1,w2):
    """ Load the galaxylines.txt list and plot the absorption (red,2)
        and emission (blue,1) lines from lambda = [w1:w2]
    """
    data = np.genfromtxt('data/galaxylines.txt',skip_header=11,dtype=str)
    
    line_wave = data[:,0]
    line_type = data[:,1]
    line_name = data[:,2]

    # select the range
    rule1 = (line_wave > '%s' % w1) == True
    rule2 = (line_wave < '%s' % w2) == True
    
    ID = np.where(rule1*rule2)
    line_wave = line_wave[ID]
    line_type = line_type[ID]
    line_name = line_name[ID]

    
    for nn in range(len(line_wave)):
        if line_type[nn] == '2':
            plt.vlines(line_wave[nn],-10,50,color='red')
            # plt.text(line_wave[nn],0,line_name[nn],rotation=-90)
        elif line_type[nn] == '1':
            plt.vlines(line_wave[nn],-10,50,color='blue')
            # plt.text(line_wave[nn],0,line_name[nn],rotation=-90)

def load_spec_z(path,obj,error_z=False,verbal=True):
    """ Load the spectroscopic redshifts and assymetric uncertainties of a galaxy
    """
    data = np.genfromtxt(path,dtype=str)
    if error_z == False:
        Name = data[:,0]
        z_spec = data[:,1]
        loc_name = np.where(obj == Name)[0][0]
        if verbal == True:
            print 'Redshift (%s): %s loaded' % (obj,np.float(z_spec[loc_name]))
        return np.float(z_spec[loc_name])
    elif error_z == True:
        Name = data[:,0]
        z_spec = data[:,1]
        z_errl = data[:,2]
        z_errh = data[:,3]

        loc_name = np.where(obj == Name)[0][0]
        if verbal == True:
            print 'Redshift (%s): %s loaded' % (obj,np.float(z_spec[loc_name]))

        z, errl, errh = np.around(z_spec[loc_name].astype(float),5), np.around(z_errl[loc_name].astype(float),5), np.around(z_errh[loc_name].astype(float),5)

        return z, errl, errh


def convert_2_rest_frame(W,F,E,z):
    """ Rest frame wavelength conversion of galaxy spectrum information
    """
    W_rest_frame = W / (1+z)
    F_rest_frame = F * (1+z)
    E_rest_frame = E * (1+z)
    return W_rest_frame, F_rest_frame, E_rest_frame

def wavelength_range(Warr):
    """ Takes 2 dimensional array and returns the min and max wavelength
    """
    wlmin = []
    wlmax = []
    for jj in range(len(Warr)):
        wlmin.append(min(Warr[jj]))
        wlmax.append(max(Warr[jj]))
    return np.min(wlmin), np.max(wlmax)

def common_wavelength(wlarr_old, wlarr_new, fluxarr_old, intp_type, fill_value = 0.):
    """ Array interpolation
    """
    from scipy import interpolate
    f = interpolate.interp1d(wlarr_old, fluxarr_old, kind=intp_type, bounds_error = False, fill_value=fill_value)
    fluxarr_new = f(wlarr_new)
    return fluxarr_new

def mask_skylines_areas(W,M,zspec):
    """ Mask the skyline regions in rest frame spectrum
    """
    bad_gaps = [10400,13600,14200,18310,19340]
    bad_gaps_rf = bad_gaps/(1+zspec)
    M_sky = np.copy(M)

    id1 = np.where(W < bad_gaps_rf[0])
    id2 = np.where((W > bad_gaps_rf[1]) & (W < bad_gaps_rf[2]))
    id3 = np.where((W > bad_gaps_rf[3]) & (W < bad_gaps_rf[4]))
    M_sky[id1] = 1
    M_sky[id2] = 1
    M_sky[id3] = 1
    return M_sky

def total_stack_wmean(F,E,M):
    """ Stacking of array F using the weights = 1/E**2 , and excluding pixels with M=1
        We use the python masking arrays np.ma.array(), so only values with flag, M=0 is included.
        For the values where all pixel that is combine have bad flag (M=1), we combine all the bad pixels in a weighted average and flag it with M=1, contrary to all the masked combinations that get M=0 
    """

    Fwm = np.zeros(F.shape[1])
    Ewm = np.zeros(F.shape[1])
    Mwm = np.zeros(F.shape[1])
    No_of_combined_spec = np.zeros(F.shape[1])
    for i, k in enumerate(F.transpose()):

        # Number of combinations - count the number of good flags in pixel i's combination
        ID_map0 = np.where(M.transpose()[i] == 0)[0]
        No_of_combined_spec[i] = len(ID_map0)

        # create masked array
        F_ma = np.ma.array(k, mask=[M.transpose()[i]])

        # weighted average
        w_err = 1/(E.transpose()[i])**2
        tmp = np.ma.average(F_ma, axis=0, weights=w_err, returned=True)

        if tmp:
            Fwm[i] = tmp[0]
            Ewm[i] = 1/np.sqrt(tmp[1])
            Mwm[i] = 0
        else:
            Fwm[i] = np.sum(k*w_err)/np.sum(w_err)
            Ewm[i] = 1/np.sqrt(np.sum(w_err))
            Mwm[i] = 1

    return Fwm, Ewm, Mwm, No_of_combined_spec


def total_stack_wmean_ewsv(F,E,M):
    """ Similar to total_stack_wmean() but here we construct the errors based on a weighted sample variance

        sigma_weighted_sample_variance**2 = sum((x_i-mu)**2*w_i)/sum(w_i)
        here mu: is the weighted average flux (Fwm[i]), x_i is the pixels to combine (F_ma)
        w_i: 1/(E.transpose()[i])**2

    """
    Fwm = np.zeros(F.shape[1])
    Ewsv = np.zeros(F.shape[1])
    Mwm = np.zeros(F.shape[1])
    No_of_combined_spec = np.zeros(F.shape[1])
    for i, k in enumerate(F.transpose()):
    
        # Number of combinations - count the number of good flags in pixel i's combination
        ID_map0 = np.where(M.transpose()[i] == 0)[0]
        No_of_combined_spec[i] = len(ID_map0)

        # create masked array
        F_ma = np.ma.array(k, mask=[M.transpose()[i]])

        # weighted average
        w_err = 1/(E.transpose()[i])**2
        tmp = np.ma.average(F_ma, axis=0, weights=w_err, returned=True)

        if tmp:
            tmp2 = np.ma.average((F_ma-tmp[0])**2, axis=0, weights=w_err, returned=True)
            
            Fwm[i] = tmp[0]
            Ewsv[i] = np.sqrt(tmp2[0])
            Mwm[i] = 0

        else:
            Fwm[i] = np.sum(k*w_err)/np.sum(w_err)
            Ewsv[i] = 1/np.sqrt(np.sum(w_err))
            Mwm[i] = 1

    return Fwm, Ewsv, Mwm, No_of_combined_spec



def total_stack_wmean_JackKnifeError(F,E,M,W):
    """ A wrapper function for running the total_stack_wmean and also calculating the jack knife error estimate.
        The jack knife error estimate sums up the difference between the total_stack and total_stack-spec(i), which is really just the noise that spec(i) adds to the stack, following this it is weighted by (N-1)/N:

        sigma_JK**2 = (N-1)/N * sum[(f-f(i))**2]

    """

    Ftot, Etot, Mtot, N_combined_spec = total_stack_wmean(F,E,M)

    Nspec = F.shape[0]
    Npixel = F.shape[1]
    array = np.zeros(shape=(Nspec,Npixel))
    array_fi = np.zeros(shape=(Nspec,Npixel))
    for i in range(Nspec):
        Fi = np.delete(F, i, 0)
        Ei = np.delete(E, i, 0)
        Mi = np.delete(M, i, 0)

        Fs_i, Es_i, Ms_i, tmp = total_stack_wmean(Fi,Ei,Mi)
        array_fi[i,:] = Fs_i

        array[i,:] = (Ftot-Fs_i)**2

    # error calculation
    sigma_JK = np.sqrt((Nspec-1)/Nspec*np.sum(array,axis=0))
    

    return Ftot, sigma_JK, Mtot, N_combined_spec



def combine_composite(W,F,E,M,comb_method):
    """ Combine the spectra using Weighted mean ('weighted_mean'), Jack Knife ('Jack_Knife'),
        Weighted Sample Variance ('weighted_sample_variance')
    """
    if comb_method == 'Jack_Knife':
        Ffinal, Efinal, Mfinal, N_combined_spec = total_stack_wmean_JackKnifeError(F,E,M,W)
        return Ffinal, Efinal, Mfinal, W, N_combined_spec
    
    elif comb_method == 'weighted_mean':
        Ffinal, Efinal, Mfinal, N_combined_spec = total_stack_wmean(F,E,M)
        return Ffinal, Efinal, Mfinal, W, N_combined_spec

    elif comb_method == 'weighted_sample_variance':
        Ffinal, Efinal, Mfinal, N_combined_spec = total_stack_wmean_ewsv(F,E,M)
        return Ffinal, Efinal, Mfinal, W, N_combined_spec



def Make_composite(Path,out_path,intp_bin,comb_method,save=False):
    """ Combining information from composite.py alongside the functions to produce the composite spectrum
    """
    n_obj = len(Path)
    Wave_arr    = []
    Flux_arr    = []
    Errs_arr    = []
    Maps_arr    = []
    target_list = []
    SN_array    = []
    z_spec_list = np.zeros(n_obj)
    for ii in range(n_obj):
        target_list.append(Path[ii].split('/')[-1].split('_')[0]) 
        Wave, Flux, Errs, M, hdf, hde, hdm = read_in_1d_fits(Path[ii])
        z_spec = load_spec_z('data/z_spec_all_v2.txt',target_list[ii])
        W_rest, F_rest, E_rest = convert_2_rest_frame(Wave,Flux,Errs,z_spec)
        z_spec_list[ii] = z_spec

        # Normalise the spectrum to the H-band photometry (4500-4730 Angstrom) 
        ID_4500_5200 = np.searchsorted(W_rest,np.array([4500,4730]))
        F_av = np.median(F_rest[ID_4500_5200[0]:ID_4500_5200[1]])
        SN_valH = np.median(F_rest[ID_4500_5200[0]:ID_4500_5200[1]]/E_rest[ID_4500_5200[0]:ID_4500_5200[1]])
        SN_array.append(SN_valH)

        Wave_arr.append(W_rest)
        Flux_arr.append((F_rest/F_av))
        Errs_arr.append((E_rest/F_av))
        Maps_arr.append(M)

    wl_lower, wl_upper = wavelength_range(Wave_arr)
    Wave_common = np.arange(wl_lower, wl_upper+intp_bin, intp_bin)
    print 'Length of interpolated wavelength axis: %s' % len(Wave_common)

    Flux_intp = np.zeros((n_obj,len(Wave_common)))
    Errs_intp = np.zeros((n_obj,len(Wave_common)))
    BPMs_intp = np.zeros((n_obj,len(Wave_common)))
    for kk in range(n_obj):
        Flux_intp[kk] = common_wavelength(Wave_arr[kk], Wave_common, Flux_arr[kk],'linear')
        Errs2_tmp = common_wavelength(Wave_arr[kk], Wave_common, Errs_arr[kk]**2,'linear',fill_value=1.0)
        Errs_intp[kk] = np.sqrt(Errs2_tmp)
        Maps_tmp = common_wavelength(Wave_arr[kk], Wave_common, Maps_arr[kk],'nearest',fill_value=1.0)
        BPMs_intp[kk] = mask_skylines_areas(Wave_common,Maps_tmp,z_spec_list[kk])


    # Combination method
    F_final, E_final, M_final, W_final, N_combine = combine_composite(Wave_common,Flux_intp,Errs_intp,BPMs_intp,comb_method)

    if save == True:
        Arr_out = np.zeros(shape=(len(F_final),5))
        Arr_out[:,0], Arr_out[:,1], Arr_out[:,2], Arr_out[:,3], Arr_out[:,4]  = W_final, F_final, E_final, M_final, N_combine
        np.savetxt(out_path,Arr_out,fmt='%g', delimiter='\t', header='Wave[A] F[norm] E[norm] Flag[0:good] N_Combinations')



def plot_absorption_lines():
    """ Plot absorption lines with corresponding name
    """ 
    # Absorption lines
    Mg_5169    = 5168.6
    Mg_5175    = 5175.4
    Mg_5183    = 5183.3
    Fe_5332    = 5331.5
    Fe_5401    = 5401.4
    Na_5890    = 5889.95
    Na_5896    = 5895.92

    CaI_4227    = 4226.7
    H_beta      = 4861.3
    Fe_4384     = 4383.6
    H_gamma     = 4340.5
    G_4304      = 4304.4
    H_delta     = 4101.7
    B4000       = 4000
    H_eps       = 3970.0
    H_3969      = 3968.5
    K           = 3933.7
    H8          = 3889.1
    H9          = 3835.4
    CN_3833     = 3833
    H10_3798    = 3798
    H11_3770    = 3769.7
    H12_3749    = 3749.2
    H13_3733    = 3733.4
 

    ## Plot absorption lines
    lines = np.array([Na_5896,Na_5890,Mg_5169, Mg_5175, Mg_5183, Fe_5332, Fe_5401, CaI_4227, H_beta, Fe_4384, H_gamma, G_4304, H_delta, B4000, H_eps, H_3969, K, H8, H9, CN_3833, H10_3798,H11_3770,H12_3749,H13_3733])
    name_list = ['Na_5896','Na_5890','Mg_5169', 'Mg_5175', 'Mg_5183', 'Fe_5332', 'Fe_5401', 'CaI_4227', 'H_beta', 'Fe_4384', 'H_gamma', 'G_4304', 'H_delta', 'B4000', 'H_eps', 'H_3969', 'K', 'H8', 'H9', 'CN_3833', 'H10_3798','H11_3770','H12_3749','H13_3733']
    for j in range(len(lines)):
        if lines[j] == 4000:
            plt.vlines(lines[j],0,5,color='pink')
            plt.text(lines[j],0,name_list[j],rotation=-90)
        else:
            plt.vlines(lines[j],0,5,color='r')
            plt.text(lines[j],0,name_list[j],rotation=-90)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@


