#!/usr/local/bin/python
#-*- coding: utf-8 -*-


from __future__ import division

__all__ = ["Plot composite spectrum"]
__version__ = "1.0.0"
__author__ = "Mikkel Stockmann (mikkelstockmann@gmail.com)"
__copyright__ = "Copyright 2018 Mikkel Stockmann"


import numpy as np
import glob
import matplotlib.pylab as pl


import pyfits as pf
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit
import sys
import methods as st





####################################################

Path_model_gallazzi = glob.glob('data/model/QG_stack_cc_bin9A.full_spec.SanV4p1z2.bestfit.fits')
f = pf.open(Path_model_gallazzi[0])
W_model = f[2].data['LAMBDA']
F_model = f[2].data['FLUX']


path_composite=glob.glob('output/CcQG_v1.txt')[0]
data = np.genfromtxt(path_composite)
W, F, E, M, N = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]



plt.figure(figsize=(14,4))
plt.subplots_adjust(bottom=0.15,right=0.98,left=0.05,top=0.95)
plt.plot(W,F,color='black')
plt.plot(W,M,color='red')
plt.fill_between(W,F-E,F+E,color='blue',alpha=0.3)
plt.plot(W_model,F_model,color='red')
st.plot_abs_em_lines(3300,6000)


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

OII = 3726.2
lines = np.array([(Na_5896+Na_5890)/2, (Mg_5169+Mg_5175+Mg_5183)/3, Fe_5332, Fe_5401, CaI_4227, H_beta, Fe_4384, H_gamma, G_4304, H_delta, B4000, H_eps, H_3969, K, H8, H9, H10_3798,H11_3770,H12_3749,H13_3733,OII])
name_list = [r'$Na$',r'$MgI$', r'$Fe[5332]$', r'$Fe[5401]$', r'$CaI$', r'$H\beta$', r'$Fe[4384]$', r'$H\gamma$', r'$G[4304]$', r'$H\delta$', r'$H\epsilon$', r'$H[3969]$', r'$K$', r'$H8$', r'$H9$', r'$H10$',r'$H11$',r'$H12$',r'$H13$',r'$[OII]$']
position = [0.4,0.4,0.5,0.4,0.4,1.3,0.4,1.3,0.5,1.3,1.1,0.2,0.3,0.9,1.0,1.1,1.2,1.3,1.4,0.7]
x_corr = [-60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-90]

color_list = ['black','black','black','black','black','red','black','red','black','red','red','black','black','red','red','red','red','red','red','blue']
line_widths = [1,1,1,1,1,10,1,10,1,10,10,1,1,10,10,10,10,10,10,10,10]


for j in range(len(position)):
	if line_widths[j] == 1:	
		plt.vlines(lines[j],0,5,color=color_list[j],alpha=0.5,linewidth=0.5)
		plt.text(lines[j]+x_corr[j],position[j],name_list[j],rotation=0,fontsize=15)
	else:
		plt.vlines(lines[j],0,5,color=color_list[j],alpha=0.5,linewidth=1)
		plt.text(lines[j]+x_corr[j],position[j],name_list[j],rotation=0,fontsize=15)

plt.ylabel('Normalised Flux',fontsize=15)
plt.xlabel('$\lambda[\AA{}]$',fontsize=15)
plt.axis([3300,6000,0.2,1.5])

outname = 'output/figure.pdf'
plt.savefig(outname)
os.system('pdfcrop %s' % outname)
os.system('rm %s' % outname)
os.system('mv %s-crop.pdf %s' % (outname.split('.pdf')[0],outname))






