#!/usr/local/bin/python
#-*- coding: utf-8 -*-

from __future__ import division

__all__ = ["Composite compact Quiescent Galaxies"]
__version__ = "1.0.0"
__author__ = "Mikkel Stockmann (mikkelstockmann@gmail.com)"
__copyright__ = "Copyright 2018 Mikkel Stockmann"

####################################################

""" Path: Path of galaxy spectra
	interp_step: Interpolation size: 3 Angstrom (9 Angstrom observed) 
	path_out: Text output of the composite spectrum: 
	comb_style: Three different combination styles: Weighted mean ('weighted_mean'),
													Jack Knife ('Jack_Knife'),
													Weighted Sample Variance ('weighted_sample_variance')

	Make_composite: Combine the data into a composite 1d optical rest-frame 1 dimensional spectrum
					Use 'save = true' to save a txt file to the 'path_out'

"""

import glob
import methods as st

Path = glob.glob('data/spec/*.fits')
path_out = 'output/CcQG_v1.txt'
comb_style = 'weighted_mean'
interp_step = 3

st.Make_composite(Path,path_out,interp_step,comb_style,save=True)


