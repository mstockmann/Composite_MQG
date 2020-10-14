# Composite_MQG


I present the reposity containing the X-Shooter z=2 compact Quiescent Galaxy composite spectrum. All the code to analyse/stack the spectra and generate the plot are in the py directory. The data presented here are from ESO observations using the X-Shooter Echelle spectrograph on the Very Large Telescope (VLT) in Atacama, Chile. The data used here (data/spec/*.fits) are from the ApJ published paper "X-shooter Spectroscopy and HST Imaging of 15 Massive Quiescent Galaxies at z > 2" and can be found here: https://iopscience.iop.org/article/10.3847/1538-4357/ab5af4/pdf


This software is publicly available on Github (https://github.com/mstockmann/Composite_MQG) so feel free to use the software for inspiration. If you're interested in the data please contact me for more information. 

The software was produced and published to allow people to access and reuse the developed methods and to strive towards and open source mindset for scientific codes. The code is meant for scientific purpose however can be utilised in other data analysis projects. 

The code is run by the main python file composite_spec.py which initiates a series of functions saved in methods.py. This saves txt output files (when specified in the code), which can be visualized via the included python plotting code, plot_composite.py. The data input files can be found in the folder data that includes each galaxy spectrum, its galaxy evolution model and the corresponding spectroscopic redshift. 

Note that no unit testing is implemented here and thus the format should specifically match the files you see in the folders below. The arXiv version of the paper can be found here: https://arxiv.org/pdf/1912.01619.pdf

- Mikkel Stockmann
