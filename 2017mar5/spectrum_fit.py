# Spectral Analysis of 2017 Keck Spectra by Ryan Rickards Vaught
# Use is to analyze a given spectrum
# Makes use of the linetools library.

def get_prop(filename):
	import re
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	import numpy
	
	#Extract the RA and DEC from the filename
	coords=[int(s) for s in re.findall(r'[-]?\d', filename)]
	targ_RA=str(coords[0])+str(coords[1])+':'+str(coords[2])+str(coords[3])
	targ_DEC=str(coords[4])+str(coords[5])+':'+str(coords[6])+str(coords[7])
	targ_coord=targ_RA, targ_DEC
	
	dtype=[('RA','U5'),('DEC','U5'),('zq',float),('z',float),('b',float),('logM',float),
		('QSO_mag',float),('Sod_line',float)]
	
	prop=numpy.genfromtxt('/Users/ryan/Desktop/keck_spectra/2017mar5/target_list.dat',
		usecols = (0,1,2,3,4,5,6,7),dtype=dtype)
		
	w = numpy.logical_and(targ_coord[0]==prop['RA'] , targ_coord[1]==prop['DEC'] )
	data=prop[w]
	radec=SkyCoord(targ_coord[0]+':00', targ_coord[1]+':00', unit=(u.hourangle, u.deg))
	return data, radec

def fit_cont(filename):

	from linetools.spectra.xspectrum1d import XSpectrum1D
	import re
	from astropy.coordinates import SkyCoord
	from linetools.isgm import abssystem as lt_absys
	from linetools.spectralline import AbsLine
	from linetools.isgm.abscomponent import AbsComponent
	from linetools import line_utils as ltlu
	from linetools.spectralline import AbsLine, SpectralLine
	import numpy
	from astropy.io import fits
	from glob import glob
	from astropy import units as u
	from linetools.lists.linelist import LineList
	import warnings
	warnings.filterwarnings('ignore')
	
	# Create a spectrum class object from the fits file
	sp=XSpectrum1D.from_file(filename)
	 

	data, radec=get_prop(filename)
	
	# Call the GUI to interactively fit the continuum
	sp.fit_continuum(kind='QSO', redshift=data['zq'])
	
	# Normalize the Continuum
	sp.normalize(co=sp.co)
	
	# Write the Normalized Continuum to a new fits file 
	sp.write_to_fits('n_'+ filename)
	
	return print( filename+' has been normalized and written to ' +'n_'+ filename) 
	
def spec_inspect(filename):
	import re
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	
	data,radec = get_prop(filename)
	
	
	# Launch the GUI xspec
	import os
	os.system('lt_xspec --norm --z '+str(float(data['z']))+' n_'+filename)
	
	return print('Closed lt_xspec')
	
def read_json(file):
	import json
	from pprint import pprint
	
	with open(file) as data_file:
		data = json.load(data_file)
	return data
         
def de_redshift(spectrum,z):
	return spectrum.wavelength/(1.0+z)
	
def spec_dered(filename):
	from linetools.lists.linelist import LineList
	from astropy import units as u
	from linetools.spectra.xspectrum1d import XSpectrum1D
	
	# Create a json file and centroid the lines
	spec_inspect(filename)
	
	# Read in the json file
	data=read_json(filename[:-5]+'.json')
	z_new=data['zabs']
	
	# Create a spectrum object, and normalize it, using the 
	# already fitted continuum.
	
	sp=XSpectrum1D.from_file('n_'+ filename)
	sp.normalize(co=sp.co)
	
	# Transform the spectra into the rest frame
	sp.wavelength=de_redshift(sp,z_new)
	sp.write_to_fits('RF_'+filename) 
	sp.plot()
	
	return print('The Spectrum has been transformed into the Lab Frame')
	
def measure_sodium_EW(filename):
	from linetools.lists.linelist import LineList
	from astropy import units as u
	from linetools.spectralline import AbsLine
	from linetools.spectra.xspectrum1d import XSpectrum1D
	import matplotlib.pyplot as plt
	
	sp=XSpectrum1D.from_file('RF_'+ filename)
	sp.normalize(co=sp.co)
	
	wvlim=[5880,6030]*u.AA
	strong=LineList('Strong')
	transitions= strong.available_transitions(wvlim, n_max_tuple= None, min_strength=0.0)
	line1=transitions['wrest'][0]
	line2=transitions['wrest'][1]
	avg_line=(line1+line2)/2.0
	
	# Plot the spectrum to get limits for EW
	fig=plt.figure()
	plt.axvline(x=line1,color='k', linestyle='--')
	plt.axhline(y=1.0,color='r', linestyle='--')
	plt.axvline(x=line2,color='k', linestyle='--')
	sp.plot(xlim=(avg_line-30,avg_line+30))
	
	S1 = AbsLine(transitions['wrest'][0]*u.AA,z=0.0)
	S1.analy['spec']= sp
	S2 = AbsLine(transitions['wrest'][1]*u.AA,z=0.0)
	S2.analy['spec']= sp
	
	#x = float(input("Enter a lower lim: "))
	#y = float(input("Enter a higher lim: "))
	
	x=5888
	y=5896
	
	S1.limits.set([x,y]*u.AA)
	S1.measure_ew(flg=1) # Measure the EW of the first line
	EW1=S1.attrib['EW'], S1.attrib['sig_EW']
	
	
	#x = float(input("Enter a lower lim: "))
	#y = float(input("Enter a higher lim: "))
	
	x=5895
	y=5905
	S2.limits.set([x,y]*u.AA)
	S2.measure_ew() # Measure the EW of the second line
	EW2=S2.attrib['EW'],S2.attrib['sig_EW']
	return EW1,EW2
    
def measure_ew_lim(filename):
	from linetools.lists.linelist import LineList
	import numpy
	from astropy import units as u
	from linetools.spectralline import AbsLine
	from linetools.spectra.xspectrum1d import XSpectrum1D
	
	sp=XSpectrum1D.from_file('RF_'+ filename)
	sp.normalize(co=sp.co)
	
	wv0=5891.5833*u.AA
	R=4013.0/.3
	d_lambda_pix=0.13333317282580992
	d_lambda_res=5891.5833/R
	N_pix=d_lambda_res/d_lambda_pix
	s_2_n=sp.get_local_s2n(wv0, npix=100)[0] #/100 S_2_N per pixel
	EW_lim=3.0*numpy.sqrt(N_pix)*d_lamba_pix/s_2_n # 1 sigma limit
	return EW_lim
	