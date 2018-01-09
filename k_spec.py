# Spectral Analysis of 2017 Keck Spectra by Ryan Rickards Vaught
# Makes use of the linetools library.

from linetools.spectra.xspectrum1d import XSpectrum1D
from astropy.coordinates import SkyCoord
from linetools.isgm import abssystem as lt_absys
from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools import line_utils as ltlu
import json
from linetools.spectralline import AbsLine, SpectralLine
import numpy
import matplotlib.pyplot as plt
from astropy.io import fits
from glob import glob
from astropy import units as u
from linetools.lists.linelist import LineList
import warnings
warnings.filterwarnings('ignore')

def get_position(iden):
    targ_ra_hr=name[48:50]
    targ_ra_min=name[50:52]
    targ_dec_hr=name[53:55]
    targ_dec_min=name[55:57]
    targ_RA=targ_ra_hr+':'+targ_ra_min
    targ_DEC=targ_dec_hr+':'+targ_dec_min
    return targ_RA,targ_DEC

 
# Now let's plot all the spectra and save each figure
fname=glob('/Users/ryan/Desktop/keck_spectra/2017mar5/*F.fits')
for name in fname:
    sp=XSpectrum1D.from_file(name)
    targ_coord=get_position(name[42:-5])
    targ_coord=get_position(name)
    
    dtype=[('RA','U5'),('DEC','U5'),('zq',float),('z',float),('b',float),('logM',float),('QSO_mag',float),('Sod_line',float)]
    prop=numpy.genfromtxt('/Users/ryan/Desktop/keck_spectra/2017mar5/target_list.dat',usecols = (0,1,2,3,4,5,6,7),dtype=dtype)
    
    # Match the target with its properties by using RA and DEC cooords
    w = numpy.logical_and(targ_coord[0]==prop['RA'] , targ_coord[1]==prop['DEC'] )
    data=prop[w]
    
    radec=SkyCoord(targ_coord[0]+':00', targ_coord[1]+':00', unit=(u.hourangle, u.deg))
    
    # Get the expected sodium transition wavelengths from the target properties
    line=data['Sod_line']
    
    #Define the wavelength range to search for redshifted Sodium lines
    wvlim=[line-30.0,line+30.0]*u.AA 
    
    strong=LineList('Strong') #Search for lines 
    transitions= strong.available_transitions(wvlim /(1+data['z']), n_max_tuple= None, min_strength=0.0)
    
    # Transform lines to the redshifted wavelength
    line1=transitions['wrest'][0]*(1+data['z'])
    line2=transitions['wrest'][1]*(1+data['z'])
    
    # plot the figure
    fig=plt.figure()
    fig.suptitle(name[42:-5]+'  z='+str(data['z']))
    plt.axvline(x=line1,color='k', linestyle='--')
    plt.axvline(x=line2,color='k', linestyle='--')
    sp.plot(xlim=(line-30,line+30))
    fig.savefig(name+'.pdf') # Save a copy of the figure
    
    
    sp.fit_continuum(kind='QSO', redshift=data['zq'])
    sp.normalize(co=sp.co)
    fig1=plt.figure()
    fig1.suptitle(name[42:-5]+'  z='+str(data['z']))
    plt.axvline(x=line1,color='k', linestyle='--')
    plt.axvline(x=line2,color='k', linestyle='--')
    sp.plot(xlim=(line-20,line+20))
    fig1.savefig('n_'+name[42:-5]+'.pdf')
    sp.write_to_fits('n_'+name[42:-5]+'.fits')  
    
   
      # Calculate the EW of 1st sodium line
    
    S1 = AbsLine(transitions['wrest'][0]*u.AA,z=float(data['z']))
    S1.analy['spec']= XSpectrum1D.from_file('n_'+name[42:-5]+'.fits')
    S1.analy['vlim'] = [-150.,150.]*u.km/u.s
    S1.attrib['coord']=radec
    x = float(input("Enter a lower lim: "))
    y = float(input("Enter a higher lim: "))
    
    S1.limits.set([x,y]*u.AA)
    S1.measure_aodm()
    S1.measure_ew() # Observer frame
    S1.measure_kin()
    N, sigN, flgN = [S1.attrib[key] for key in ['N','sig_N','flag_N']] 
    print(S1.attrib)
    abscomp = AbsComponent.from_abslines([S1])
    abscomp.stack_plot()
    
    # Calculate the EW of 2nd sodium line
    S2 = AbsLine(transitions['wrest'][1]*u.AA, z=float(data['z']))
    S2.attrib['coord']=radec
    S2.analy['vlim'] = [-150.,150.]*u.km/u.s
    S2.analy['spec']= XSpectrum1D.from_file('n_'+name[42:-5]+'.fits')
    x = float(input("Enter a lower lim: "))
    y = float(input("Enter a higher lim: "))
    S2.limits.set([x,y]*u.AA)
    S2.measure_ew() # Observer frame
    N, sigN, flgN = [S2.attrib[key] for key in ['N','sig_N','flag_N']] 
    
    

    


print('DONE')            
#adict = specline.to_dict()
#with io.open(outfil, 'w', encoding='utf-8') as f:
#   f.write(unicode(json.dumps(tdict, sort_keys=True,
#      indent=4, separators=(',', ': '))))