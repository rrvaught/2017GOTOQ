
import math
import numpy as np
import matplotlib.pyplot as pl
pl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.ticker import MaxNLocator
from astropy.table import Table
import json
from scipy import constants

#
def interp_json(infil,fit_cmp):

    #infil = '../../data/J1534+5015_model.json'
    #fit_cmp = 'z-0.00005_NaI'

    with open(infil) as data_file:
        data=json.load(data_file)

    systems=['z0.00000_MW']
    #components=[]

    cmp_dict = 0
    for cmp in data["cmps"]:
        systems.append(str(cmp))
        #print(data["cmps"][str(cmp)]['wrest'])
        if(fit_cmp == str(cmp)):
            cmp_dict=data["cmps"][str(cmp)]

    if(cmp_dict==0):
        print("Your selected system is not in this json file!")
        print("Pick one of the following:")
        print(systems)

        #cmp_data+=str(cmp)+';'
        #cmp_data+=str(cmp_dict["Reliability"])+';'
        #if cmp_dict["Comment"]=='':
        #    cmp_data+='None;'
        #else:
        #    cmp_data+=str(cmp_dict["Comment"])+';'
        #cmp_data+=str(cmp_dict["Nfit"])+';'
        #cmp_data+=str(cmp_dict["bfit"])+';'
        #components.append(cmp_data)
        #print(cmp_dict['zfit'])
        #print(cmp_dict['vlim'])

    if(cmp_dict==0):
        return {'zfit':0.0, 'vlim':0.0}
    else:
        
        zfit = cmp_dict['zfit']
        vlim = cmp_dict['vlim']

        lam0 = (1.0+zfit)*cmp_dict['wrest']
        lam0red = (1.0+zfit)*5897.5581

        print("igm_guesses results for this component:")
        print(lam0,lam0red,cmp_dict['Nfit'],cmp_dict['bfit'])    
        #print(cmp_dict.keys())    
        print("All systems in json file:")
        print(systems)
        #print(components)

    
        return {'zfit':zfit, 'vlim':vlim}
