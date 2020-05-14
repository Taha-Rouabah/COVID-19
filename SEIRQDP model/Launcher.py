#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:05:56 2020

@author: rouabah

This is a quick sript to launch the programs for Italy
"""

from Experiment import Experiment
import numpy as np
import os
from datetime import datetime


print("\n___________________ INFO___________________ \n")
## Comments about the run
# dd/mm/YY H:M:S
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(" date and time =", dt_string)

print("\n This run simulates using data of Algeria, 01/03-29/04 with nExp=4, mG=500 \n ")
#print("\n This run simulates using data of Italy, 22/02-28/03 with nExp=8, mG=2k \n ")

nExp = 4
mGen = 500
location = 'Algeria'
# Input directory
#os.chdir('/Users/rouabah/MacWork/MyResearch/COVID19/Numerics/Data/Inputs')

# Programs directory
#os.chdir('/Users/rouabah/MacWork/MyResearch/COVID19/Numerics/Programs/09042020/v007.3')

# Run programs
print("\n___________________ START___________________ \n")

Al_exp=Experiment(location,180,180)
Al_exp.run(nExp,mGen)

# Show results
Al_exp.show_results()


# Save variables
# np.savetxt('Italy_Q_cQ_E_I_R_D_Rn_exp1000_mG1000_LSODA.txt', np.c_[Italy_exp.get_BaseExperiment().Q, Italy_exp.get_BaseExperiment().cQ, Italy_exp.get_BaseExperiment().E, Italy_exp.get_BaseExperiment().I, Italy_exp.get_BaseExperiment().R, Italy_exp.get_BaseExperiment().D, Italy_exp.get_BaseExperiment().Rn])

# Save parameters


# #############################################
# # Added by Taha
#         # Writting to a file
#         fid=open('Parmeters.txt','a');
#         fid.write(f'{Fore.GREEN}_______________________ Typical Values _______________________{Style.RESET_ALL}')
#         fid.write(f'{Fore.RED} Typical values are chosen to give the closet value to the mean BRN{Style.RESET_ALL}' )
#         fid.write('')
#         fid.write('BRN: ', round(self.typical_BRN,3));
#         if self.peak_found:
#             fid.write('PeakTime: ', self.city.first_day + timedelta(round(self.typical_peakTime)),\
#                   ' (', round(self.typical_peakTime,2),')');
#             fid.write('Peak:', round(self.typical_peak));
#         else:
#             fid.write("A peak could not be found in the chosen time period.")
#         fid.write('Protection rate: ', round(self.typical_parameters.alpha,3))
#         fid.write('Transmission rate: ', round(self.typical_parameters.beta,3))
#         fid.write('Latent period: ', round(1/self.typical_parameters.gamma,3))
#         fid.write('Infection period: ', round(1/self.typical_parameters.delta,3))
#         fid.write('Recovering rate: ', round(self.typical_parameters.lembda,3))
#         fid.write('Fatality rate: ', round(self.typical_parameters.lembda,3))
#         fid.write('Number of exposed people in (',self.city.first_day,') = ', round(self.typical_parameters.E0))
#         fid.write('Number of infective people in (',self.city.first_day,') = ', round(self.typical_parameters.I0))
#         fid.close();
# #############################################


# Analysis part

"""
# Suceptible to total population ratio
import matplotlib.pyplot as plt
StoN = []
for s in Al_exp.get_BaseExperiment().S:
    StoN.append(s/N)
    
plt.plot(StoN)
"""





