#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: rouabah
This is a quick sript to launch the SEIQRDP program
"""

from seiqrdp_model.Experiment import Experiment
from datetime import datetime

if __name__ == '__main__':
    print("\n___________________ INFO___________________ \n")

    # time info dd/mm/YY H:M:S
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(" date and time =", dt_string)

    # Parameters
    nExp = 12
    mGen = 20
    nDays = 90

    location = 'Italy'

    # location = 'Spain'

    print(f'\n This run simulates using internet available data of {location}, for {nDays} days, with {nExp} experiments and {mGen} maximum generations. \n ')


    # Run programs
    print("\n___________________ START___________________ \n")

    SEIQRDP_exp=Experiment(location,nDays, useMP=True)
    SEIQRDP_exp.run(nExp,mGen)

    # Show results
    SEIQRDP_exp.show_results()

    print("\n___________________ Done !___________________ \n")


    now = datetime.now()
    dt_ending = now.strftime("%d/%m/%Y %H:%M:%S")
    print(" date and time =", dt_ending)
