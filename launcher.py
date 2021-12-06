#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: rouabah (original script), belaloui (corrections and doc.)

This is a quick example sript to launch the SEIQRDP fitting and simulation.

The region for which the algorithm is applied can be changed (don't forget
to change the population size!)

NOTE: Some Python IDEs (like Spyder) don't work well with multiprocessing.
      Run in a dedicated console if the script doesn't work.

PARAMETERS
----------
region_name: The region's name. (currently countries only.)
pop_size: The size of the region's population.
n_experiments: The number of runs. The results are then averaged. Must be > 1
                to estimate the prediction errors.
max_gen: The maximum number of generations the genetic algorithm goes through.
n_sim_days: The number of simulated days for the epidemic.
train_data_size: The number of days of data to train (fit) on.
use_mp: If True, all CPUs are used by the program. If False, only one is used.

--------------
Paper: Rouabah, Tounsi, Belaloui. Genetic algorithm with cross-validation-based
       epidemic model and application to the early diffusion of
       COVID-19 in Algeria, DOI: 10.1016/j.sciaf.2021.e01050.
"""
__version__ = '2.0'

# Everything can be done with an Experiment object.
from experiment import Experiment
from datetime import datetime

if __name__ == '__main__':
    # Specifying the region (country name) and the population size.
    region_name = 'Italy'
    pop_size = 60483054

    print("\n___________________ INFO ___________________\n")

    # time info dd/mm/YY H:M:S
    print("Date and time =", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

    # Parameters:
    n_experiments = 64  # Number of fittings: results are then averaged.
    max_gen = 20  # The maximum number of iterations in the algorithm.
    n_sim_days = 220  # Simulated days
    train_data_size = 30  # Training data size.
    use_mp = True  # To use multiprocessing (all CPUs) or not.

    # Information.
    print(f'\nSimulation of the epidemic (using online data) in {region_name}'
          f' with N={pop_size},'
          f' for {n_sim_days} days, with {n_experiments} experiments'
          f' and {max_gen} maximum generations. \n'
          f'It uses the first {train_data_size} days of the available data for'
          f' training.')

    print("\n___________________ START ___________________ \n")

    # Setting up the 'experiment'.
    # if use_mp=True, all the available CPUs are used.
    seiqrdp_exp = Experiment(region_name, pop_size, n_sim_days,
                             train_data_size, use_mp)
    # Start
    seiqrdp_exp.run(n_experiments, max_gen)

    # Show results: 'all' -> curves and data.
    seiqrdp_exp.show_results('all')

    print("\n___________________ Done !___________________ \n")
