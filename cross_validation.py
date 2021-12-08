# -*- coding: utf-8 -*-
"""
@author: belaloui

This script computes the average fitness value for a validation dataset at
each number of generations 'g'. Analyzing the fitness versus the number g
shows when the genetic algorithm starts to overfit the training data.
An optimum generations number 'g_opt' should give the best (lowest) fitness
value for the validation dataset. Therefore, when using this GA for
predictive purposes, one must stop the algorithm after 'g_opt' generations
to get a good model for prediction.
This is not necessary when using the GA only for fitting the model's
parameters, i.e., to estimate the epidemic parameters from epidemic data.
"""
__version__ = '2.0'

# from region_data import Region
from seiqrdp_suite import load_data
from seiqrdp_suite import GeneticFit
from seiqrdp_solver import SEIQRDPSolver
import matplotlib.pyplot as plt
import numpy as np

# Data

region_name = 'Italy'
pop_size = 60483054

data_size = 30
# Cross-validation
train_data_size = 22
validation_data_size = 8
# Genetic algorithm
init_n_gens = 1
max_n_gens = 100
#
n_fits = 64


def main():
    generations = [g for g in range(init_n_gens, max_n_gens+1)]
    fitnesses = [0]*len(generations)
    # Load all dataset
    region = load_data(region_name, pop_size)

    # Cutoff to the desired size (all data)
    region.cutoff_data(data_size)
    print(f'Loaded all (training + validation) dataset for {region.name} '
          f'from {region.first_day} to {region.last_day}.')

    # Taking out the validation data.
    cq_realdata = region.rcQ[-validation_data_size:]

    # Cutoff to the training data only
    region.cutoff_data(train_data_size)

    print(f'Using training dataset for {region.name} '
          f'from {region.first_day} to {region.last_day}.')

    for _ in range(n_fits):
        print(f'Fitting {_+1} of {n_fits}.')
        # First fit (from new random population)
        fit = GeneticFit(region=region, elite_size=10, max_gen=init_n_gens)
        fit.run()
        solution = fit.elite

        # Predicting data points (n=validation_data_size)
        model = SEIQRDPSolver(solution, region)
        model.solve_model(t_i=0, t_f=data_size)
        cq_prediction = model.cQ[-validation_data_size:]

        # first result
        fitnesses[0] += fit.n_rmse(cq_realdata, cq_prediction)

        # Evolve 1 generation at a time, and log the fitness.
        for g in range(1, len(generations)):
            # Evolve for 1 more generation.
            fit.run(resume=True, n_gens=1)
            solution = fit.elite

            # Predicting data points (n=validation_data_size)
            model = SEIQRDPSolver(solution, region)
            model.solve_model(t_i=0, t_f=data_size)
            cq_prediction = model.cQ[-validation_data_size:]

            # Adding the fitness value to the appropriate generation
            # number. The resulting sum over all fitnesses for generation g
            # will be divided by the number of fits later.
            fitnesses[g] += fit.n_rmse(cq_realdata, cq_prediction)

    # Computing the average fitness per generation.
    fitnesses = [f/n_fits for f in fitnesses]

    # Saving data to csv file.
    np.savetxt(fname='fit_vs_gen.csv', X=np.c_[generations, fitnesses],
               fmt=('%d', '%1.5f'),  delimiter=',')

    # Showing and saving plot.
    plt.plot(generations, fitnesses)
    plt.savefig(fname='fit_vs_gen.svg', format='svg')


if __name__ == '__main__':
    main()
