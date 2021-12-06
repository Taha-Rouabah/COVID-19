"""
Version 8.0
***********
This version is object oriented and uses 5 methods via
SEIQRDPSolver to solve the SEIQRDP model.
It supports multiprocessing (with experiment module or any other method
                             by 'adding' experiments together: exp1 + exp2).
This version is not optimized.

In this version (8.0.alpha):
    - The dataset is now handled properly. No more global variable LOCATIONS.
    - Only the data for the population/region in question is loaded.
    - Dataset is downloaded with user's permission if no data file is found.
    - Corrected the fitness function (nRMSE as described in the paper).
    - The fitness function has been re-written to be readable.
    - Got rid of Comparator class. Fitness is handled in the GA class.
    - Got rid of Simulator class. SEIQRDPSolver does the same job.
    - *MAJOR* Dictionaries are now the data structure for the model's
        parameters. This will allow easily integration of other models.
    - Got rid of Configurator class and other parameter's handling functions.
        The parameters are accessed by their names in the new dictionary
        structure.
    - The model's search space is no longer hardcoded in this file. It was
        moved to the model's file.
    - Corrected the 95% confidence interval computation. (was off by 1%)
    - GeneticFit no longer fits at initialization. Requires a call to run().
    - Data is loaded more properly and much faster.
    - Renamed a lot of stuff (name conventions were not respected before)
    - Some code refactoring. It runs significantly faster.
"""
__version__ = '8.0.a'

# Math
import numpy as np
import scipy.stats as st
from math import sqrt
from random import random
from random import choice as rand_choice
from itertools import combinations
# Results
import matplotlib.pyplot as plt
from datetime import timedelta
# Data
import pandas as pd
from os.path import isfile, isdir
from os import mkdir
# Utils
from copy import deepcopy  # Is needed to create hard copies
# Internal
from region_data import Region
import seiqrdp_solver as slv


# ###############    DATA    ################
def update_data():
    """
    Downloads data from
    https://raw.githubusercontent.com/datasets/
    """
    dataset = pd.read_csv('https://raw.githubusercontent.com/datasets/'
                          'covid-19/master/data/'
                          'time-series-19-covid-combined.csv')
    if not isdir('Data/'):
        mkdir('Data/')
    dataset.to_csv(r'Data/time-series-19-covid-combined.csv',
                   index=False)


def load_data(region_name, pop_size):
    """
        Loading data from csv file.
        Returns the Region object.
    """
    print("Dataset is loading...")

    # If data is missing, prompt for data download.
    if not isfile('Data/time-series-19-covid-combined.csv'):
        print("Dataset not found.")
        usr_input = input('Download data from github repo? (y/any) :\n')
        if usr_input == 'y' or usr_input == 'Y':
            print("Dataset is being updated...")
            update_data()
    # Load data from file.
    dataframe = pd.read_csv('Data/time-series-19-covid-combined.csv')

    region = Region.from_dataframe(dataframe, region_name, pop_size)
    return region
# ###############    DATA    ################


class GeneticFit:
    """
        This class implements an evolutionary genetic algorithm that
        finds the best fitting SEIR parameters for a 'region' (pop + cases).
        At each generation we keep the best 'elite_size' parameters and breed
        them with a 'mutation_rate'.
        We stop at 'max_gen' generations.
        The best result (parameters) is 'self.elite'.
        Set 'verbose' = True to see the progress.
    """

    def __init__(self, region, elite_size, max_gen, mutation_rate=0.4,
                 search_space=slv.search_space, method='LSODA', verbose=False):
        self.region = region
        self.elite_size = elite_size
        self.max_gen = max_gen
        self.mutation_rate = mutation_rate
        self.search_space = search_space

        # I0 and E0 first guesses.
        self.search_space['I0'] = [region.rcQ[0], region.rcQ[0]*5]
        self.search_space['E0'] = [region.rcQ[0], region.rcQ[0]*5]

        self.method = method
        self.verbose = verbose

        (self.lembda, self.kappa) = self.compute_lambda_kappa()

        self.result = []
        self.elite = []

    def run(self):
        self.result = self.evolve()
        self.elite = self.result[0]

    def compute_lambda_kappa(self):
        """
            Determine Lambda and Kappa from real data.
            Note: 'lembda' is used instead of the Python reserved name lambda.
        """
        lembda = []
        kappa = []
        cQ = np.array(self.region.rcQ)
        R = np.array(self.region.rR)
        D = np.array(self.region.rD)
        Q = cQ - R - D

        for n in range(len(Q)-1):
            if (R[n] == 0) and (D[n] == 0):  # v007.4
                continue
            lembda.append((R[n+1]-R[n])/Q[n])
            kappa.append((D[n+1]-D[n])/Q[n])

        Lembda = np.mean(np.array(lembda))
        Kappa = np.mean(np.array(kappa))
        return (Lembda, Kappa)

    def evolve(self):
        """
        Evolves the initial population of parameters over a number of
        generations to fit the epidemic model.

        Returns:
        (parameters, fitness)
            A tuple of the best parameters set and its fitness value.

        """
        # Computing the total population size from the elite size.
        pop_size = self.elite_size*(self.elite_size-1)

        # Computing the search space's widths (in which params can vary)
        width = {}
        for param in self.search_space:
            width[param] = self.search_space[param][1]\
                - self.search_space[param][0]

        population = []
        elite = []

        # Generating the initial population
        for i in range(pop_size):
            parameters = {'lambda': self.lembda, 'kappa': self.kappa}
            for param in self.search_space:
                parameters[param] = random()*width[param]\
                    + self.search_space[param][0]

            # Adding the set of parameters to the population along with
            # its fitness.
            population.append([parameters, self.fitness(parameters)])

        # Selecting the elite
        elite = self.generate_elite(population, self.elite_size)

        # Main loop: iterates up to max_gen generations.
        for g in range(self.max_gen):
            # Passing the elite into the next population
            new_population = self.generate_next_population(elite, width)
            # Generate the next elite.
            elite = self.generate_elite(new_population, self.elite_size)

            if self.verbose:
                print(f'Generation {g+1} of {self.max_gen}.'
                      f'Best: {elite[0][1]}\n')
        # Return a tuple (params, fitness) of the best elite.
        return ([elite[0][0], elite[0][1]])

    def generate_next_population(self, elite, width):
        """
        Generates the next generation's population from the previous ones'
        elite.

        Parameters:
            elite :
                The elite parameters which survive and breed.
            width :
                List of widths of the search space; by how much the parameters
                can mutate.

        Returns:
            new_population : list
                The population of the next generation.
        """
        # Passing the elite directly into the next population
        new_population = deepcopy(elite)

        # Breed and mutate, and add to the new population
        for par1, par2 in combinations(elite, 2):
            offspring = self.cross(par1[0], par2[0],
                                   width, self.search_space)
            # Add the offspring to the new population along with
            # its fitness.
            new_population.append([offspring, self.fitness(offspring)])

        return new_population

    def sort_population(self, pop):
        """
            Sorts the population of the genetic pool by
            their fitting score (2nd element in population list)
            Note: get rid of lambda and replace it with a externally
                    declared function (performance issue)
        """
        # reverse = None (Sorts in Ascending order)
        pop.sort(key=lambda x: x[1])
        return pop

    def cross(self, par1, par2, width, search_space):
        """
            Generates an offspring by breeding 2 parents, with
            a 'mutation_chance' = 0-1.
        """
        offspring = {'lambda': self.lembda, 'kappa': self.kappa}
        for param in search_space:
            if random() < self.mutation_rate:
                #   Mutation: the gene is random
                offspring[param] = random()*width[param]+search_space[param][0]
            else:
                #   Breeding: the gene is from a single random parent
                offspring[param] = rand_choice([par1[param], par2[param]])
        return offspring

    def generate_elite(self, population, elite_size):
        """
            Sorts the population by their scores and
            returns 'elite_size' best elements.
        """
        # Selecting the elite
        self.sort_population(population)
        return population[:elite_size]

    def fitness(self, genome):
        """ Computes the fitness of a genome.
            genome = [6 values], no kappa/lambda.
        """
        simulation = slv.SEIQRDPSolver(genome, self.region, mode='fit')

        simulation.solve_model(0, len(self.region.rcQ), 1, self.method)

        sqr_sum = 0
        for a, b in zip(self.region.rcQ, simulation.cQ):
            sqr_sum += (a-b)**2
        n_rmse = sqrt(np.mean(sqr_sum))/np.mean(self.region.rcQ)
        return n_rmse


class Analyzer:
    """
        This class measures BRN, peak_time, Peak values for a given 'Region'
        and a set of SEIR parameters.
        It gives also the data of simulation for n_sim_days period after
        the first day.
    """

    def __init__(self, region, parameters, n_sim_days):
        self.region = region
        self.parameters = parameters
        self.n_sim_days = n_sim_days

        simulation = slv.SEIQRDPSolver(parameters, region)
        simulation.solve_model(0, self.n_sim_days)

        active_cases = simulation.Q

        peak_active_infections = max(simulation.Q)
        peak_day = simulation.Q.index(peak_active_infections)
        self.peak_is_outofrange = False
        if peak_day == len(simulation.Q)-1:
            self.peak_is_outofrange = True

        self.peak = active_cases[peak_day-1]
        self.peak_time = peak_day
        self.brn = parameters['beta']/parameters['delta']
        self.fit = simulation

    def get_epidemic_startdate(self):
        """
        This can be used to seek the probable first day when the
        first exposition occured.
        """

        solver = slv.SEIQRDPSolver(self.parameters, self.region)

        step = -30
        solver.solve_model(0, step, 1, 'LSODA')
        E = solver.E

        while E[0] > 1:
            step += -30
            # Simulates the previous 30 days
            solver.solve_model(step+30, step, 1, 'LSODA')
            E = solver.E
        starting_point = 0

        while E[starting_point] < 0:
            starting_point += 1

        starting_point += step
        starting_date = self.region.first_day\
            + timedelta(days=starting_point)
        return(starting_date)


class BaseExperiment:
    """
        A BaseExperiment perfoms a genetic fitting + simulation
        and computes the model variables.
        It is suited for parallel computing since 2 instences of
        BaseExeriment can be added together to form a 3rd one:
            exp3 = exp1 + exp2
        Example:
            >> import Algeria_SEIR_COVID2019_Object as SEIR
            >> SEIR.load_data()
            >> ita_exp = SEIR.BaseExperiment(region=SEIR.LOCATIONS['Italy'],
                                             n_sim_days=100, n_experiments=4)
            >> ita_exp.run()
            >> ita_exp.compute_results()
            >> ita_exp.show_results('all')
    """

    def __init__(self, region, n_sim_days, n_experiments,
                 elite_size=10, max_gen=50, mutation_rate=0.4,
                 search_space=slv.search_space, method='LSODA', verbose=False):

        self.region = region
        self.n_sim_days = n_sim_days
        self.n_experiments = n_experiments
        self.elite_size = elite_size
        self.max_gen = max_gen
        self.mutation_rate = mutation_rate
        self.search_space = search_space
        self.method = method
        self.verbose = verbose

    def init_states(self):

        # typical values
        self.typical_brn = 0
        self.typical_experiment_index = 0
        self.typical_parameters = {}
        self.typical_peak = 0
        self.typical_peak_time = 0
        self.typical_peak_found = True

        # observables mean values
        self.brn = 0
        self.peak = 0
        self.peak_time = 0
        self.peak_found = True

        # observables covariances
        self.err_brn = 0
        self.err_peak = 0
        self.err_peak_time = 0

        # Elite parameters
        self.elites = []

        # All exps. data.
        self.raw_brn, self.raw_Rn, self.raw_peak, self.raw_peak_time,\
            self.raw_S, self.raw_E, self.raw_I, self.raw_Q, self.raw_cQ,\
            self.raw_R, self.raw_D, self.raw_P\
            = [], [], [], [], [], [], [], [], [], [], [], []

        # mean compartmental variables
        self.S, self.E, self.I, self.Q, self.cQ, self.R, self.D,\
            self.P, self.Rn\
            = [], [], [], [], [], [], [], [], []

        # compartmental variables covariances
        self.err_S, self.err_E, self.err_I, self.err_Q, self.err_cQ,\
            self.err_R, self.err_D, self.err_P, self.err_Rn\
            = [], [], [], [], [], [], [], [], []

        # days list of simulation
        self.days = []

    def run(self):
        """
            Runs the simulations WITHOUT computing the final results.
        """
        # Here we perform n_experiments of experiments
        # (fitting + calculating observables)
        self.init_states()

        for n in range(self.n_experiments):
            # Fit real data and obtaining the n'th elite parameter
            if self.verbose:
                print("Fitting data: ", n+1, " out of ",
                      self.n_experiments, ' experiments.')
            fit = GeneticFit(self.region, self.elite_size, self.max_gen,
                             self.mutation_rate, self.search_space,
                             self.method, self.verbose)
            fit.run()

            # Calculate observables of n'th elite parameter
            if self.verbose:
                print("Calculating epidemic data.")
            observables = Analyzer(self.region, fit.elite, self.n_sim_days)
            # Assign what you have obtained to lists
            self.elites.append(fit.elite)
            self.raw_brn.append(observables.brn)
            self.raw_peak.append(observables.peak)
            self.raw_peak_time.append(observables.peak_time)
            self.raw_S.append(observables.fit.S)
            self.raw_E.append(observables.fit.E)
            self.raw_I.append(observables.fit.I)
            self.raw_Q.append(observables.fit.Q)
            self.raw_cQ.append(observables.fit.cQ)
            self.raw_R.append(observables.fit.R)
            self.raw_D.append(observables.fit.D)
            self.raw_P.append(observables.fit.P)
            # R(t) = brn * S/N  v007.4
            self.raw_Rn.append([observables.brn * s/self.region.N
                                for s in self.raw_S[n]])
            if observables.peak_is_outofrange:
                self.peak_found = False

    def compute_results(self):
        """
            Computes the final results of the simulations.
        """
        # Calculate the {mean values} and {typical values}
        # of observables from n_experiments of experiments
        if self.verbose:
            print("I am calculating expected values")

        self.brn = np.mean(np.array(self.raw_brn))

        self.typical_brn = self.raw_brn[0]
        for n in range(1, self.n_experiments):
            if (self.raw_brn[n]-self.brn)**2 < (self.typical_brn-self.brn)**2:
                self.typical_brn = self.raw_brn[n]
                self.typical_experiment_index = n
            else:
                continue

        self.typical_parameters = self.elites[self.typical_experiment_index]

        self.peak = np.mean(np.array(self.raw_peak))
        self.typical_peak = self.raw_peak[self.typical_experiment_index]

        self.peak_time = np.mean(np.array(self.raw_peak_time))
        self.typical_peak_time\
            = self.raw_peak_time[self.typical_experiment_index]

        # Calculate the mean values of compartmental variables
        for i in range(self.n_sim_days):
            self.S.append(np.mean(np.array(self.raw_S)[:, i]))
            self.E.append(np.mean(np.array(self.raw_E)[:, i]))
            self.I.append(np.mean(np.array(self.raw_I)[:, i]))
            self.Q.append(np.mean(np.array(self.raw_Q)[:, i]))
            self.cQ.append(np.mean(np.array(self.raw_cQ)[:, i]))
            self.R.append(np.mean(np.array(self.raw_R)[:, i]))
            self.D.append(np.mean(np.array(self.raw_D)[:, i]))
            self.P.append(np.mean(np.array(self.raw_P)[:, i]))
            self.Rn.append(np.mean(np.array(self.raw_Rn)[:, i]))

        # Calculate the error of observables
        # Critical value for 95% confidence interval.
        crit_value = st.t.ppf((1 + 0.95) / 2, self.n_sim_days - 1)
        if self.verbose:
            print("Calculating errors")
        self.err_brn = st.sem(np.array(self.raw_brn)) * crit_value
        self.err_peak = st.sem(np.array(self.raw_peak)) * crit_value
        self.err_peak_time = st.sem(np.array(self.raw_peak_time)) * crit_value

        # Calculate the std error on comartmental variables  v007.4
        #  * crit_value means sigma2 : 95% confidence
        for i in range(self.n_sim_days):
            self.days.append(i)
            self.err_S.append(st.sem(np.array(self.raw_S)[:, i]) * crit_value)
            self.err_E.append(st.sem(np.array(self.raw_E)[:, i]) * crit_value)
            self.err_I.append(st.sem(np.array(self.raw_I)[:, i]) * crit_value)
            self.err_Q.append(st.sem(np.array(self.raw_Q)[:, i]) * crit_value)
            self.err_cQ.append(st.sem(np.array(self.raw_cQ)[:, i])*crit_value)
            self.err_R.append(st.sem(np.array(self.raw_R)[:, i]) * crit_value)
            self.err_D.append(st.sem(np.array(self.raw_D)[:, i]) * crit_value)
            self.err_P.append(st.sem(np.array(self.raw_P)[:, i]) * crit_value)
            self.err_Rn.append(st.sem(np.array(self.raw_Rn)[:, i])*crit_value)

    def show_results(self, *curves):
        """
            Shows the curves and prints the results of
            n_experiments experiments.
        """
        print('___________________ Plots ___________________')
        self.show_curves(*curves)

        self.show_result_values()
        return

    def show_curves(self, *curves):
        """
            Shows the curves only.
            Available curves:
                cq: Cumulative quarantined.
                q: active quarantined.
                rn: Rn(t)
            Choosing 'all' or no arguments shows all the curves.
        """
        # Errors handling
        args = ['cq', 'q', 'rn', 'all']
        n_plots = 0
        for arg in curves:
            if arg not in args:
                print(f'Warning: \'{arg}\' is not a valid curve.')
            else:
                n_plots += 1

        draw_all = (len(curves) == 0) or ('all' in curves)
        if draw_all:
            n_plots = 3

        # Plot cQ and Q with error bars
        fig, plots = plt.subplots(nrows=n_plots, ncols=1, sharex=False,
                                  figsize=[8, (13/3)*n_plots])

        if n_plots == 1:
            plots = (plots,)

        last_plot = -1
        if 'cq' in curves or draw_all:
            last_plot += 1
            cq = plots[last_plot]
            cq.set_title('Cumulative Number of cases')
            cq.plot(self.region.rcQ, color='red', marker="o",
                    label='Real data (fitting)')
            cq.errorbar(self.days, self.cQ, yerr=self.err_cQ, errorevery=6,
                        ecolor='black', label='Model')
            cq.legend(loc="lower right")
            cq.set_xlabel(f'Days since {self.region.first_day}')
            cq.set_ylabel('cQ(t)')
            cq.grid(True)

        if 'q' in curves or draw_all:
            last_plot += 1
            q = plots[last_plot]
            q.set_title('Number of active cases')
            q.errorbar(self.days, self.Q, yerr=self.err_Q, errorevery=6,
                       ecolor='black')
            q.set_xlabel(f'Days since {self.region.first_day}')
            q.set_ylabel('Q(t)')
            q.grid(True)

        if 'rn' in curves or draw_all:
            last_plot += 1
            rn = plots[last_plot]
            rn.set_title('Reproduction Number')
            rn.errorbar(self.days, self.Rn, yerr=self.err_Rn, errorevery=6,
                        ecolor='black', label='R(t)')
            # Drawing Rn = 1
            rn.plot(self.days, [1]*len(self.Rn), label='R=1', color='red',
                    linestyle='dashed')
            rn.legend(loc="upper right")
            rn.set_xlabel(f'Days since {self.region.first_day}')
            rn.set_ylabel('R(t)')
            rn.grid(True)

        fig.tight_layout()
        plt.show()

    def show_result_values(self):
        """
            Shows the numerical values of the experiment.
        """
        # Showing mean values
        print('_______________________ Mean Values _______________________')
        print('BRN: ', round(self.brn, 3), "\u00B1", round(self.err_brn, 3))

        if self.peak_found:
            print('peak_time: ',
                  self.region.first_day + timedelta(round(self.peak_time)),
                  ' (', round(self.peak_time, 2), ") \u00B1",
                  round(self.err_peak_time, 2))
            print('Peak:', int(self.peak), "\u00B1", int(self.err_peak))
        else:
            print("A peak could not be found in the chosen time period.")
        # Totals (mean)
        print('')
        print(f'Total infected (t={self.n_sim_days}):',
              int(self.cQ[-1]), "\u00B1", int(self.err_cQ[-1]))
        print(f'Total recovered (t={self.n_sim_days}):',
              int(self.R[-1]), "\u00B1", int(self.err_R[-1]))
        print(f'Total deaths (t={self.n_sim_days}):',
              int(self.D[-1]), "\u00B1", int(self.err_D[-1]))

        # Showing typical values
        print('_______________________ Typical Values _______________________')

        print('Typical values are chosen to give the closest BRN '
              'to the mean BRN.\n')

        print('BRN: ', round(self.typical_brn, 3))
        if self.peak_found:
            print('peak_time: ', self.region.first_day
                  + timedelta(round(self.typical_peak_time)),
                  ' (', round(self.typical_peak_time, 2), ')')
            print('Peak:', int(self.typical_peak))
        else:
            print("A peak could not be found in the chosen time period.")

        # Totals (typical)
        print('')
        print(f'Total infected (t={self.n_sim_days}):',
              int(self.raw_cQ[self.typical_experiment_index][-1]))
        print(f'Total recovered (t={self.n_sim_days}):',
              int(self.raw_R[self.typical_experiment_index][-1]))
        print(f'Total deaths (t={self.n_sim_days}):',
              int(self.raw_D[self.typical_experiment_index][-1]))

        print('')
        print('Protection rate: ',
              round(self.typical_parameters['alpha'], 3))
        print('Transmission rate: ',
              round(self.typical_parameters['beta'], 3))
        print('Latent period: ',
              round(1/self.typical_parameters['gamma'], 3))
        print('Infectious period: ',
              round(1/self.typical_parameters['delta'], 3))
        print('Recovery rate: ',
              round(self.typical_parameters['lambda'], 3))
        print('Fatality rate: ',
              round(self.typical_parameters['kappa'], 3))
        print('Number of exposed people in (', self.region.first_day, ') = ',
              round(self.typical_parameters['E0']))
        print('Number of infectious people in (', self.region.first_day,
              ') = ', round(self.typical_parameters['I0']))

    def __add__(self, other):
        """
            Defines the + operator for BaseExperiment.
            NOTE : This is just a hack for now!
            To Do: Make this more robust
        """
        newExp = deepcopy(self)
        newExp.elites += other.elites
        newExp.raw_brn += other.raw_brn
        newExp.raw_peak += other.raw_peak
        newExp.raw_peak_time += other.raw_peak_time
        newExp.raw_S += other.raw_S
        newExp.raw_E += other.raw_E
        newExp.raw_I += other.raw_I
        newExp.raw_Q += other.raw_Q
        newExp.raw_cQ += other.raw_cQ
        newExp.raw_R += other.raw_R
        newExp.raw_D += other.raw_D
        newExp.raw_P += other.raw_P
        newExp.raw_Rn += other.raw_Rn
        newExp.n_experiments += other.n_experiments

        # observables mean values
        newExp.brn = 0
        newExp.peak = 0
        newExp.peak_time = 0
        newExp.peak_found = newExp.peak_found and other.peak_found

        # observables covariances
        newExp.err_brn = 0
        newExp.err_peak = 0
        newExp.err_peak_time = 0

        # mean compartmental variables
        newExp.S, newExp.E, newExp.I, newExp.Q, newExp.cQ, newExp.R,\
            newExp.D, newExp.P, newExp.Rn = [], [], [], [], [], [], [], [], []

        # compartmental variables covariances
        newExp.err_S, newExp.err_E, newExp.err_I, newExp.err_Q, newExp.err_cQ,\
            newExp.err_R, newExp.err_D, newExp.err_P, newExp.err_Rn = [], [],\
            [], [], [], [], [], [], []

        newExp.days = []

        return newExp
