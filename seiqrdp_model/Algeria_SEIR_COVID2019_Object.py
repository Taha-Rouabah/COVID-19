"""
Version 7.5
***********
This version is object oriented and uses 5 methods via
SEIR_Solver to solve the SEIR model.
It supports multiprocessing (with Experiment module or any other method
                             by 'adding' experiments together: exp1 + exp2).
This version is not optimized.
"""
__version__ = '7.5'

from math import sqrt
from random import random
import matplotlib.pyplot as pl
from colorama import Fore, Style
import numpy as np
import scipy.stats as sp
from copy import deepcopy  # Is needed to create hard copies
from datetime import timedelta
import pandas as pd
from os.path import isdir
from os import mkdir

import SEIQRDP_model.SEIR_Solver as slv
from SEIQRDP_model.DataCollector import WebDataReader

# Our SEIR parameters are bound to these intervals
# They are: [0:alpha, 1:beta, 2:gamma, 3:delta, 4:I0, 5:E0 ]
ref_intervals = [[0, 0.15], [0.2, 1.5], [0.3, 0.6],
                 [1/7.5, 1/3.5], [0, 50], [0, 50]]


# ###############    DATA    ################

# PLEASE update load_data() when you add new regions.
LOCATIONS = {}

def update_data():
    """
    Uploads data from 
    https://raw.githubusercontent.com/datasets/
    """
    dataset = pd.read_csv('https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv')
    dataset.to_csv(r'Data/time-series-19-covid-combined.csv',
                   index=False)
    
    
def load_data():
    """
        Loading data from csv file into the LOCATIONS dict.
    """
    # Load data from source
    print("Dataset is loading...")
    if not isdir('Data/'):
        mkdir('Data/')
        update_data()

    dataset = pd.read_csv('Data/time-series-19-covid-combined.csv')
    


    # Austria
    austria = WebDataReader(dataset, 'Austria', custom_threshold=1).Region
    austria.N = 8994481  # Worldometer.info
    LOCATIONS[austria.name] = austria

    # Italy
    italy = WebDataReader(dataset, 'Italy', custom_threshold=1).Region
    italy.N = 60483054  # Worldometer.info
    LOCATIONS[italy.name] = italy  # 30

    # Spain
    spain = WebDataReader(dataset, 'Spain', custom_threshold=1).Region
    spain.N = 46754778  # Worldometer.info
    LOCATIONS[spain.name] = spain  # 24

    # Hubei
    hubei = WebDataReader(dataset, 'China', 'Hubei', custom_threshold=1).Region
    hubei.N = 59170000  # statista.com
    LOCATIONS[hubei.name] = hubei

    print(list(LOCATIONS.keys()))

    return

# ###############    DATA    ################


class Configurator:
    """
        Maps a list of SEIR parameters to a human-readable object where
        the parameters are refered to by their names.
    """

    def __init__(self, p=[0, 0, 0, 0, 0, 0, 0, 0]):
        if (len(p) == 8):
            self.alpha = p[0]
            self.beta = p[1]
            self.gamma = p[2]
            self.delta = p[3]
            self.lembda = p[4]
            self.kappa = p[5]
            self.I0 = p[6]
            self.E0 = p[7]
        else:
            print("Wrong configuration")

    def __repr__(self):
        return("{%s,%s,%s,%s,%s,%s,%s,%s}" %
               (self.alpha, self.beta, self.gamma, self.delta,
                self.lembda, self.kappa, self.I0, self.E0))

    def __str__(self):
        return("{%s,%s,%s,%s,%s,%s,%s,%s}" %
               (self.alpha, self.beta, self.gamma, self.delta,
                self.lembda, self.kappa, self.I0, self.E0))

    def values(self):
        return([self.alpha, self.beta, self.gamma, self.delta,
                self.lembda, self.kappa, self.I0, self.E0])


class Simulator:
    """
        Solves the SEIR model in a 'location' for 'parameters' during 'nDays'.
        In this version we use the SEIRSolver tool.
    """

    def __init__(self, location, parameters, nDays, method='LSODA', mode='sim'):
        """
            method = 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA' // DOP??
        """
        solver = slv.SEIRSolver(parameters, location, mode)
        solver.modelSol(0, nDays, 1, method)
        self.S = solver.S
        self.E = solver.E
        self.I = solver.I
        self.Q = solver.Q
        self.cQ = solver.cQ
        self.R = solver.R
        self.D = solver.D
        self.P = solver.P


class Comparator:
    """
        This class provides a tool to compare two set of data
        by calculating the normalized mean-square deviation between them.
        It is analog of harmonic oscillator potential between them.
    """

    def __init__(self, location, parameters, _method='LSODA'):

        simulation = Simulator(location, parameters, len(location.rcQ), method=_method,
                               mode='fit')
        self.energy = self.energy(location.rcQ, simulation.cQ)

    def energy(self, f, g):
        erg = 0
        for a, b in zip(f, g):  # v007.4
            erg += (a-b)**2
        return(sqrt(erg)/sum(f))  # Compute sum outside of comparator ?


class GeneticFit:
    """
        This class implements an evolutionary genetic algorithm that
        finds the best fitting SEIR parameters for a 'location' (pop + cases).
        At each generation we keep the best 'eliteSize' parameters and breed
        them with a 'mutationRate'.
        We stop at 'maxGen' generations.
        The best result (parameters) is 'self.elite'.
        Set 'verbose' = True to see the progress.
    """

    def __init__(self, location, eliteSize, maxGen, mutationRate=0.4,
                 intervals=ref_intervals, method='LSODA', verbose=False):
        self.location = location
        self.eliteSize = eliteSize
        self.maxGen = maxGen
        self.mutationRate = mutationRate
        self.intervals = intervals

        # Estimate I0
        self.intervals[4][0] = location.rcQ[0]
        self.intervals[4][1] = location.rcQ[0]*5

        # Estimate E0
        self.intervals[5][0] = location.rcQ[0]
        self.intervals[5][1] = location.rcQ[0]*5

        self.method = method
        self.verbose = verbose

        (self.lembda, self.kappa) = self.compute_lembda_kappa()

        self.evolution = self.evolve()
        self.elite = Configurator(self.evolution[0])
        self.energy = self.evolution[1]

    def reorder_paremeters(self, parameters):
        return([parameters[0], parameters[1], parameters[2], parameters[3],
                self.lembda, self.kappa, parameters[4], parameters[5]])

    def compute_lembda_kappa(self):
        """
            Determine Lambda and Kappa from real data.
            Note: 'lembda' is used instead of the Python reserved name lambda.
        """
        lembda = []
        kappa = []
        cQ = np.array(self.location.rcQ)
        R = np.array(self.location.rR)
        D = np.array(self.location.rD)
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
        popSize = self.eliteSize*(self.eliteSize-1)

        width = []
        for i in range(0, len(self.intervals)):
            width.append(self.intervals[i][1]-self.intervals[i][0])

        population = []
        elite = []

        popSize = self.eliteSize*(self.eliteSize-1)

        # Generating the initial population
        for i in range(popSize):
            parameters = []
            for k in range(len(self.intervals)):
                parameters.append(random()*width[k]+self.intervals[k][0])

            # compare the simulated scenario with real one
            config = Configurator(self.reorder_paremeters(parameters))
            comparisonValue = Comparator(self.location, config, self.method).energy

            # Put that into population
            # Fitting the current parameters
            population.append([parameters, comparisonValue])

        # Selecting the elite
        elite = self.generate_elite(population, self.eliteSize)

        for g in range(self.maxGen):
            if self.verbose:
                print(f"Gen {g+1} of {self.maxGen}")
            newPop = deepcopy(elite)
            for i in range(self.eliteSize):
                for j in range(self.eliteSize):
                    if j == i:
                        continue
                    # Breed and mutate
                    child = self.cross(elite[i][0], elite[j][0],
                                       self.mutationRate, width,
                                       self.intervals)
                    config = Configurator(self.reorder_paremeters(child))

                    #  compare the simulated scenario with real one
                    comparisonValue = Comparator(self.location, config,
                                                 self.method).energy
                    newPop.append([child, comparisonValue])

            elite = self.generate_elite(newPop, self.eliteSize)
            if self.verbose:
                print(f'Best: {elite[0][1]}\n')

        return ([self.reorder_paremeters(elite[0][0]), elite[0][1]])

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

    def cross(self, par1, par2, mutationChance, width, intervals):
        """
             Generates a child by breeding 2 parents, with
            a 'mutationChance' = 0-1.
        """
        child = []
        for i in range(len(par1)):
            if random() < mutationChance:
                #   Mutation
                child.append(random()*width[i]+intervals[i][0])
            else:
                #   Breeding
                rand = round(random())
                child.append(rand*par1[i] + (1-rand)*par2[i])
        return child

    def generate_elite(self, population, eliteSize):
        """
            Sorts the population by their scores and
            returns 'eliteSize' best elements.
        """
        elite = []
        # Selecting the elite
        self.sort_population(population)
        for i in range(eliteSize):
            elite.append(population[i])
        return elite


class Analyzer:
    # It is a tool to measure BRN, peakTime, peak
    # for a specific parameters configuration and Region
    """
        This class measures BRN, peakTime, Peak values for a given 'Region'
        and a set of SEIR parameters.
        It gives also the data of simulation for nDays period after
        the first day.
    """
    def __init__(self, location, parameters, nDays):
        self.location = location
        self.parameters = parameters
        self.nDays = nDays

        simulation = Simulator(location, parameters, nDays)
        activeCases = simulation.Q

        peak = max(simulation.Q)
        peakDay = simulation.Q.index(peak)
        self.peak_is_outofrange = False
        if peakDay == len(simulation.Q)-1:
            self.peak_is_outofrange = True

        simulation = Simulator(location, parameters, nDays)

        self.peak = activeCases[peakDay-1]
        self.peakTime = peakDay
        self.BRN = parameters.beta/parameters.delta
        self.data = simulation

    def find_starting_point(self):
        """
        This can be used to seek the probable first day when the
        first exposition occured.
        """
        solver = slv.SEIRSolver(self.parameters, self.location)

        step = -30
        solver.modelSol(0, step, 1, 'LSODA')
        E = solver.E

        while E[0] > 1:
            step += -30
            # Simulates the previous 30 days
            solver.modelSol(step+30, step, 1, 'LSODA')
            E = solver.E
        starting_point = 0

        while E[starting_point] < 0:
            starting_point += 1

        starting_point += step
        starting_date = self.location.first_day + timedelta(days=starting_point)
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
            >> ita_exp = SEIR.BaseExperiment(location=SEIR.LOCATIONS['Italy'],
                                             nDays=100, nExp=4)
            >> ita_exp.run()
            >> ita_exp.compute_results()
            >> ita_exp.show_results('all')
    """

    def __init__(self, location, nDays, nExp,
                 eliteSize=10, maxGen=50, mutationRate=0.4,
                 intervals=ref_intervals, method='LSODA', verbose=False):

        self.location = location
        self.nDays = nDays
        self.nExp = nExp
        self.eliteSize = eliteSize
        self.maxGen = maxGen
        self.mutationRate = mutationRate
        self.intervals = intervals
        self.method = method
        self.verbose = verbose

    def init_states(self):

        # typical values
        self.typical_BRN = 0
        self.typical_experiment_index = 0
        self.typical_parameters = Configurator()
        self.typical_peak = 0
        self.typical_peakTime = 0
        self.typical_peak_found = True

        # observables mean values
        self.BRN = 0
        self.peak = 0
        self.peakTime = 0
        self.peak_found = True

        # observables covariances
        self.errBRN = 0
        self.errPeak = 0
        self.errPeakTime = 0

        # Elite parameters
        self.elites = []

        # All exps. data.
        self.raw_BRN, self.raw_Rn, self.raw_peak, self.raw_peakTime,\
            self.raw_S, self.raw_E, self.raw_I, self.raw_Q, self.raw_cQ,\
            self.raw_R, self.raw_D, self.raw_P\
            = [], [], [], [], [], [], [], [], [], [], [], []

        # mean compartmental variables
        self.S, self.E, self.I, self.Q, self.cQ, self.R, self.D,\
            self.P, self.Rn\
            = [], [], [], [], [], [], [], [], []

        # compartmental variables covariances
        self.errS, self.errE, self.errI, self.errQ, self.errcQ, self.errR,\
            self.errD, self.errP, self.errRn\
            = [], [], [], [], [], [], [], [], []

        # days list of simulation
        self.days = []

    def run(self):
        """
            Runs the simulations WITHOUT computing the final results.
        """
        # Here we perform nExp of experiments
        # (fitting + calculating observables)
        self.init_states()

        for n in range(self.nExp):
            # Fit real data and obtaining the n'th elite parameter
            if self.verbose:
                print("I am fitting real Data ", n+1, "th time of ", self.nExp)
            data = GeneticFit(self.location, self.eliteSize, self.maxGen,
                              self.mutationRate, self.intervals, self.method,
                              self.verbose)

            # Calculate observables of n'th elite parameter
            if self.verbose:
                print("I am calculating observables")
            observables = Analyzer(self.location, data.elite, self.nDays)
            # Assign what you have obtained to lists
            self.elites.append(data.elite)
            self.raw_BRN.append(observables.BRN)
            self.raw_peak.append(observables.peak)
            self.raw_peakTime.append(observables.peakTime)
            self.raw_S.append(observables.data.S)
            self.raw_E.append(observables.data.E)
            self.raw_I.append(observables.data.I)
            self.raw_Q.append(observables.data.Q)
            self.raw_cQ.append(observables.data.cQ)
            self.raw_R.append(observables.data.R)
            self.raw_D.append(observables.data.D)
            self.raw_P.append(observables.data.P)
            # R(t) = BRN * S/N  v007.4
            self.raw_Rn.append([observables.BRN * s/self.location.N
                                for s in self.raw_S[n]])
            if observables.peak_is_outofrange:
                self.peak_found = False

    def compute_results(self):
        """
            Computes the final results of the simulations.
        """
        # Calculate the {mean values} and {typical values}
        # of observables from nExp of experiments
        if self.verbose:
                print("I am calculating expected values")
        self.BRN = np.mean(np.array(self.raw_BRN))

        self.typical_BRN = self.raw_BRN[0]
        for n in range(1, self.nExp):
            if (self.raw_BRN[n]-self.BRN)**2 < (self.typical_BRN-self.BRN)**2:
                self.typical_BRN = self.raw_BRN[n]
                self.typical_experiment_index = n
            else:
                continue

        self.typical_parameters = self.elites[self.typical_experiment_index]

        self.peak = np.mean(np.array(self.raw_peak))
        self.typical_peak = self.raw_peak[self.typical_experiment_index]

        self.peakTime = np.mean(np.array(self.raw_peakTime))
        self.typical_peakTime\
            = self.raw_peakTime[self.typical_experiment_index]

        # Calculate the mean values of compartmental variables
        for i in range(self.nDays):
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
        if self.verbose:
                print("I am calculating uncertainties")
        self.errBRN = sp.sem(np.array(self.raw_BRN)) * 1.96
        self.errPeak = sp.sem(np.array(self.raw_peak)) * 1.96
        self.errPeakTime = sp.sem(np.array(self.raw_peakTime)) * 1.96

        # Calculate the std error on comartmental variables  v007.4
        #  * 1.96 means sigma2 : 95% confidence
        for i in range(self.nDays):
            self.days.append(i)
            self.errS.append(sp.sem(np.array(self.raw_S)[:, i]) * 1.96)
            self.errE.append(sp.sem(np.array(self.raw_E)[:, i]) * 1.96)
            self.errI.append(sp.sem(np.array(self.raw_I)[:, i]) * 1.96)
            self.errQ.append(sp.sem(np.array(self.raw_Q)[:, i]) * 1.96)
            self.errcQ.append(sp.sem(np.array(self.raw_cQ)[:, i]) * 1.96)
            self.errR.append(sp.sem(np.array(self.raw_R)[:, i]) * 1.96)
            self.errD.append(sp.sem(np.array(self.raw_D)[:, i]) * 1.96)
            self.errP.append(sp.sem(np.array(self.raw_P)[:, i]) * 1.96)
            self.errRn.append(sp.sem(np.array(self.raw_Rn)[:, i]) * 1.96)

    def show_results(self, *curves):
        """
            Shows the curves and prints the results of nExp of experiments
        """
        print(f'{Fore.GREEN}'
              f'___________________ Plots ___________________'
              f'{Style.RESET_ALL}')
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
        fig, plots = pl.subplots(nrows=n_plots, ncols=1, sharex=False,
                                 figsize=[8, (13/3)*n_plots])

        if n_plots == 1:
            plots = (plots,)

        last_plot = -1
        if 'cq' in curves or draw_all:
            last_plot += 1
            cq = plots[last_plot]
            cq.set_title('Cumulative Number of cases')
            cq.plot(self.location.rcQ, color='red', marker="o",
                    label='Real data (fitting)')
            cq.errorbar(self.days, self.cQ, yerr=self.errcQ, errorevery=6,
                        ecolor='black', label='Model')
            cq.legend(loc="lower right")
            cq.set_xlabel(f'Days since {self.location.first_day}')
            cq.set_ylabel('cQ(t)')
            cq.grid(True)

        if 'q' in curves or draw_all:
            last_plot += 1
            q = plots[last_plot]
            q.set_title('Number of active cases')
            q.errorbar(self.days, self.Q, yerr=self.errQ, errorevery=6,
                       ecolor='black')
            q.set_xlabel(f'Days since {self.location.first_day}')
            q.set_ylabel('Q(t)')
            q.grid(True)

        if 'rn' in curves or draw_all:
            last_plot += 1
            rn = plots[last_plot]
            rn.set_title('Reproduction Number')
            rn.errorbar(self.days, self.Rn, yerr=self.errRn, errorevery=6,
                        ecolor='black', label='R(t)')
            # Drawing Rn = 1
            rn.plot(self.days, [1]*len(self.Rn), label='R=1', color='red',
                    linestyle='dashed')
            rn.legend(loc="upper right")
            rn.set_xlabel(f'Days since {self.location.first_day}')
            rn.set_ylabel('R(t)')
            rn.grid(True)

        fig.tight_layout()
        pl.show()

    def show_result_values(self):
        """
            Shows the numerical values of the experiment.
        """
        # Showing mean values
        print(f'{Fore.GREEN}'
              f'_______________________ Mean Values _______________________'
              f'{Style.RESET_ALL}')
        print('BRN: ', round(self.BRN, 3), "\u00B1", round(self.errBRN, 3))

        if self.peak_found:
            print('PeakTime: ',
                  self.location.first_day + timedelta(round(self.peakTime)),
                  ' (', round(self.peakTime, 2), ") \u00B1",
                  round(self.errPeakTime, 2))
            print('Peak:', int(self.peak), "\u00B1", int(self.errPeak))
        else:
            print("A peak could not be found in the chosen time period.")
        # Totals (mean)
        print('')
        print(f'Total infected (t={self.nDays}):',
              int(self.cQ[-1]), "\u00B1", int(self.errcQ[-1]))
        print(f'Total recovered (t={self.nDays}):',
              int(self.R[-1]), "\u00B1", int(self.errR[-1]))
        print(f'Total deaths (t={self.nDays}):',
              int(self.D[-1]), "\u00B1", int(self.errD[-1]))

        # Showing typical values
        print(f'{Fore.GREEN}'
              f'_______________________ Typical Values _______________________'
              f'{Style.RESET_ALL}')

        print(f'{Fore.RED}'
              f'Typical values are chosen to give the closest BRN '
              f'to the mean BRN.'
              f'{Style.RESET_ALL}')

        print('')
        print('BRN: ', round(self.typical_BRN, 3))
        if self.peak_found:
            print('PeakTime: ',
                  self.location.first_day + timedelta(round(self.typical_peakTime)),
                  ' (', round(self.typical_peakTime, 2), ')')
            print('Peak:', int(self.typical_peak))
        else:
            print("A peak could not be found in the chosen time period.")

        tei = self.typical_experiment_index
        # Totals (typical)
        print('')
        print(f'Total infected (t={self.nDays}):',
              int(self.raw_cQ[tei][-1]))
        print(f'Total recovered (t={self.nDays}):',
              int(self.raw_R[tei][-1]))
        print(f'Total deaths (t={self.nDays}):',
              int(self.raw_D[tei][-1]))

        print('')
        print('Protection rate: ', round(self.typical_parameters.alpha, 3))
        print('Transmission rate: ', round(self.typical_parameters.beta, 3))
        print('Latent period: ', round(1/self.typical_parameters.gamma, 3))
        print('Infection period: ', round(1/self.typical_parameters.delta, 3))
        print('Recovering rate: ', round(self.typical_parameters.lembda, 3))
        print('Fatality rate: ', round(self.typical_parameters.kappa, 3))
        print('Number of exposed people in (', self.location.first_day, ') = ',
              round(self.typical_parameters.E0))
        print('Number of infective people in (', self.location.first_day, ') = ',
              round(self.typical_parameters.I0))

    def __add__(self, other):
        """
            Defines the + operator for BaseExperiment.
            NOTE : This is just a hack for now!
            To Do: Make this more robust
        """
        newExp = deepcopy(self)
        newExp.elites += other.elites
        newExp.raw_BRN += other.raw_BRN
        newExp.raw_peak += other.raw_peak
        newExp.raw_peakTime += other.raw_peakTime
        newExp.raw_S += other.raw_S
        newExp.raw_E += other.raw_E
        newExp.raw_I += other.raw_I
        newExp.raw_Q += other.raw_Q
        newExp.raw_cQ += other.raw_cQ
        newExp.raw_R += other.raw_R
        newExp.raw_D += other.raw_D
        newExp.raw_P += other.raw_P
        newExp.raw_Rn += other.raw_Rn
        newExp.nExp += other.nExp

        # observables mean values
        newExp.BRN = 0
        newExp.peak = 0
        newExp.peakTime = 0
        newExp.peak_found = newExp.peak_found and other.peak_found

        # observables covariances
        newExp.errBRN = 0
        newExp.errPeak = 0
        newExp.errPeakTime = 0

        # mean compartmental variables
        newExp.S, newExp.E, newExp.I, newExp.Q, newExp.cQ, newExp.R,\
            newExp.D, newExp.P, newExp.Rn = [], [], [], [], [], [], [], [], []

        # compartmental variables covariances
        newExp.errS, newExp.errE, newExp.errI, newExp.errQ, newExp.errcQ,\
            newExp.errR, newExp.errD, newExp.errP, newExp.errRn = [], [], [],\
            [], [], [], [], [], []

        newExp.days = []

        return newExp
