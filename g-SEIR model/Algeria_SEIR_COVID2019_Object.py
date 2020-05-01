"""
Version 007.1 = 004.2 + 006 + additions
*********
This version is object oriented and uses 5 methods via
SEIR_Solver to solve the SEIR model.
It supports multiprocessing

Taha:
    - modif on DATA part
    1- 5/04 modif on line 535
"""

from math import sqrt
from random import random
import matplotlib.pyplot as pl
import numpy as np
from copy import deepcopy #Is needed to create hard copies
from datetime import date
from datetime import timedelta
import pandas as pd

import SEIR_Solver as slv
from DataCollector import Region, WebDataReader

"""Our SEIR parameters are bound to these intervals
They are: [0:alpha, 1:beta, 2:gamma, 3:delta, 4:I0, 5:E0 ]"""
ref_intervals = [[0,0.15],[0.2,1.5],[0.3,0.6],[1/7.5,1/3.5],[0,50],[0,50]] 


################    TAHA - DATA    ################
#  (change input-data file's paths and names)

import os        
# Go to input-data directory 
os.chdir('/Users/rouabah/MacWork/MyResearch/COVID19/Numerics/Data/Inputs') 
sk_h = 42
sk_f = 0

# Algeria
Algeria = Region()
Algeria.name = 'Algeria'
Algeria.N = 43851044  #worldometers.info/population
Algeria.first_day = date(2020,3,1)
Algeria.rcQ = list(np.genfromtxt('Algeria_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skip_header=sk_h, skip_footer=sk_f, usecols=(5))) 
Algeria.rR = list(np.genfromtxt('Algeria_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skip_header=sk_h, skip_footer=sk_f, usecols=(6)))
Algeria.rD = list(np.genfromtxt('Algeria_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skip_header=sk_h, skip_footer=sk_f, usecols=(7))) 

# Algeria Center  (Algiers+Blida)
AlgeriaCenter = Region()
AlgeriaCenter.name = 'AlgeriaCenter'
AlgeriaCenter.N = 2200000  #worldometers.info/population # Algiers: 1,977,663, Blida: 182,447
AlgeriaCenter.first_day = date(2020,3,12)
AlgeriaCenter.rcQ = list(np.genfromtxt('Algeria_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skip_header=53, skip_footer=0, usecols=(6))*(0.45)) 
AlgeriaCenter.rR = list(np.genfromtxt('Algeria_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skip_header=53, skip_footer=0, usecols=(6))*(0.46))     # 1(Confirmed)-0.53(death)-0.01(active)
AlgeriaCenter.rD =list(np.genfromtxt('Algeria_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skip_header=53, skip_footer=0, usecols=(6))*(0.53)) 

# France
France = Region()
France.name = 'France'
France.N = 65273511  #worldometers.info/population
France.first_day = date(2020,2,15) # Since non-zero cQ,R & D.
France.rcQ = list(np.loadtxt('France_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=25, usecols=(5)))
France.rR = list(np.loadtxt('France_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=25, usecols=(6)))
France.rD = list(np.loadtxt('France_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=25, usecols=(7)))

# Germany
Germany = Region()
Germany.name = 'Germany'
Germany.N = 83783942  #worldometers.info/population
Germany.first_day = date(2020,3,9) # Since non-zero cQ,R & D.
Germany.rcQ = list(np.loadtxt('Germany_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=48, usecols=(5)))
Germany.rR = list(np.loadtxt('Germany_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=48, usecols=(6)))
Germany.rD = list(np.loadtxt('Germany_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=48, usecols=(7)))

# Hubei
Hubei = Region()
Hubei.name = 'Hubei'
Hubei.N = 59170000  #statista.com
Hubei.first_day = date(2020,1,22) # All non-zero cQ,R & D.
Hubei.rcQ = list(np.loadtxt('Hubei_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=1, usecols=(5)))
Hubei.rR = list(np.loadtxt('Hubei_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=1, usecols=(6)))
Hubei.rD = list(np.loadtxt('Hubei_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=1, usecols=(7)))

# Italy
Italy = Region()
Italy.name = 'Italy'
Italy.N = 60483054  #Worldometer.info
Italy.first_day = date(2020,2,22) # Since non-zero cQ,R & D.
Italy.rcQ = list(np.loadtxt('Italy_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=34, usecols=(5)))[:35]
Italy.rR = list(np.loadtxt('Italy_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=34, usecols=(6)))[:35]
Italy.rD = list(np.loadtxt('Italy_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=34, usecols=(7)))[:35]

# Spain
Spain = Region()
Spain.name = 'Spain'
Spain.N = 46754778  #worldometers.info/population
Spain.first_day = date(2020,3,3) # Since non-zero cQ,R & D.
Spain.rcQ = list(np.loadtxt('Spain_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=42, usecols=(5)))
Spain.rR = list(np.loadtxt('Spain_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=42, usecols=(6)))
Spain.rD =list(np.loadtxt('Spain_since22012020.txt',dtype='int',delimiter=',', comments='#' ,skiprows=42, usecols=(7)))

# PLEASE update this dictionary when you add new regions
Locations = {'Italy':Italy, 'Hubei':Hubei, 'Algeria':Algeria, 'AlgeriaCenter':AlgeriaCenter, 'Germany':Germany, 'Spain':Spain, 'France':France}


# ################    DATA    ################

# # PLEASE update this dictionary when you add new regions
# Locations = {}

# # Download data from web source
# dataset = pd.read_csv('time-series-19-covid-combined.csv')
# print("Dataset is loaded!")

# # Austria
# Austria = WebDataReader(dataset, 'Austria').Region
# Austria.N = 8994481  #Worldometer.info
# Locations[Austria.name] = Austria

# # Italy
# Italy = WebDataReader(dataset, 'Italy').Region
# Italy.N = 60483054  #Worldometer.info
# Locations[Italy.name] = Italy

# # Hubei
# Hubei = WebDataReader(dataset, 'China', 'Hubei').Region
# Hubei.N = 59170000  #statista.com
# Locations[Hubei.name] = Hubei

# # Algeria # BE CAREFUL WITH THE DATA FROM time-series-19-covid-combined.csv !!!
# #Algeria = WebDataReader(dataset, 'Algeria').Region
# #Algeria.N = 43665440 #Worldometer.info
# #Locations[Algeria.name] = Algeria

# # Spain
# Spain =  WebDataReader(dataset, 'Spain').Region
# Spain.N = 46754778  #worldometers.info/population
# Locations[Italy.name] = Spain


# print(list(Locations.keys()))


################    DATA    ################

class Configurator:
    """
        Maps a list of SEIR parameters to a human-readable object
    """
    def __init__(self,p=[0,0,0,0,0,0,0,0]):
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
        return("{%s,%s,%s,%s,%s,%s,%s,%s}" % (self.alpha, self.beta,self.gamma,self.delta,self.lembda,self.kappa,self.I0,self.E0))
    
    def __str__(self):
        return("{%s,%s,%s,%s,%s,%s,%s,%s}" % (self.alpha, self.beta,self.gamma,self.delta,self.lembda,self.kappa,self.I0,self.E0))
    
    def values(self):
        return([self.alpha, self.beta, self.gamma, self.delta, self.lembda, self.kappa , self.I0, self.E0])


class Simulator:
    """
        Solves the SEIR model in a 'city' for 'parameters' during 'nDays'.
        In this version we use the SEIRSolver tool.
    """
    def __init__(self, city, parameters, nDays, method = 'LSODA', mode = 'sim'): 
        """
            method = 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA' // DOP??
        """
        solver = slv.SEIRSolver(parameters, city, mode)
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
    def __init__(self,city,parameters, _method):
        
        simulation = Simulator(city, parameters, len(city.rcQ), method = _method, mode = 'fit')
        self.energy = self.energy(city.rcQ,simulation.cQ)
    
            
    def energy(self, f, g):
        
        erg = 0;
        for i in range(0,len(f)):
            erg = erg + (f[i]-g[i])**2
        return(sqrt(erg)/sum(f))

class GeneticFit:
    """
        This class implements an evolutionary genetic algorithm that
        finds the best fitting SEIR parameters for a 'city' (pop + cases).
        At each generation we keep the best 'eliteSize' parameters and breed
        them with a 'mutationRate'.
        We stop at 'maxGen' generations.
        The best result (parameters) is 'self.elite'.
        Set 'verbose' = True to see the progress.
    """         
    
    def __init__(self, city, eliteSize, maxGen, mutationRate = 0.4, intervals = ref_intervals, method = 'LSODA', verbose = False):
        self.city = city
        self.eliteSize = eliteSize
        self.maxGen = maxGen
        self.mutationRate = mutationRate
        self.intervals = intervals
        
        # Estimate I0
        self.intervals[4][0] = city.rcQ[0]
        self.intervals[4][1] = city.rcQ[0]*5
        
        # Estimate E0
        self.intervals[5][0] = city.rcQ[0]
        self.intervals[5][1] = city.rcQ[0]*5
        
        self.method = method
        self.verbose = verbose
        
        (self.lembda, self.kappa) = self.compute_lembda_kappa()
        
        # ??
        self.evolution = self.evolve()
        self.elite = Configurator(self.evolution[0])
        self.energy = self.evolution[1]
        
    def reorder_paremeters(self, parameters):
        return([parameters[0], parameters[1], parameters[2], parameters[3], self.lembda, self.kappa, parameters[4],parameters[5]])
        
    def compute_lembda_kappa(self):
        """
            Determine Lembda and Kappa from real data
        """
        lembda = []
        kappa = []
        cQ = np.array(self.city.rcQ)
        R = np.array(self.city.rR)
        D = np.array(self.city.rD)
        Q = cQ - R - D
        
        l = len(Q)-1
        for n in range(l):
            lembda.append((R[n+1]-R[n])/Q[n])
            kappa.append((D[n+1]-D[n])/Q[n])
        #self.l = lembda
        #self.k = kappa
        Lembda = np.mean(np.array(lembda))
        Kappa = np.mean(np.array(kappa))
        return(Lembda,Kappa)

        
    def evolve(self):
        popSize = self.eliteSize*(self.eliteSize-1)
        
        
        width = []
        for i in range(0,len(self.intervals)):
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
            comparisonValue = Comparator(self.city, config, self.method).energy
            
            # Put that into population
            population.append([parameters,comparisonValue]) # fitting the current parameters
            
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
                    child = self.cross(elite[i][0], elite[j][0], self.mutationRate, width, self.intervals)
                    config = Configurator(self.reorder_paremeters(child))
                    
                    #  compare the simulated scenario with real one
                    comparisonValue = Comparator(self.city, config, self.method).energy
                    newPop.append([child,comparisonValue])
                    
            elite = self.generate_elite(newPop, self.eliteSize)
            if self.verbose:
               print(f'Best: {elite[0][1]}\n')
              
        return([self.reorder_paremeters(elite[0][0]),elite[0][1]] )
            
    def sort_population(self,sub_li): 
        """
            Sorts the population of the genetic pool by
            their fitting score (2nd element in population list)
        """
        # reverse = None (Sorts in Ascending order) 
        sub_li.sort(key = lambda x: x[1]) 
        return sub_li 
    
    
    def cross(self, par1, par2, mutationChance, width, intervals):
        """
             Generates a child by breeding 2 parents, with
            a 'mutationChance' = 0-1.
        """
        child = []
        for i in range(len(par1)):
            if random()<mutationChance:
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
            
# It is a tool to measure BRN, peakTime, peak 
# for a specific parameters configuration and Region
class Analyzer:
    """
        This class measures BRN, peakTime, Peak values for a given 'Region' and a set of SEIR parameters
        It gives also the data of simulation for nDays period after the first day
    """
    def __init__(self, city, parameters, nDays):
        self.city = city
        self.parameters = parameters
        self.nDays = nDays
        
        simulation = Simulator(city, parameters, nDays)
        activeCases = simulation.Q
        
        peak = max(simulation.Q)
        peakDay = simulation.Q.index(peak)
        self.peak_is_outofrange = False
        if peakDay == len(simulation.Q)-1:
            self.peak_is_outofrange = True
        
        simulation = Simulator(city,parameters,nDays)
        
        self.peak = activeCases[peakDay -1]
        self.peakTime = peakDay
        self.BRN = parameters.beta/parameters.delta
        self.data = simulation
    
    def find_starting_point(self):
        """
        This can be used to seek the probable first day that first exposition is occured    
        """
        solver = slv.SEIRSolver(self.parameters, self.city)
        
        step = -30;
        solver.modelSol(0, step, 1, 'LSODA')
        E = solver.E

        while E[0]>1:
            step += -30;
            solver.modelSol(step+30, step, 1, 'LSODA') #simulates the previous 30 days
            E = solver.E
        starting_point = 0

        while E[starting_point]<0 :
            starting_point += 1
        
        starting_point += step
        starting_date = self.city.first_day + timedelta(days = starting_point)
        return(starting_date)        
        
class BaseExperiment:
    """
        This is class perform nExp number of experiments
        (genetic fitting + simulating + measuring observables)
    """
    def __init__(self, city, nDays, nExp,
                 eliteSize = 10, maxGen = 500, mutationRate = 0.4,
                 intervals = ref_intervals, method = 'LSODA', verbose = False):
        
        self.city = city
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
        self.raw_BRN, self.raw_Rn, self.raw_peak, self.raw_peakTime, self.raw_S,\
        self.raw_E, self.raw_I, self.raw_Q, self.raw_cQ, \
        self.raw_R, self.raw_D, self.raw_P = [],[],[],[],[],[],[],[],[],[],[],[]
        
        # mean compartmental variables
        self.S, self.E, self.I, self.Q, self.cQ, self.R, self.D, self.P, self.Rn =\
        [],[],[],[],[],[],[],[], []
        
        # compartmental variables covariances 
        self.errS, self.errE, self.errI, self.errQ, self.errcQ, self.errR,\
        self.errD, self.errP, self.errRn = [],[],[],[],[],[],[],[],[]
        
        # days list of simulation
        self.days = [] 
    
    def run(self):        
        # Here we perform nExp of experiments (fitting + calculating observables)
        self.init_states()
        
        for n in range(self.nExp):            
            # Fit real data and obtaining the n'th elite parameter 
            if self.verbose:
                print("I am fitting real Data ", n+1,"th time of ", self.nExp)
            data = GeneticFit(self.city, self.eliteSize, self.maxGen, self.mutationRate, self.intervals, self.method, self.verbose);
            
            # Calculate observables of n'th elite parameter
            if self.verbose:
                print("I am calculating observables")
            observables = Analyzer(self.city, data.elite, self.nDays)
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
            #self.raw_Rn.append(observables.data.Rn)
            self.raw_Rn.append(self.raw_BRN[n]*(1-self.elites[n].alpha)**np.linspace(0,self.nDays,self.nDays+1))
            if observables.peak_is_outofrange:
                self.peak_found = False
    
    def compute_results(self):
        # Calculate the {mean values} and {typical values} of observables from nExp of experiments
        if self.verbose:
                print("I am calculating expected values")
        self.BRN = np.mean(np.array(self.raw_BRN))
        
        self.typical_BRN = self.raw_BRN[0]
        for n in range(1,self.nExp):
            if (self.raw_BRN[n]-self.BRN)**2 < (self.typical_BRN-self.BRN)**2:
                self.typical_BRN = self.raw_BRN[n]
                self.typical_experiment_index = n
            else:
                continue
        
        self.typical_parameters = self.elites[self.typical_experiment_index]
        
        self.peak = np.mean(np.array(self.raw_peak))
        self.typical_peak = self.raw_peak[self.typical_experiment_index]
        
        self.peakTime = np.mean(np.array(self.raw_peakTime))
        self.typical_peakTime = self.raw_peakTime[self.typical_experiment_index]
        
        # Calculate the mean values of compartmental variables
        for i in range(self.nDays):
            self.S.append(np.mean(np.array(self.raw_S)[:,i]))
            self.E.append(np.mean(np.array(self.raw_E)[:,i]))
            self.I.append(np.mean(np.array(self.raw_I)[:,i]))
            self.Q.append(np.mean(np.array(self.raw_Q)[:,i]))
            self.cQ.append(np.mean(np.array(self.raw_cQ)[:,i]))
            self.R.append(np.mean(np.array(self.raw_R)[:,i]))
            self.D.append(np.mean(np.array(self.raw_D)[:,i]))
            self.P.append(np.mean(np.array(self.raw_P)[:,i]))
            self.Rn.append(np.mean(np.array(self.raw_Rn)[:,i]))
            
        
        # # Calculate the covariance of observables
        # if self.verbose:
        #         print("I am calculating uncertainties")  
        # self.errBRN = np.std(np.array(self.raw_BRN))
        # self.errPeak = np.std(np.array(self.raw_peak))
        # self.errPeakTime = np.std(np.array(self.raw_peakTime))
        
        # # Calculate the covariance of comartmental variables
        # for i in range(self.nDays):
        #     self.days.append(i)
        #     self.errS.append(np.std(np.array(self.raw_S)[:,i]))
        #     self.errE.append(np.std(np.array(self.raw_E)[:,i]))
        #     self.errI.append(np.std(np.array(self.raw_I)[:,i]))
        #     self.errQ.append(np.std(np.array(self.raw_Q)[:,i]))
        #     self.errcQ.append(np.std(np.array(self.raw_cQ)[:,i]))
        #     self.errR.append(np.std(np.array(self.raw_R)[:,i]))
        #     self.errD.append(np.std(np.array(self.raw_D)[:,i]))
        #     self.errP.append(np.std(np.array(self.raw_P)[:,i]))
        #     self.errRn.append(np.std(np.array(self.raw_Rn)[:,i]))
            
        # Calculate the covariance of observables
        if self.verbose:
                print("I am calculating uncertainties")  
        self.errBRN = np.std(np.array(self.raw_BRN)) /self.nExp
        self.errPeak = np.std(np.array(self.raw_peak)) /self.nExp
        self.errPeakTime = np.std(np.array(self.raw_peakTime)) /self.nExp
        
        # Calculate the covariance of comartmental variables
        for i in range(self.nDays):
            self.days.append(i)
            self.errS.append(np.std(np.array(self.raw_S)[:,i]) /self.nExp)
            self.errE.append(np.std(np.array(self.raw_E)[:,i]) /self.nExp)
            self.errI.append(np.std(np.array(self.raw_I)[:,i]) /self.nExp)
            self.errQ.append(np.std(np.array(self.raw_Q)[:,i]) /self.nExp)
            self.errcQ.append(np.std(np.array(self.raw_cQ)[:,i]) /self.nExp)
            self.errR.append(np.std(np.array(self.raw_R)[:,i]) /self.nExp)
            self.errD.append(np.std(np.array(self.raw_D)[:,i]) /self.nExp)
            self.errP.append(np.std(np.array(self.raw_P)[:,i]) /self.nExp)
            self.errRn.append(np.std(np.array(self.raw_Rn)[:,i]) /self.nExp)
        
        
    def show_results(self):
        from colorama import Fore, Style
        #To do: add Rn(t) TODAY + Daily new cases as default results
        #       show_results should become dynamic
        """
            Show the plots and print the results of nExp of experiments
        """
        print(f'{Fore.GREEN}_______________________ Plots _______________________{Style.RESET_ALL}')
        # Plot cQ and Q with error bars
        fig, (cq,q,rn) = pl.subplots(nrows = 3, ncols=1, sharex=False, figsize=[6,10])


        cq.set_title('Cumulative Number of cases')
        cq.plot(self.city.rcQ, color = 'red', label='Real data')
        cq.errorbar( self.days, self.cQ, yerr=self.errcQ, errorevery=6, label = 'Prediction')
        cq.legend(loc="lower right")
        cq.set_xlabel(f'Days since {self.city.first_day}')
        cq.set_ylabel('cQ(t)')
        cq.grid(True)
        
        q.set_title('Number of active cases')
        q.errorbar( self.days, self.Q, yerr=self.errQ, errorevery=6)
        q.set_xlabel(f'Days since {self.city.first_day}')
        q.set_ylabel('Q(t)')
        q.grid(True)
        
        rn.set_title('Reproduction Number')
        rn.errorbar( self.days, self.Rn, yerr=self.errRn, errorevery=6, label='R(t)')
        #Drawing Rn = 1
        rn.plot(self.days, [1]*len(self.Rn), label='R=1', color='red', linestyle='dashed')
        rn.legend(loc="upper right")
        rn.set_xlabel(f'Days since {self.city.first_day}')
        rn.set_ylabel('R(t)')
        rn.grid(True)
        
        fig.tight_layout()
        pl.show()
        
        # Showing mean values
        print(f'{Fore.GREEN}_______________________ Mean Values _______________________{Style.RESET_ALL}')
        print('BRN: ', round(self.BRN,3),"\u00B1",round(self.errBRN,3));
        if self.peak_found:
            print('PeakTime: ', self.city.first_day + timedelta(round(self.peakTime)),\
                  ' (', round(self.peakTime,2),") \u00B1",round(self.errPeakTime,2));
            print('Peak:', round(self.peak),"\u00B1",round(self.errPeak))
        else:
            print("A peak could not be found in the chosen time period.")
            
        # Showing typical values   
        print(f'{Fore.GREEN}_______________________ Typical Values _______________________{Style.RESET_ALL}')
        print(f'{Fore.RED} Typical values are chosen to give the closet value to the mean BRN{Style.RESET_ALL}' )
        print('')
        print('BRN: ', round(self.typical_BRN,3));
        if self.peak_found:
            print('PeakTime: ', self.city.first_day + timedelta(round(self.typical_peakTime)),\
                  ' (', round(self.typical_peakTime,2),')');
            print('Peak:', round(self.typical_peak));
        else:
            print("A peak could not be found in the chosen time period.")
        print('Protection rate: ', round(self.typical_parameters.alpha,3))
        print('Transmission rate: ', round(self.typical_parameters.beta,3))
        print('Latent period: ', round(1/self.typical_parameters.gamma,3))
        print('Infection period: ', round(1/self.typical_parameters.delta,3))
        print('Recovering rate: ', round(self.typical_parameters.lembda,3))
        print('Fatality rate: ', round(self.typical_parameters.lembda,3))
        print('Number of exposed people in (',self.city.first_day,') = ', round(self.typical_parameters.E0))
        print('Number of infective people in (',self.city.first_day,') = ', round(self.typical_parameters.I0))
        
        
    def __add__(self, other):
        """
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
        newExp.D, newExp.P, newExp.Rn = [],[],[],[],[],[],[],[], []
        
        # compartmental variables covariances 
        newExp.errS, newExp.errE, newExp.errI, newExp.errQ, newExp.errcQ,\
        newExp.errR, newExp.errD, newExp.errP, newExp.errRn = [],[],[],[],[],[],[],[],[]
        
        newExp.days = []
        
        return newExp
