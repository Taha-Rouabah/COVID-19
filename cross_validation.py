"""
This script helps finding the best maxGen range to avoid under/over-fitting.
See paper (section 2.2: the overfitting problem)

Training set = Total data MINUS validation data
Validation set = last 'VALIDATION_PERIOD' days of data.
"""
__version__ = '0.1'

from seiqrdp_model.Algeria_SEIR_COVID2019_Object import Comparator, BaseExperiment
import seiqrdp_model.Algeria_SEIR_COVID2019_Object as SEIR
import matplotlib.pyplot as plt
import random
from numpy import sqrt
import time

# Constants
REGION = 'Italy'
VN_RATIO = 1/4  # v/n, where v: validation set size, n: total data set size.
VALIDATION_PERIOD = 0  # in days (will be computed later).
INIT_MAXGEN = 5  # from ...
MAX_MAXGEN = 30  # to
STEP_MAXGEN = 5  # with a step of ...
SEED = 2020  # pick anything you like.

# We need a Comparator which measures the fitness for
# the validation set. Thus we introduce SpecialComparator.


class SpecialComparator(Comparator):
    """
        A special Comparator which fits only for the validation
        period.
        Inherits from Comparator.
    """

    def energy(self, f, g):  # We over-write the original energy function.
        """ Fitness function """
        erg = 0
        for a, b in zip(f[-VALIDATION_PERIOD:],
                        g[-VALIDATION_PERIOD:]):
            erg += (a-b)**2
        return(sqrt(erg)/sum(f))


# Variables
max_gens = []
fits = []

# Loading data and setting the training/validation set...
SEIR.load_data()
tot_region = SEIR.LOCATIONS[REGION]
# Resizing the training set
data_size = len(SEIR.LOCATIONS[REGION].rcQ)
VALIDATION_PERIOD = round(data_size * VN_RATIO)  # Computing the size of 'v'
SEIR.LOCATIONS[REGION].rcQ = SEIR.LOCATIONS[REGION].rcQ[: -VALIDATION_PERIOD]
SEIR.LOCATIONS[REGION].rR = SEIR.LOCATIONS[REGION].rR[: -VALIDATION_PERIOD]
SEIR.LOCATIONS[REGION].rD = SEIR.LOCATIONS[REGION].rD[: -VALIDATION_PERIOD]

print('Data size and validation set size:\n'
      'n =', data_size, '   v =', VALIDATION_PERIOD)

# Setting the timer.
init_time = time.time()

# Evaluating the fitness in the validation set for each maxGen.
for cur_maxgen in range(INIT_MAXGEN, MAX_MAXGEN+STEP_MAXGEN,
                        STEP_MAXGEN):
    # Resetting the random numbers generator.
    random.seed(SEED)

    # Simulating with the given parameters...
    ita_exp = BaseExperiment(SEIR.LOCATIONS[REGION], nDays=data_size, nExp=1,
                             maxGen=cur_maxgen)
    ita_exp.run()

    # Obtaining the fitness.
    fit_comp = SpecialComparator(tot_region, ita_exp.elites[0])
    print('maxGen:', cur_maxgen, '\nFitness:', fit_comp.energy, '\n')
    # Saving data.
    max_gens.append(cur_maxgen)
    fits.append(fit_comp.energy)

print('Done in', int(time.time() - init_time), 'seconds.\n')
# Showing results
plt.plot(max_gens, fits, label=f'Validation period = {VALIDATION_PERIOD}')
plt.legend()
plt.show()
