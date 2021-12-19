"""
Version 0.1.3
Created on Wed Apr  1 23:27:02 2020
This version works only with Algeria_SEIR_COVID2019_Object v008.0 and above.
A multiprocessed module for Windows/Linux/Mac.


BUGS:
    In 'with open("temp_exp_0.data", "rb") as filehandle:'
        If the simulation wasn't carried, this will not find the file and
        throw an error. Must handle this case.

In this version (0.1.3):
    - We can now cutoff the training data to only a limited number of days.
    - Works with seiqrdp_suite 8.0 ONLY.
    - Renamed variables.
"""
__version__ = '0.1.3'


from time import time

# Multiprocessing packages
import pickle
import os
import multiprocessing
import psutil

# SEIQRDP model
import seiqrdp_model.seiqrdp_suite as seiqrdp


def SEIR_Worker(ID, region, n_sim_days, n_experiments, _max_gen, _n_processes):
    """
        This worker function creates a BaseExperiment object and
        runs the simulations.
    """
    experi = seiqrdp.BaseExperiment(region, n_sim_days, n_experiments,
                                    max_gen=_max_gen, method='LSODA',
                                    verbose=False)
    experi.run()
    # Saving the results in a temporary file.
    with open(f"temp_exp_{ID}.data", "wb") as filehandle:
        pickle.dump(experi, filehandle)

    return


class Experiment():
    """
    This is the multithreaded Experiment class.
    """

    def __init__(self, region, population_size,
                 n_sim_days, train_data_size=-1, use_mp=True):
        # Initializing a void experiment. Just to reserve the object name.
        self.ex = None
        # Do we use multiprocessing
        self.use_multiprocessing = use_mp

        self.region = seiqrdp.load_data(region_name='Italy',
                                        pop_size=population_size)

        if train_data_size > 0:
            self.region.cutoff_data(train_data_size)
            print(f'Cutting off the data to the first {train_data_size} days.')
            print(f'Loaded dataset for {self.region.name} '
                  f'from {self.region.first_day} to {self.region.last_day}.')

        self.n_sim_days = n_sim_days
        self.n_experiments = 0
        self.max_gen = 0

        return

    def run(self, n_experiments, max_gen):
        """
        NOTE:
        *****
        Calling MultiExperiment.run() on Spyder/Windows seems OK.
        If there's a problem, try
            if __name__ == '__main__':
                Experiment.run()
        instead.
        """
        # Params
        self.n_experiments = n_experiments
        self.max_gen = max_gen
        # Setting the number of precesses
        n_processes = 1
        if self.use_multiprocessing:
            # Getting the number of threads available
            # True: counts the physical + logical threads,
            # False: counts only the physical ones
            n_processes = psutil.cpu_count(logical=True)

        # Handeling the surplus of resources (more CPUs than experiments)
        if n_processes > self.n_experiments:
            n_processes = self.n_experiments

        # All the processes are managed here.
        # Info
        print(f'Starting {self.n_experiments} EXPs on'
              f' {n_processes} cores ...\n')

        # For time measurement
        t0 = time()

        # Initiating the list of processes/threads
        processes = []

        # Computing the left over experiments after dividing
        # the work equally on n_processes
        undividable_count = n_experiments % n_processes

        # Filling the last list with 'n_processes' processes.
        for i in range(n_processes):
            # Splitting up the work in the most equal manner
            # among the processes
            exp_by_process = n_experiments//n_processes
            if undividable_count > 0:
                exp_by_process += 1
                undividable_count -= 1

            p = multiprocessing.Process(target=SEIR_Worker,
                                        args=[i, self.region,
                                              self.n_sim_days, exp_by_process,
                                              self.max_gen, n_processes])
            p.start()
            processes.append(p)

        # Starting the processes
        for process in processes:
            process.join()

        print(f"Calculation time = {round(time()-t0, 1)}s \n")

        # Getting the returned results from each process
        with open("temp_exp_0.data", "rb") as filehandle:
            self.ex = pickle.load(filehandle)
        os.remove("temp_exp_0.data")
        for ID in range(1, n_processes):
            with open(f"temp_exp_{ID}.data", "rb") as filehandle:
                self.ex += pickle.load(filehandle)
            os.remove(f"temp_exp_{ID}.data")

        # The results of the global experiment are processed
        self.ex.compute_results()
        print(f'{self.n_experiments} experiments done!')

    def show_results(self, *curves):
        if self.ex is None:
            print('Warning (class:Experiment):'
                  ' No experiment has been run yet!')
            return
        self.ex.show_results(*curves)
        return

    def show_curves(self, *curves):  # v0.1.1
        if self.ex is None:
            print('Warning (class:Experiment):'
                  ' No experiment has been run yet!')
            return
        self.ex.show_curves(*curves)
        return

    def show_result_values(self):  # v0.1.1
        if self.ex is None:
            print('Warning (class:Experiment):'
                  ' No experiment has been run yet!')
            return
        self.ex.show_result_values()
        return

    def get_BaseExperiment(self):
        """
            Get the equivalent non multiprocessed BaseExperiment object.
        """
        return self.ex
