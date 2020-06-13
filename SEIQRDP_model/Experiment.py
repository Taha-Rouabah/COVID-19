"""
Version 0.1.2
Created on Wed Apr  1 23:27:02 2020

This version works only with Algeria_SEIR_COVID2019_Object v007.4 and above.

A multiprocessed module for Windows/Linux/Mac.

"""
__version__ = '0.1.2'


from time import time

# Multiprocessing packages
import pickle
import os
import multiprocessing
import psutil

# SEIR Model
import SEIQRDP_model.Algeria_SEIR_COVID2019_Object as SEIR


def SEIR_Worker(ID, location, nDays, nExp, _maxGen, _n_processes):
    """
        This worker function creates a BaseExperiment object and
        runs the simulations.
    """
    experi = SEIR.BaseExperiment(location, nDays, nExp, maxGen=_maxGen,
                                 method='LSODA', verbose=False)
    experi.run()
    # Saving the results in a temporary file.
    with open(f"temp_exp_{ID}.data", "wb") as filehandle:
        pickle.dump(experi, filehandle)

    return


class Experiment():
    """
    This is the multithreaded Experiment class.
    """

    def __init__(self, location, nDays, useMP=True):
        # Initializing a void experiment. Just to reserve the object name.
        self.ex = None
        # Do we use or not multiprocessing?
        self.useMultiprocessing = useMP
        SEIR.load_data()
        self.location = SEIR.LOCATIONS[location]
        self.nDays = nDays
        self.nExp = 0
        self.maxGen = 0

        return

    def run(self, nExp, maxGen):
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
        self.nExp = nExp
        self.maxGen = maxGen
        # Setting the number of precesses
        n_processes = 1
        if self.useMultiprocessing:
            # Getting the number of threads available
            # True: counts the physical + logical threads,
            # False: counts only the physical ones
            n_processes = psutil.cpu_count(logical=True)

        # Handeling the surplus of resources (more CPUs than experiments)
        if n_processes > self.nExp:
            n_processes = self.nExp

        # All the processes are managed here.
        # Info
        print(f"Starting {self.nExp} EXPs on {n_processes} cores ...\n")

        # For time measurement
        t0 = time()

        # Initiating the list of processes/threads
        processes = []

        # Computing the left over experiments after dividing
        # the work equally on n_processes
        undividable_count = nExp % n_processes

        # Filling the last list with 'n_processes' processes.
        for i in range(n_processes):
            # Splitting up the work in the most equal manner
            # among the processes
            procNExp = nExp//n_processes
            if undividable_count > 0:
                procNExp += 1
                undividable_count -= 1

            p = multiprocessing.Process(target=SEIR_Worker,
                                        args=[i, self.location,
                                              self.nDays, procNExp,
                                              self.maxGen, n_processes])
            p.start()
            processes.append(p)

        # Starting the processes
        for process in processes:
            process.join()

        print(f"Calculation time = {round(time()-t0, 1)}s \n")

        # Getting the returned results from each process
        with open("temp_exp_0.data", "rb") as filehandle:
            self.ex = pickle.load(filehandle)
        os.remove(f"temp_exp_0.data")
        for ID in range(1, n_processes):
            with open(f"temp_exp_{ID}.data", "rb") as filehandle:
                self.ex += pickle.load(filehandle)
            os.remove(f"temp_exp_{ID}.data")

        # The results of the global experiment are processed
        self.ex.compute_results()
        print(f'{self.nExp} experiments done!')

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
