"""
v002
"""
__version__ = '2.0.1'

import numpy as np
from scipy.integrate import solve_ivp


class SEIRSolver:
    """
        SEIRSolver is a Solver class for our SEIR+Q model.
        It must be initialized with the SEIR 'parameters' (model parameters +
        initial cond.) and the 'location' for which we want to solve the model.
        It uses 2 modes: 'sim' for simulation and 'fit' for fitting only.
        'fit' is faster but only solves for S, E, I, cQ
    """

    def __init__(self, parameters, location, mode='sim'):
        """
            mode = 'sim', 'fit'
        """
        if mode != 'sim' and mode != 'fit':
            print("ERROR (SEIRSolver): Incorrect mode!"
                  "Must be 'sim' or 'fit'.")
        self.mode = mode
        self.parameters = parameters
        self.location = location

        # Initial conditions
        self.I = [self.parameters.I0]  # infections
        self.E = [self.parameters.E0]  # transmitions
        self.Q = [self.location.rcQ[0]]  # Quarentines
        self.cQ = [self.Q[0]]  # cumulative quarentines
        self.S = [self.location.N-self.I[0]-self.E[0]]  # susceptibles
        self.R = [0]  # Recovred
        self.D = [0]  # Deaths
        self.P = [0]  # insusceptibles (by protective measures)
        self.cQ = [self.location.rcQ[0]]
        # Empty results
        self.res = []

    def modelEqts(self, t, x):
        """
            ODEs as taken from Peng et al. 18-feb-2020
            SEIR Model + Q
            x (elements) = 0:I 1:E 2:Q 3:S 4:R 5:D 6:P 7:cQ
            This is our chosen SEIR model.
        """
        # The SEIR model
        dSdt = - (self.parameters.beta*x[3]*x[0])/self.location.N\
            - self.parameters.alpha * x[3]
        dEdt = (self.parameters.beta*x[3]*x[0])/self.location.N\
            - self.parameters.gamma * x[1]
        dIdt = self.parameters.gamma * x[1] - self.parameters.delta * x[0]
        dcQdt = self.parameters.delta * x[0]
        # Simulation mode (complete model)
        if self.mode == 'sim':
            dQdt = self.parameters.delta * x[0]\
                - self.parameters.lembda * x[2]\
                - self.parameters.kappa * x[2]
            dRdt = self.parameters.lembda * x[2]
            dDdt = self.parameters.kappa * x[2]
            dPdt = self.parameters.alpha * x[3]
            return [dIdt, dEdt, dQdt, dSdt, dRdt, dDdt, dPdt, dcQdt]
        # Fitting mode (minimal model)
        elif self.mode == 'fit':
            return [dIdt, dEdt, 0, dSdt, 0, 0, 0, dcQdt]

    def modelSol(self, t_i, t_f, npointsperday=1, method='LSODA'):
        """
            Solves the SEIR model from day 't_i' to day 't_f'.
            Generates 'npointsperday' per day.
            'method' = 'RK45', 'RK23', 'Radau', 'BDF', 'LSODA'
        """
        t_eval = np.linspace(t_i, t_f, (t_f-t_i)*npointsperday + 1)

        # Putting the parameters in a list. Needed by solve_ivp()
        initialCond = [self.I[0], self.E[0], self.Q[0], self.S[0],
                       self.R[0], self.D[0], self.P[0], self.cQ[0]]

        # ODE solver
        self.res = solve_ivp(self.modelEqts,
                             (t_i, t_f), initialCond, method, t_eval)

        # Updating the model's states
        self.I = list(self.res.y[0, :])
        self.E = list(self.res.y[1, :])
        self.Q = list(self.res.y[2, :])
        self.S = list(self.res.y[3, :])
        self.R = list(self.res.y[4, :])
        self.D = list(self.res.y[5, :])
        self.P = list(self.res.y[6, :])
        self.cQ = list(self.res.y[7, :])

        return

    def getRawResults(self):
        return self.res
