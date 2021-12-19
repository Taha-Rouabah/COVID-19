"""
v002:
    - Region class handles its own data efficiently now.
    - Can create a Region from a dataframe.
    - Got rid of WebDataReader class.
"""
__version__ = '2.0'

from datetime import date
from datetime import timedelta
from dateutil.parser import parse


class DataError(Exception):
    """ Exception raised when there are no data to fit that is above
        the threshold.
    """

    def __init__(self, message):
        self.message = message
        pass


class Region:
    """
        Defines the city/province/country of the epidemic by its
        population and daily cumulative cases curve.
    """

    def __init__(self, name='', population_size=0, first_day=date(2020, 1, 1),
                 cumulative_cases=[], recovered_cases=[], death_cases=[]):
        # Name and data.
        self.name = name
        self.N = population_size
        self.rcQ = cumulative_cases
        self.rR = recovered_cases
        self.rD = death_cases
        # Time frame.
        self.first_day = first_day
        self.last_day = self.first_day + timedelta(days=len(self.rcQ))

    def cutoff_data(self, data_size_cutoff):
        """ Cuts off the (training) data to a certain number of days. """
        self.rcQ = self.rcQ[:data_size_cutoff]
        self.rR = self.rR[:data_size_cutoff]
        self.rD = self.rD[:data_size_cutoff]
        self.last_day = self.first_day + timedelta(days=data_size_cutoff)

    @staticmethod
    def from_dataframe(dataframe, region_name, population_size, threshold=0):
        """ Constructs and returns a Region object from a dataframe.
        """

        # Creating initial Region object.
        region = Region(name=region_name, population_size=population_size)

        # Filtering to only region of interest.
        df = dataframe.loc[dataframe['Country/Region'] == region_name]

        # Populating the region object
        region.name = region_name
        region.N = population_size
        region.first_day = parse(dataframe['Date'].iloc[0]).date()
        region.rcQ = df['Confirmed'].tolist()
        region.rR = df['Recovered'].tolist()
        region.rD = df['Deaths'].tolist()

        # Find the first day above threshold for each data.
        first_rcQ = next((i for i, x in enumerate(region.rcQ)
                          if x > threshold), None)
        first_rR = next((i for i, x in enumerate(region.rR)
                         if x > threshold), None)
        first_rD = next((i for i, x in enumerate(region.rD)
                         if x > threshold), None)

        # If at least one data column is not above threshold, raise error.
        if None in (first_rcQ, first_rR, first_rD):
            raise DataError('Could not find data points'
                            ' above specified threshold!')

        # First day to be above all thresholds.
        first_ind = max(first_rcQ, first_rR, first_rD)

        # Dropping data below the threshold.
        region.rcQ = region.rcQ[first_ind:]
        region.rR = region.rR[first_ind:]
        region.rD = region.rD[first_ind:]
        # Re-computing the time boundaries of the data.
        region.first_day = region.first_day + timedelta(days=first_ind+1)
        region.last_day = region.first_day + timedelta(days=len(region.rcQ))

        return region
