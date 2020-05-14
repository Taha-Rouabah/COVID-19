"""
v001
"""

import pandas as pd
from datetime import date

class Region:
    """
        Defines the city/province/country of the epidemic by its
        population and daily cumulative cases curve.
    """
    def __init__(self, name = 'Void', nCitizen = 0, first_day = date(2020,1,1), cumulativeCases = [], recoveredCases = [], deathCases = []):
        self.name = name
        self.N = nCitizen
        self.first_day = first_day
        self.rcQ = cumulativeCases
        self.rR = recoveredCases
        self.rD = deathCases


class WebDataReader:
    """    
        This class creates Region object with necessary data 
        inputs: Country (str)
                State (str)
        output: (Region object, dataframe object)        
    """
    def __init__(self, dataset, country, state='nan'):
        self.dataset = dataset
        self.country = country
        self.state = state
        (self.Region,self.data) = self.form(self.country,self.state)
    
    def form(self, country, state):
        data = self.pick(country,state)
        
        # Initialize the province
        Province = Region()
        
        # province's name
        if state == 'nan':
            Province.name = country
        else:
            Province.name = state
        
        # First day of data
        year = int(data.loc[0, 'Date'][0] + data.loc[0, 'Date'][1] + data.loc[0, 'Date'][2] + data.loc[0, 'Date'][3])
        month = int(data.loc[0, 'Date'][5] + data.loc[0, 'Date'][6])
        day = int(data.loc[0, 'Date'][8] + data.loc[0, 'Date'][9])
        
        Province.first_day = date(year,month,day)
        
        Province.rcQ = list(data.loc[:,'Confirmed'].values)
        Province.rR = list(data.loc[:,'Recovered'].values)
        Province.rD = list(data.loc[:,'Deaths'].values)
        return(Province,data)
    
    def pick(self, country, state = 'nan'):    
        """        
            This function picks the relevent country/state's data from the whole dataset        
        """
        # initialize a dataframe
        data = pd.DataFrame()
        
        if state == 'nan':
            index=0
            head = False
            end = False
                
            while not end:
                
                info = self.dataset.loc[index, 'Confirmed'] + self.dataset.loc[index, 'Recovered'] + self.dataset.loc[index, 'Deaths']
                if self.dataset.loc[index, 'Country/Region'] == country and info > 0:
                    head = True
                    data = data.append(self.dataset.iloc[index,:],ignore_index=True)
                    data
                elif head == True:
                    end = True
            
                index += 1
                if index > self.dataset.shape[0]:
                    end =True
        
        else:
            index=0
            head = False
            end = False
                
            while not end:
                
                info = self.dataset.loc[index, 'Confirmed'] + self.dataset.loc[index, 'Recovered'] + self.dataset.loc[index, 'Deaths']
                if self.dataset.loc[index, 'Country/Region'] == country and self.dataset.loc[index, 'Province/State'] == state and info > 0:
                    head = True
                    data = data.append(self.dataset.iloc[index,:],ignore_index=True)
                elif head == True:
                    end = True
            
                index += 1
                if index > self.dataset.shape[0]:
                    end =True
        data = data.drop(columns =['Country/Region','Province/State','Lat', 'Long'],axis=1)
        
        return(data)
        
