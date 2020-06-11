# SEIQRDP_model

A COVID-19 modeling that uses a generalized SEIR model.

A tailored set of codes has been developed under Python 3.7 to:

1. download data from the online COVID-19 data repository  
2. calibrate the model by applying a genetic fitting algorithm to the real data
3. perform forecast modeling
4. and provides parallel computing techniques.
<!---
(determine the optimum fitting depth using cross-validation method)
-->

## Architecture

![The architecture of the module](/covid-19/SEIQRDP_model/images/Diagram.png)


## How to use it?

The module is full documented. Feel free to read all necessary information
about the classes and methods. 

You are welcome to check *launcher.py* file which shows an example of
how to use this module.

## Data

Data is provided by [github user](https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv)
