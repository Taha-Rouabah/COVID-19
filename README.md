# COVID-19 mathematical modeling

This repository is aimed to contain the programs used for dynamical modeling of COVID-19 spread with SEIQRDP compartmental model.

This reposetory will also contain an agent-based program for epidemic spread simulation.

# SEIQRDP_model

This module includes computational tools to fit data, simulate and calibrate
parameters of the compartmental epedimiological SEIQRDP model using genetic 
optimization algorithm and parallel computing techniques.

## Architecture

![The architecture of the module](/images/Diagram.png)


## How to use it?

The module is full documented. Feel free to read all necessary information
about the classes and methods. 

You are welcome to check *launcher.py* file which shows an example of
how to use this module.

## Data

*World wide data is provided by [github user](https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv)
*Data of Algeria is provided by [Algerian official Health Ministry](http://covid19.sante.gov.dz/carte/) and [CDTA](https://covid19.cdta.dz/dashboard/production/index.php#)
