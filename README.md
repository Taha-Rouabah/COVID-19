# COVID-19 mathematical modeling

This repository is aimed to contain the programs used for dynamical modeling of COVID-19 spread with SEIQRDP compartmental model. It is described in [this paper](https://www.researchgate.net/publication/341655849_Early_dynamics_of_COVID-19_in_Algeria_a_model-based_study) by Mohamed Taha Rouabah, Abdellah Tounsi, and Nacer Eddine Belaloui.

This repository will also contain an agent-based program for epidemic spread simulation in the near future.

## seiqrdp_model

This package includes computational tools to fit data, simulate and calibrate
parameters of the compartmental epidemiological SEIQRDP model using a genetic
fitting algorithm and parallel computing techniques.

### Architecture

![The architecture of the module](https://github.com/Taha-Rouabah/COVID-19/raw/master/images/SEIRQDPDiag.png)

### How to install

#### Using pip _(preferred way)_

Execute the following line in your command-line interface:

```
pip install seiqrdp-model
```


#### Using setup.py

Download the repository then execute the following line in your command-line interface:

```
python setup.py install
```

> You need to [cd](<https://en.wikipedia.org/wiki/Cd_(command)>) inside the repository folder.

### Quick start

The module is fully documented and equiped with all necessary information
about the classes and methods.
Moreover, you are welcome to check `launcher.py` file which shows an example of
how to use this module.

## Data

- World wide data is sourced and maintained by the [John Hopkins University Center for Systems Science and Engineering](https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv)
- Data for Algeria is provided by [Algerian Health Ministry](http://covid19.sante.gov.dz/carte/) and [CDTA](https://covid19.cdta.dz/dashboard/production/index.php#)

## License

This work is licensed under the [Apache 2.0 License](https://github.com/Taha-Rouabah/COVID-19/blob/master/LICENSE)
