# Estimating the Zika outbreak size in Cuba using travel and local data

The goal of this analysis is to estimate local incidence of Zika in Cuba. We fitted models to estimate local and travel cases for other countries in the Americas and used the posterior relationship between local and travel cases from all other countries in the data set to infer local cases in Cuba. Here, we outline how to run the analysis to recreate these results from the Grubaugh *et al.* manuscript. 

## Getting started

The R scripts are written assuming you have the following folder structure:

```
cuba_outbreak
│   README.md
└─── code
└─── output
└─── data
```
Where all of the MCMC code, analysis code, and processing code is in the 'code' folder. All of the figures generated and resulting .csv files feed to the 'output' folder. All requisite data is in the 'data' folder. All of the R scripts are written assuming you are inside of the 'code' folder. 

### Software and packages

We used R version 3.4.1, "Single Candle". 

\newline R packages necessary for these analyses:

* fda
* VGAM
* mvtnorm
* spatstat

```
install.packages(c('fda', 'VGAM', 'mvtnorm', 'spatstat'))
```

## Analysis

### 1. Processing data
First, you will need to process the data to format it into one data frame with travel cases, local cases, travel population, and local population for each country in the data set. To do this, you can source the data processing file: 

```
souce('process_data.R')
```
This script will produce **../data/zika_processed_16-17.csv**, which will be used in all subsequent scripts. 

### 2. Running MCMC scripts

You will need to run an MCMC script for each country considered in the analysis, which can be done by opening the R script and running the MCMC file. These countries include: Antigua and Barbuda, Bahamas, Bolivia, Brazil, Cayman Islands, Colombia, Costa Rica, Dominica, Dominican Republic, Ecuador, El Salvador, Grenada, Guadeloupe, Guatemala, Guyana, Haiti, Honduras, Jamaica, Martinique, Mexico, Nicaragua, Panama, Puerto Rico, Panama, Trinidad and Tobago, Venezuela, British Virgin Islands, and US Virgin Islands.

Alternatively, you can just source the file, as seen below:

```
source('mcmc_bahamas.R')
```

And repeat for each country with at least a correlation of 0.25 between local and travel cases:

```
source('mcmc_virgin-islands-US.R')
```
These MCMC scripts will produce an .RData file with posterior estimates (**../output/posterior_country.RData**) for the model parameters and the acceptance rate of the MCMC algorithm. The script will also produce two figures:

* **../output/country_case_rate.pdf**: case rate estimates for local cases and travel cases 
* **../output/country_acceptance_rate.pdf**: acceptance rate for iteration steps in the MCMC chain


### 3. Estimating local Zika incidence in Cuba

To estimate local Zika incidence in Cuba using the posterior estimates of travel cases and local cases of the other countries considered in the analysis, we source the R script: 

```
source('process_mcmc_output.R')
```
The script will produce three csv files:

* **../output/local_cases_posterior.csv**: 1,000 samples from posterior for monthly local cases in Cuba (24 month time series)
* **../output/local_cases_log_IQR.csv**: log10 IQR for monthly local cases in Cuba (0.25, 0.75, 0.5)
* **../output/local_cases_log_95.csv**: log10 95% posterior predictive interval for monthly local cases in Cuba (0.025, 0.975, 0.5)

This script will also produce supplementary figures (Fig. S2 and Fig. S3)

* **../output/hist_total_travel.pdf**: Fig. S2, posterior distribution of total travel cases for each country in the analysis with empiricial value
* **../output/hist_total_local_mu.pdf**: Fig. S3, posterior distribution of total local cases for each country in the analysis with empiricial value, with posterior distribution of total local cases in Cuba in blue

## Authors

* **Rachel J. Oidtman**
* T. Alex Perkins

## Acknowledgments
This README.md is part of a larger project from the manuscript:

Grubaugh N.D., Saraf S., Gangavarapu K., Watts A., Tan A.L., ... , Anderson K.G. International travelers and genomics uncover a 'hidden' Zika outbreak. 