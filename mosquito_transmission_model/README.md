* Summary: Monthly Aedes aegypti transmission potential in Cuba.*

To assess yearly and seasonal variations in Ae. aegypti transmission potential for dengue and Zika virus in Cuba, we used a temperature-dependent model of transmission using a previously developed R0 framework (https://doi.org/10.1371/journal.pntd.0005568). Using hourly temperature data obtained from OpenWeatherMap (https://openweathermap.org/), we calculated monthly mean temperatures in Havana, Cuba. These temperatures were used to calculate monthly R0 as estimated by Mordecai et al. (https://figshare.com/s/b79bc7537201e7b5603f). Doing so for 5,000 samples from the posterior of temperature-R0 relationships and normalizing between 0 and 1 yielded a description of relative Ae. aegypti transmission potential per month in Havana, Cuba during 2014-2017.

#### Data
OpenWeatherMap_Havana_2014-2017_MontlyTemp.csv
Hourly temperature data from Havana, Cuba, 2014-2017, was collected from https://openweathermap.org/, and converted to monthly average, min, max, and standard deviations. Months from 2014-2017 were numbered 1-58.

#### R Code
R0_Cuba.R
R code to calculate monthly mean temperatures and calculate R0.

Aegypti_DENV_model-output-informative.Rsave
R model to calculate monthly R0 values based on monthly mean temperatures.

#### Raw output data
R0_Cuba.txt 
Raw output data file with 5,000 monthly R0 calculations for Havana, Cuba during 2014-2017.

#### Processed output data
R0_Cuba_output.xlsx
Processed output data with calculated monthly average R0 values, standard deviation, N, minimum R0, maximum R0, 95% confidence intervals, and normalized R0 values. 
