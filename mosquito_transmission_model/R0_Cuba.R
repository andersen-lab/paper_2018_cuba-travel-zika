# read in temperature data and convert to Celsius
w = read.csv('OpenWeatherMap_Havana_2014-2017_MontlyTemp.csv',sep=';')
tempC = w$DailyTemp_Avg - 273.15
tempC = tempC - tempC %% 0.1

# read in Mordecai et al. posterior samples
load('Aegypti_DENV_model_outputs-informative.Rsave')

# get indices for their vectors
temp.ind = match(tempC,temp)
temp.ind[c(2)] <- 178
temp.ind[c(25)] <- 163
temp.ind[c(26)] <- 155
temp.ind[c(38)] <- 178
temp.ind[c(49)] <- 155
temp.ind[c(51)] <- 168

R0.Cuba = R0[temp.ind,]

# export table with R0 values
write.table(R0.Cuba, file="R0.Cuba.txt")

