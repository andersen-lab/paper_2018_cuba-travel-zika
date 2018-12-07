## processing Cuba data to put into a usable format

# =============================================================== # 
## read in data

dat.loc = read.csv('../data/zika_local_16-17.csv')
dat.trav = read.csv('../data/zika_travel_16-17.csv')
pop.loc = read.csv('../data/zika_local_pop_16-17.csv')
pop.trav = read.csv('../data/zika_travel_pop_16-17.csv', stringsAsFactors = F)
countries = as.character(dat.loc$Local.cases[2:nrow(dat.loc)])


# =============================================================== # 
## processing data

# remove commas from pop.trav
for(ii in 2:ncol(pop.trav)){
  pop.trav[,ii] = gsub( ',', '', pop.trav[,ii])
  pop.trav[,ii] = as.numeric(pop.trav[,ii])
}

# create new data frame
dat.new = data.frame(
  country = rep(countries, 12 * 2),
  month = rep(rep(1:12, each = length(countries))),
  year = c(rep(2016, length(1:12) * length(countries)),
           rep(2017, length(1:12) * length(countries))),
  local.cases = rep(NA, length(1:12) * length(countries) * 2),
  local.pop = rep(NA, length(1:12) * length(countries) * 2),
  travel.cases = rep(NA, length(1:12) * length(countries) * 2),
  travel.pop = rep(NA, length(1:12) * length(countries) * 2)
)

# allocate local cases
for(ii in 1:length(countries)){
  # local cases
  dat.new$local.cases[which(dat.new$country == countries[ii])] = 
    dat.loc[which(dat.loc$Local.cases == countries[ii]),][-1]
  
  # local pop
  dat.new$local.pop[which(dat.new$country == countries[ii])] = 
    pop.loc[which(pop.loc$Population == countries[ii]),][-1]
  
  # travel cases
  dat.new$travel.cases[which(dat.new$country == countries[ii])] = 
    dat.trav[which(dat.trav$Travel.cases == countries[ii]), ][-1]
  
  # travel pop
  if(ii != 17){
    dat.new$travel.pop[which(dat.new$country == countries[ii])] = 
      pop.trav[which(pop.trav$Total.travel.volume == countries[ii]), ][-1]
  }else{
    dat.new$travel.pop[which(dat.new$country == countries[ii])] = 
      pop.trav[which(pop.trav$Total.travel.volume == 'Cuba - IATA + charter'), ][-1]
  }
}
dat.new$local.cases = unlist(dat.new$local.cases)
dat.new$local.pop = unlist(dat.new$local.pop)
dat.new$travel.cases = unlist(dat.new$travel.cases)
dat.new$travel.pop = unlist(dat.new$travel.pop)


# remove commas from pop.trav
write.csv(dat.new, row.names = F, file = '../data/zika_processed_16-17.csv')
# =============================================================== # 