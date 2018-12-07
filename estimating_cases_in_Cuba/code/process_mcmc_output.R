# processing mcmc output and determining local and travel cases for Cuba, and other countries

# =============================================================== #
## read in data
dat = read.csv(file = '../data/zika_processed_16-17.csv', stringsAsFactors = F)
dat[is.na(dat)] = 0


countries = unique(dat$country)

## get rid of white space in country names
countries_squish = countries
countries_squish[2] = 'Antigua_Barbuda'
countries_squish[10] = 'Bonaire_SaintEustatius_andSaba'
countries_squish[13] = 'CaymanIslands'
countries_squish[16] = 'CostaRica'
countries_squish[20] = 'DominicanRep'
countries_squish[22] = 'ElSalvador'
countries_squish[23] = 'FrenchGuiana'
countries_squish[37] = 'PuertoRico'
countries_squish[38] = 'StBarthelemy'
countries_squish[39] = 'StKitts_Nevis'
countries_squish[40] = 'StLucia'
countries_squish[41] = 'StMaarten_DutchPart'
countries_squish[42] = 'StVincent_Grenad'
countries_squish[44] = 'Trinidad_Tobago'
countries_squish[45] = 'Turks_Caicos'
countries_squish[48] = 'VirginIslands_GB'
countries_squish[49] = 'VirginIslands_US'


require(fda)
require(VGAM)
require(mvtnorm)
require(spatstat)

## determine how many total travel cases for each country
total_travel = rep(NA, length(countries))
names(total_travel) = countries
for(cc in countries){
  total_travel[cc] = sum(dat$travel.cases[which(dat$country == cc)])
}

## determine correlations between local and travel cases
cor_by_country = rep(NA, length(countries))
l = 1
for(cc in countries){
  cor_by_country[l] = cor(dat$local.cases[which(dat$country == cc)], dat$travel.cases[which(dat$country == cc)])
  l = l + 1
}

## keep countries with at least 0.25 correlation between local and travel cases
countries_keep = countries[which(cor_by_country > 0.25)]
countries_keep_squish = countries_squish[which(cor_by_country > 0.25)]

## loading posteriors for mu
for(ff in 1:length(countries_keep)){
  f = paste0('../output/posterior_', countries_keep[ff], '.RData')
  load(f)
  assign(paste0('mu_', countries_keep_squish[ff]), param_chain[, 1])
  assign(paste0('param_chain_', countries_keep_squish[ff]), param_chain)
  rm(param_chain)
}


# =============================================================== #
## data storage

months_2016 = c(rep(1, 31), rep(2, 29), rep(3, 31), rep(4, 30), rep(5, 31),
                rep(6, 30), rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31),
                rep(11, 30), rep(12, 31))
months_2017 = c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4, 30), rep(5, 31),
                rep(6, 30), rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31),
                rep(11, 30), rep(12, 31))
days_2016 = c(1:31, 1:29, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31,
              1:30, 1:31)
days_2017 = c(1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31,
              1:30, 1:31)

spline_df = data.frame(year = c(rep(2016, 366), rep(2017, 365)),
                       month = c(months_2016, months_2017),
                       day = c(days_2016, days_2017), 
                       country = rep(NA, length(c(days_2016, days_2017))),
                       daily.spline = rep(NA, length(c(days_2016, days_2017))))


# =============================================================== #
## functions

# generate spline function
generate_spline = function(y_in, spline_len = nrow(spline_df), bs_in = spl_basis, browse = F){
  if(browse) browser()
  
  fd_params = y_in
  fd_bs = fd(
    coef = fd_params,
    basisobj = bs_in
  )
  
  spline_tmp = predict(fd_bs, 1:spline_len)
  return(spline_tmp)
}


# generate basis function
generate_basis = function(k_in, num_years = 2, spline_len = nrow(spline_df)){
  basis_tmp = create.bspline.basis(
    rangeval = c(1, spline_len),
    nbasis = k_in * num_years,
    norder = k_in
  )
  return(basis_tmp)
}


# logit fxn
logit = function(x_in){
  return(log(x_in) - log(1-x_in))
}


# inverse logit fxn
inv_logit = function(x_in){
  return(exp(x_in) / (exp(x_in) + 1))
}


# function to produce index for input into align_index function
produce_ind_align = function(spline_df_in = spline_df, browse = F){
  if(browse) browser()
  yr_tmp = c(rep(2016, 12), rep(2017, 12))
  month_tmp = c(1:12, 1:12)
  
  ind_tmp = list()
  for(ff in 1:(length(unique(spline_df_in$month))* 2)){
    ind_tmp[[ff]] = which(spline_df_in$year == yr_tmp[ff] &
                            spline_df_in$month == month_tmp[ff])
  }
  
  return(ind_tmp)
}


# align dates function
align_dates = function(spline_in, dat_in = dat, ind_in = ind_align, browse = F){
  if(browse) browser()
  spline_monthly = sapply(1:(length(unique(dat_in$month)) * 2), function(ff) mean(spline_in[ind_in[[ff]]]))
  return(spline_monthly)
}


# write 95% CI
write_CI = function(param_posterior, lwr = 0.025, upr = 0.975){
  
  mat = matrix(0, ncol = ncol(param_posterior), nrow = 3)
  
  for(ii in 1:ncol(param_posterior)){
    mat[1, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[1]
    mat[2, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[2]  
    mat[3, ii] = median(param_posterior[, ii], na.rm = T)
  }
  return(mat)
}


# =============================================================== #
## put mus into a matrix

mus = ls()[substr(ls(),0,3)=='mu_']

mu = matrix(0,length(mu_Brazil),length(mus))
for(ii in 1:length(mus)){
  eval(parse(text=paste('
                        mu[,ii] = ', mus[ii], sep='')))
}

colnames(mu) = sapply(strsplit(mus,'_'),function(ii)ii[2])

mu_means = colMeans(mu)


# =============================================================== #
## calculate cumulative incidence for local and travel cases in data

travel_cum = rep(NA, length(mus))
local_cum = rep(NA, length(mus))
for(ii in 1:length(mus)){
  travel_cum[ii] = sum(dat$travel.cases[which(dat$country == countries_keep[ii])]) /
    mean(dat$travel.pop[which(dat$country == countries_keep[ii])])

  local_cum[ii] = sum(dat$local.cases[which(dat$country == countries_keep[ii])]) /
    mean(dat$local.pop[which(dat$country == countries_keep[ii])])
}

names(local_cum) = sapply(strsplit(mus,'_'),function(ii)ii[2])
names(travel_cum) = sapply(strsplit(mus,'_'),function(ii)ii[2])


# =============================================================== #
## calculate cumulative incidence for local and travel cases in data

# inputs for generating spline; k = knots per year, spl_basis = spline basis
k = 4
spl_basis = generate_basis(k_in = k)
ind_align = produce_ind_align()


# daily travel spline
spline_array_travel = array(NA, c(nrow(param_chain_Antigua_Barbuda), nrow(spline_df), length(countries_keep)))
for(cc in 1:length(countries_keep_squish)){
  spline_matrix_travel = matrix(NA, nrow = nrow(param_chain_Antigua_Barbuda), ncol = nrow(spline_df))
  
  for(ss in 1:nrow(spline_matrix_travel)){
    eval(parse(text = paste0(
      'x = generate_spline(y_in = param_chain_', countries_keep_squish[cc], '[ss, 2:9])'
    )))
    # x = generate_spline(y_in = param_chain[ss, 2:8])
    spline_matrix_travel[ss, ] = x
  }
  spline_array_travel[,,cc] = spline_matrix_travel
}


# produce monthly travel spline
spline_array_travel_monthly = array(NA, c(nrow(param_chain_Antigua_Barbuda), 24, length(countries_keep)))

for(cc in 1:length(countries_keep)){
  for(ss in 1:nrow(spline_matrix_travel)){
    tmp = align_dates(spline_in = spline_array_travel[ss,, cc])
    spline_array_travel_monthly[ss,,cc] = tmp
  }
}


# calculate cumulative incidence
mat_cum_travel = matrix(NA, nrow = nrow(param_chain_Antigua_Barbuda), ncol = length(countries_keep))
for(ii in 1:length(countries_keep_squish)){
  for(ss in 1:nrow(spline_matrix_travel)){
    mat_cum_travel[ss, ii] = sum(inv_logit(spline_array_travel_monthly[ss,,ii]) * 
                                   dat$travel.pop[which(dat$country == countries_keep[ii])]) / 
      mean(dat$travel.pop[which(dat$country == countries_keep[ii])])
  }
}
colnames(mat_cum_travel) = sapply(strsplit(mus,'_'),function(ii)ii[2])


# calculate daily local spline
spline_array_local = array(NA, c(nrow(param_chain_Antigua_Barbuda), nrow(spline_df), length(countries_keep)))
for(cc in 1:length(countries_keep_squish)){
  spline_matrix_local = matrix(NA, nrow = nrow(param_chain_Antigua_Barbuda), ncol = nrow(spline_df))
  
  for(ss in 1:nrow(spline_matrix_local)){
    # spline_matrix_local[ss, ] = spline_array_travel[ss,,cc] + param_chain_Antigua_Barbuda[ss, 1]
    
    eval(parse(text = paste0(
      'x = spline_matrix_travel[ss, ] = spline_array_travel[ss,,cc] + param_chain_', countries_keep_squish[cc], '[ss, 1]'
    )))
    # x = generate_spline(y_in = param_chain[ss, 2:8])
    spline_matrix_local[ss, ] = x
  }
  spline_array_local[,,cc] = spline_matrix_local
}


# produce monthly local spline
spline_array_local_monthly = array(NA, c(nrow(param_chain_Antigua_Barbuda), 24, length(countries_keep)))

for(cc in 1:length(countries_keep)){
  for(ss in 1:nrow(spline_matrix_travel)){
    tmp = align_dates(spline_in = spline_array_local[ss,, cc])
    spline_array_local_monthly[ss,,cc] = tmp
  }
}


# calculate cumulative incidence
mat_cum_local = matrix(NA, nrow = nrow(param_chain_Antigua_Barbuda), ncol = length(countries_keep))
for(ii in 1:length(countries_keep_squish)){
  for(ss in 1:nrow(spline_matrix_travel)){
    mat_cum_local[ss, ii] = sum(inv_logit(spline_array_local_monthly[ss,,ii]) * 
                                  dat$local.pop[which(dat$country == countries_keep[ii])]) / 
      mean(dat$local.pop[which(dat$country == countries_keep[ii])])
  }
}
colnames(mat_cum_local) = sapply(strsplit(mus,'_'),function(ii)ii[2])


# =============================================================== #
## figure out which countries have the true value in the 95% CI

# fxn to calculate p value 
calc_p = function(posterior_AR_in, dat_in = dat, ctry_in, trav = F, loc = F, browse = F){
  if(browse) browser()
  if(trav == T){
    v = sum(dat_in$travel.cases[which(dat_in$country == ctry_in)])
    dens = posterior_AR_in * mean(dat_in$travel.pop[which(dat_in$country == ctry_in)])
  }else{
    v = sum(dat_in$local.cases[which(dat_in$country == ctry_in)])
    dens = posterior_AR_in * mean(dat_in$local.pop[which(dat_in$country == ctry_in)])
  }
  quant = sum(dens <= v) / length(dens)
  return(quant)
}

# determining which countries have true value in 95% CI for local and travel cases
local_pValue = rep(NA, length(countries_keep))
for(cc in 1:length(countries_keep)){
  local_pValue[cc] = calc_p(posterior_AR_in = mat_cum_local[,cc], ctry_in = countries_keep[cc], loc = T)
}
travel_pValue = rep(NA, length(countries_keep))
for(cc in 1:length(countries_keep)){
  travel_pValue[cc] = calc_p(posterior_AR_in = mat_cum_travel[,cc], ctry_in = countries_keep[cc], trav = T)
}
which(local_pValue < 0.025 | local_pValue > 0.975)
which(travel_pValue < 0.025 | travel_pValue > 0.975)
countries_keep_pValue = countries_keep[-which(travel_pValue < 0.025 | travel_pValue > 0.975)]
mat_cum_local = mat_cum_local[, -which(travel_pValue < 0.025 | travel_pValue > 0.975)]
mat_cum_travel = mat_cum_travel[, -which(travel_pValue < 0.025 | travel_pValue > 0.975)]
mu = mu[, -which(travel_pValue < 0.025 | travel_pValue > 0.975)]

mat_cum_local_all_countries = as.numeric(mat_cum_local)
mat_cum_travel_all_countries = as.numeric(mat_cum_travel)
mu_all_countries = as.numeric(mu)


# =============================================================== #
## probabilistically get ARt vs. mu

# load in travel only spline estimates for Cuba 
load('../output/posterior_travelOnly_Cuba.RData')
param_chain_Cuba = param_chain; rm(param_chain)


# create probability distributions for all other countries posteriors 
w = as.owin(list(xrange = inv_logit(c(min(mu_all_countries), max(mu_all_countries))),
                 yrange = c(min(mat_cum_travel_all_countries), max(mat_cum_travel_all_countries))))
xy = cbind(inv_logit(mu_all_countries), mat_cum_travel_all_countries)
m = density(as.ppp(X = xy, W = w), dimyx = 128 * 10)
pr = matrix(NA, nrow = m$dim[1], ncol = m$dim[1])
for(ii in 1:m$dim[1]){
  pr[ii, ] = m$v[ii,] / sum(m$v[ii,])
}


# function to input travel AR and output local AR
approx_cdf = function(ARt_in, pr_in = pr, dens_in = m, x_in = runif(1,0,1), browse = F){
  if(browse) browser()
  
  i = which.min(abs(dens_in$yrow - ARt_in))
  fxn = approx(cumsum(pr_in[i, ]), n = nrow(pr_in))
  y = fxn$y
  j = fxn$x[which.min(abs(y - x_in))]
  
  return(dens_in$xcol[j])
}


# calculating travel spline and travel AR
spline_matrix_travel_Cuba = matrix(NA, nrow = nrow(param_chain_Cuba), ncol = nrow(spline_df))
spline_matrix_travel_monthly_Cuba = matrix(NA, nrow = nrow(param_chain_Cuba), ncol = 24)
ARt_Cuba = rep(NA, nrow(spline_matrix_travel_Cuba))
for(ss in 1:nrow(spline_matrix_travel_Cuba)){
  spline_matrix_travel_Cuba[ss, ] = generate_spline(y_in = param_chain_Cuba[ss, ])
  spline_matrix_travel_monthly_Cuba[ss, ] = align_dates(spline_in = spline_matrix_travel_Cuba[ss,])
  ARt_Cuba[ss] = sum(inv_logit(spline_matrix_travel_monthly_Cuba[ss, ]) * dat$travel.pop[which(dat$country == 'Cuba')]) / 
    mean(dat$travel.pop[which(dat$country == 'Cuba')])
}


# calculate ARl for Cuba
spline_matrix_local_Cuba = matrix(NA, nrow = nrow(param_chain_Cuba), ncol = nrow(spline_df))
spline_matrix_local_monthly_Cuba = matrix(NA, nrow = nrow(param_chain_Cuba), ncol = 24)
mu_Cuba = rep(NA, nrow(param_chain_Cuba))
ARl_Cuba_from_mu = rep(NA, nrow(param_chain_Cuba))
for(ss in 1:nrow(spline_matrix_local_Cuba)){
  mu_Cuba[ss] = approx_cdf(ARt_in = ARt_Cuba[ss], browse = F)
  spline_matrix_local_Cuba[ss, ] = spline_matrix_travel_Cuba[ss, ] + mu_Cuba[ss] 
  spline_matrix_local_monthly_Cuba[ss, ] = align_dates(spline_in = spline_matrix_local_Cuba[ss, ])
  ARl_Cuba_from_mu[ss] = sum(inv_logit(spline_matrix_local_monthly_Cuba[ss, ]) * dat$local.pop[which(dat$country == 'Cuba')]) / 
    mean(dat$local.pop[which(dat$country == 'Cuba')])
}

# =============================================================== #
# producing figures and output files

local_Cuba_IQR = write_CI(param_posterior = inv_logit(spline_matrix_local_monthly_Cuba), lwr = 0.25, upr = 0.75) * 
  mean(dat$local.pop[which(dat$country == 'Cuba')])
local_Cuba_95 = write_CI(param_posterior = inv_logit(spline_matrix_local_monthly_Cuba)) * 
  mean(dat$local.pop[which(dat$country == 'Cuba')])
spline_to_plot = inv_logit(spline_matrix_local_monthly_Cuba)
log_spline_to_plot = log10((spline_to_plot) * mean(dat$local.pop[which(dat$country == 'Cuba')]) + 1)
log_spline_to_plot_IQR = write_CI(param_posterior = log_spline_to_plot, lwr = 0.25, upr = 0.75)
log_spline_to_plot_95 = write_CI(param_posterior = log_spline_to_plot)


# save estimates 
local_cases_posterior = spline_to_plot * mean(dat$local.pop[which(dat$country == 'Cuba')])
colnames(local_cases_posterior) = c('01_16','02_16','03_16','04_16','05_16','06_16','07_16','08_16','09_16',
                                    '10_16','11_16','12_16','01_17','02_17','03_17','04_17','05_17','06_17',
                                    '07_17','08_17','09_17','10_17','11_17','12_17')
write.csv(local_cases_posterior, file = '../output/local_cases_posterior.csv',
          row.names = F)

rownames(log_spline_to_plot_IQR) = c('lwr', 'upr', 'median')
write.csv(t(log_spline_to_plot_IQR), file = '../output/local_cases_log_IQR.csv',
          row.names = F)

rownames(log_spline_to_plot_95) = c('lwr', 'upr', 'median')
write.csv(t(log_spline_to_plot_95), file = '../output/local_cases_log_95.csv',
          row.names = F)


# fig S2
pdf('../output/hist_total_travel.pdf', height = 7, width = 8)
layout(matrix(1:25, 5, 5))
par(mar = c(3,2,1,0.5), oma = c(1,1.5,3,0))
for(ii in 1:length(countries_keep_pValue)){
  if(sum(which(mat_cum_travel[,ii] > 0.1)) == 0){
    hist(mat_cum_travel[, ii] * mean(dat$travel.pop[which(dat$country == countries_keep_pValue[ii])]), 100,
         xlab = '', ylab = '', main = countries_keep_pValue[ii], border = rgb(0,0,0,0.3),
         col = rgb(0,0,0,0.5))
  }else{
    tmp = mat_cum_travel[,ii]
    tmp = tmp[-which(tmp > 0.1)]
    hist(tmp * mean(dat$travel.pop[which(dat$country == countries_keep_pValue[ii])]), 100,
         xlab = '', ylab = '', main = countries_keep_pValue[ii], border = rgb(0,0,0,0.3),
         col = rgb(0,0,0,0.5))
  }
  abline(v = sum(dat$travel.cases[which(dat$country == countries_keep_pValue[ii])]), col = 'red')
}
hist(ARt_Cuba * mean(dat$travel.pop[which(dat$country == 'Cuba')]), 100,
     main = 'Cuba', xlab = '', ylab = '',
     border = rgb(0,0,1,0.3), col = rgb(0,0,1,0.5))
abline(v = sum(dat$travel.cases[which(dat$country == 'Cuba')]), col = 'red')
mtext(side = 1, text = 'Travel cases', outer = T)
mtext(side = 2, text = 'Frequency', outer = T)
mtext(side = 3, text = 'Posterior predictions of total travel incidence', outer = T, line = 0.5)
dev.off()


# fig S3
pdf('../output/hist_total_local_mu.pdf', height = 7, width = 8)
layout(matrix(1:25, 5, 5))
par(mar = c(3,2,1,0.5), oma = c(1,1.5,3,0))
for(ii in 1:length(countries_keep_pValue)){
  hist(mat_cum_local[, ii] * mean(dat$local.pop[which(dat$country == countries_keep_pValue[ii])]), breaks = 50, 
       xlab = '', ylab = '', main = countries_keep_pValue[ii], border = rgb(0,0,0,0.3),
       col = rgb(0,0,0,0.5))
  abline(v = sum(dat$local.cases[which(dat$country == countries_keep_pValue[ii])]), col = 'red')
}
hist(ARl_Cuba_from_mu * mean(dat$local.pop[which(dat$country == 'Cuba')]), 50,
     main = 'Cuba', xlab = '', ylab = '',
     border = rgb(0,0,1,0.3), col = rgb(0,0,1,0.5))
mtext(side = 1, text = 'Local cases', outer = T)
mtext(side = 2, text = 'Frequency', outer = T)
mtext(side = 3, text = 'Posterior predictions of total local incidence', outer = T, line = 0.5)
dev.off()

