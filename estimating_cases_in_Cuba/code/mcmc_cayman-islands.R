# fitting GAMs to Zika travel data and local data

# =============================================================== #
## read in data
dat = read.csv(file = '../data/zika_processed_16-17.csv', stringsAsFactors = F)
dat[is.na(dat)] = 0


countries = unique(dat$country)

require(fda)
require(VGAM)
require(mvtnorm)
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

# function to product index for input into align_index function
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


# log likeliehood function for spline vs. data
# inputs: dat_in = dat; mu_in = multiplier on local spline;
# spline_in = mean monthly spline from align_dates fxn; ctry = country of interest
ll = function(spline_in, mu_in, ctry = cc, dat_in = dat, browse = F){
  if(browse) browser()
  # dat_in, mu_in, spline_in, ctry, ind_in){
  
  travel_ll = dbinom(x = dat_in$travel.cases[which(dat_in$country == ctry)],
                     size = dat_in$travel.pop[which(dat_in$country == ctry)],
                     prob = inv_logit(spline_in), log = T)
  
  local_ll = dbinom(x = dat_in$local.cases[which(dat_in$country == ctry)],
                    size = dat_in$local.pop[which(dat_in$country == ctry)],
                    prob = inv_logit(spline_in + mu_in), log = T)
  
  # travel_ll[which(is.infinite(travel_ll))] = (-.Machine$double.xmax / 1e100)
  # local_ll[which(is.infinite(local_ll))] = (-.Machine$double.xmax / 1e100)
  
  return(0.01 * sum(travel_ll, local_ll, na.rm = T))
  # return(0.001 * sum(travel_ll, na.rm = T))
}


# logit fxn
logit = function(x_in){
  return(log(x_in) - log(1-x_in))
}


# inverse logit fxn
inv_logit = function(x_in){
  return(exp(x_in) / (exp(x_in) + 1))
}


# proposal draw
proposal_draw = function(param_old, sigma_in){
  return(rmvnorm(1, param_old, makePositiveDefinite(sigma_in)))
}


# proposal probability
proposal_prob = function(param_old, param_new, sigma_in){
  return(dmvnorm(param_new, param_old, makePositiveDefinite(sigma_in), log = T))
}


# =============================================================== #
## determine correlations of 

cor_by_country = rep(NA, length(countries))
l = 1
for(cc in countries){
  cor_by_country[l] = cor(dat$local.cases[which(dat$country == cc)], dat$travel.cases[which(dat$country == cc)])
  l = l + 1
}

countries_keep = countries[which(cor_by_country > 0.25)]


# =============================================================== #
## mcmc inputs

cc = countries_keep[5]

# knots per year
k = 4
spl_basis = generate_basis(k_in = k)
ind_align = produce_ind_align()

# initial params
local_travel_ratio = (mean(dat$local.cases[which(dat$country == cc)]) / mean(dat$local.pop[which(dat$country == cc)])) / 
  (mean(dat$travel.cases[which(dat$country == cc)]) / mean(dat$travel.pop[which(dat$country == cc)]))

initial_mu = log(local_travel_ratio)
initial_y = logit(c(0.0001, 0.0002, 0.0005, 0.0002, 0.0001, 0.0001, 0.0001, 0.0001))

initial_param = c(initial_mu, initial_y)

sigma_param = diag(rep(1, length(initial_param)))

# =============================================================== #
## mcmc

f = paste0('../output/posterior_', cc, '.RData')

if(file.exists(f)){
  load(f)
  
}else{
  #################################
  # start of burn in loop 
  
  burn_in = 1000
  
  p_spline = generate_spline(y_in = initial_param[2:length(initial_param)])
  monthly_p_spline = align_dates(spline_in = p_spline)
  LL_current = ll(spline_in = monthly_p_spline, mu_in = initial_param[1])
  
  # initalize param_current
  param_current = initial_param
  
  for(ii in 2:burn_in){
    param_proposed = proposal_draw(param_current, sigma_param)
    
    p_spline = generate_spline(y_in = param_proposed[2:length(initial_param)])
    monthly_p_spline = align_dates(spline_in = p_spline)
    LL_proposed = ll(spline_in = monthly_p_spline, mu_in = param_proposed[1])
    
    acceptance_prob = LL_proposed - LL_current + 
      proposal_prob(param_current, param_proposed, sigma_param) - 
      proposal_prob(param_proposed, param_current, sigma_param)
    
    if(runif(1) < exp(acceptance_prob)){
      param_current = param_proposed
      LL_current = LL_proposed
    }
  }
  param_burn = param_current
  
  
  #################################
  ## adaptive loop
  p_spline = generate_spline(y_in = param_burn[2:length(initial_param)])
  monthly_p_spline = align_dates(spline_in = p_spline)
  LL_current = ll(spline_in = monthly_p_spline, mu_in = param_burn[1])
  
  # define number of times to repeat adaptive sampling and chain length
  loop_reps = 3
  chain_length = 1000
  
  # allocate array to store all of the param values
  param_chain_cov = array(0, dim = c(chain_length, length(initial_param), loop_reps))
  param_chain_cov[1,,1] = param_burn
  
  acceptance_rate_cov = rep(NA, loop_reps)
  
  param_adapt = matrix(NA, nrow = loop_reps + 1, ncol = length(initial_param))
  param_adapt[1, ] = param_burn
  
  for(rr in 1:loop_reps){
    param_current = param_adapt[rr, ]
    accepted = rejected = 0
    
    for(ii in 2:chain_length){
      
      param_proposed = proposal_draw(param_current, sigma_param)
      
      p_spline = generate_spline(y_in = param_proposed[2:length(initial_param)])
      monthly_p_spline = align_dates(spline_in = p_spline)
      LL_proposed = ll(spline_in = monthly_p_spline, mu_in = param_proposed[1])
      
      acceptance_prob = LL_proposed - LL_current + 
        proposal_prob(param_current, param_proposed, sigma_param) - 
        proposal_prob(param_proposed, param_current, sigma_param)
      
      if(runif(1) < exp(acceptance_prob)){
        param_current = param_proposed
        LL_current = LL_proposed
        accepted = accepted + 1
      }else{
        rejected = rejected + 1
      }
      
      param_chain_cov[ii,,rr] = param_current
    }
    param_adapt[rr + 1, ] = param_current
    acceptance_rate_cov[rr] = accepted / (accepted + rejected)
    
    sigma_param = (2.38 ^ 2) * cov(param_chain_cov[,,rr]) / length(initial_param)
  }
  
  param_cov = param_current
  
  
  #################################
  ## final mcmc chain
  
  p_spline = generate_spline(y_in = param_cov[2:length(initial_param)])
  monthly_p_spline = align_dates(spline_in = p_spline)
  LL_current = ll(spline_in = monthly_p_spline, mu_in = param_cov[1])
  
  # define number of interations of MCMC and frequency of recording
  chain_length = 2000
  chain_freq = 20
  
  # allocate vectors to stor chains for each parameter
  param_chain = matrix(NA, nrow = chain_length, ncol = length(initial_param))
  
  acceptance_rate = numeric(0)
  
  for(ii in 1:chain_length){
    accepted = rejected = 0
    
    for(jj in 1:chain_freq){
      param_proposed = proposal_draw(param_current, sigma_param)
      
      p_spline = generate_spline(param_proposed[2:length(initial_param)])
      monthly_p_spline = align_dates(spline_in = p_spline)
      LL_proposed = ll(spline_in = monthly_p_spline, mu_in = param_proposed[1])
      
      acceptance_prob = LL_proposed - LL_current + 
        proposal_prob(param_current, param_proposed, sigma_param) - 
        proposal_prob(param_proposed, param_current, sigma_param)
      
      if(runif(1) < exp(acceptance_prob)){
        param_current = param_proposed
        LL_current = LL_proposed
        accepted = accepted + 1
      }else{
        rejected = rejected + 1
      }
    }
    
    param_chain[ii,] = param_current
    acceptance_rate = c(acceptance_rate, (accepted / (accepted + rejected)))
    
    # write to file
    save(acceptance_rate, param_chain, file = f)
  }
  
}


# =============================================================== #
## produce output


spline_matrix_travel = matrix(NA, nrow = nrow(param_chain), ncol = nrow(spline_df))
for(ss in 1:nrow(spline_matrix_travel)){
  spline_matrix_travel[ss, ] = generate_spline(y_in = param_chain[ss, 2:length(initial_param)])
}

spline_matrix_local = matrix(NA, nrow = nrow(param_chain), ncol = nrow(spline_df))
for(ss in 1:nrow(spline_matrix_local)){
  spline_matrix_local[ss, ] = spline_matrix_travel[ss, ] + param_chain[ss, 1] 
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


spline_travel_CI = inv_logit(write_CI(spline_matrix_travel))
spline_local_CI = inv_logit(write_CI(spline_matrix_local))

f = paste0('../output/', cc, '_case_rate.pdf')
pdf(f, width = 6, height = 5)
layout(matrix(1:2, 2, 1))
par(mar = c(1, 2, 0.5, 0.5), oma = c(2,1,0,0))
plot(dat$travel.cases[which(dat$country == cc)] / dat$travel.pop[which(dat$country == cc)], type = 'l',
     col = rgb(1,0,0,0.75), lwd = 3, xaxt = 'n', ylim = c(0, 0.002))
lines(dat$local.cases[which(dat$country == cc)] / dat$local.pop[which(dat$country == cc)], type = 'l',
      col = rgb(0,0,1,0.75), lwd = 3)


plot(-100, -100, xlim = c(0, 731), ylim = c(0, 0.002))
polygon(x = c(1:731, rev(1:731)),
        y = c(spline_local_CI[1,], rev(spline_local_CI[2,])), border = NA,
        col = rgb(0,0,1,0.3))
lines(spline_local_CI[3,], col = rgb(0,0,1,0.75), lwd = 4)

polygon(x = c(1:731, rev(1:731)),
        y = c(spline_travel_CI[1,], rev(spline_travel_CI[2,])), border = NA,
        col = rgb(1,0,0,0.3))
lines(spline_travel_CI[3,], col = rgb(1,0,0,0.75), lwd = 4)

mtext(side = 1, text = 'weeks since jan. 1 2016', line = 1.75)
mtext(side = 2, text = 'case rate per population', outer = T)
dev.off()


f = paste0('../output/', cc, '_acceptance_rate.pdf')
pdf(f, width = 6, height = 3)
par(mar = c(2,2,0.5,0.5), oma = c(1,1,0,0))
layout(1)
plot(acceptance_rate, type = 'l', xlab = '', ylab = '')
mtext(side = 1, text = 'iterations', line  = 1.75)
mtext(side = 2, text = 'acceptance rate', line  = 1.75)
dev.off()
