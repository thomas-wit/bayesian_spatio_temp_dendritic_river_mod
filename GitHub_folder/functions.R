## Function that calcluates gelmen Ruben diagnostic given three lists of parameter chains 
gelmanrubin_table <- function(list1,list2,list3){
  for(name in names(list1)){
    if(!is.null(dim(list1[[name]]))){
      print(name)
      print(gelman.diag(mcmc.list(mcmc(t(list1[[name]])),mcmc(t(list2[[name]])),mcmc(t(list3[[name]])))))
    } else {
      print(name)
      print(gelman.diag(mcmc.list(mcmc(list1[[name]]),mcmc(list2[[name]]),mcmc((list3[[name]])))))
    }
  }
}

## This calclulates a table with fatery lewis diags, and effective samp size given a vector of MCMC chain
MCMCdiag <- function(chain){
  param_name <- chain[1]
  chain <- as.numeric(chain[-c(1,length(chain))])
  rd1 <- raftery.diag(chain)$resmatrix
  rd2 <- t(as.matrix(raftery.diag(chain)$params))
  colnames(rd1) <- c('Burn in', 'Samples Desired', 'Samples Needed w/out Cor','Dependence Factor')
  colnames(rd2) <- c('Accuracy', 'Probability', 'Quantile')
  cbind(Param_name = param_name,Eff_Samp_Size = effectiveSize(chain),rd1,rd2)
}

## This helps create a table that checks if a sampled parameter was between two points
CI_table_helper <- function(chain){
  truval <- as.numeric(chain[length(chain)])
  param_name <- chain[1]
  chain <- as.numeric(chain[-c(1,length(chain))])
  lwr <- quantile(chain,.025)
  upr <- quantile(chain,.975)
  wthin <- sum(lwr < truval & upr > truval)
  cbind(param_name,truval,lwr,upr,wthin)
}

# Function to plot posterior marginals of each parameter
plot_parameter <- function(param, param_name, beg, n, true_value=NA) {
  plot(param[beg:n], type = 'l', main = param_name, xlab = "Index", ylab = "Value")
  if(!is.na(true_value)){
    abline(h = true_value, col = 'red')
  }
  plot(density(param[beg:n]), main = paste("Density of", param_name),
       xlab = paste("Value-True Val",ifelse(!is.na(true_value),true_value,'none')), ylab = "Density")
  
  if(!is.na(true_value)){
    abline(v = true_value, col = 'red')
  }
}

# Function to make all values < 0 to equal min of col-used in data cleaning
fix_0s <- function(data_set){
  for(i in 1:ncol(data_set)){
    minposval <- min(data_set[data_set[,i]>0,i])
    data_set[data_set[,i]<0,i] <- minposval
  }
  data_set
}

## A more stable way of inverting a matrix that will then multiply another
better_solve <- function(A,B,C){
  
  (L = t(chol(B)))
  
  Linv_tA = forwardsolve(L, t(A))
  Linv_C = forwardsolve(L, C)
  
  crossprod(Linv_tA, Linv_C)
}

## Used for testing things out
onedimdist <- function(x1, x2){
  abs(x2-x1)
}

## Calcluates Euclidean distance given two columns of coordinates
twodimdist <- function(coords){ ## Takes in two col matrix or datafram of coords
  x1 <- coords[,1]
  y1 <- coords[,2]
  xs <- outer(x1, x1, FUN = function(x1,x2) (x1-x2)^2)
  ys <- outer(y1, y1, FUN = function(y1,y2) (y1-y2)^2)
  dists <- sqrt(xs + ys)
  dists
}

mean_func <- function(x){0}

mean_funct2 <- function(coords){
  rep(0, nrow(coords))
}

###### A few options for exponential functions in gaussian process making
cov_sqrd_exponential <- function(h, range, sigma_squared) {
  result <- sigma_squared * exp(-1 * (h^2) / (range^2))
  return(result)
}

cov_exponential <- function(h, range, sigma_squared) {
  result <- sigma_squared * exp(-1 * (h) / (range))
  return(result)
}

## Creates a covariance matrix given row of point coordinates (one dimensional)
covfunct <- function(t, range = 1, sigma_squared = 1){
  hs <- outer(t, t, FUN = onedimdist)
  
  cov_sqrd_exponential(hs, range, sigma_squared)
}

## Creates a covariance matrix given rows of point coordinates (two dimensional)
covfunct2 <- function(coords, range = .1, tau_squared = 1){
  hs <- twodimdist(coords)
  # cov_exponential(hs, .8, 1)
  cov_sqrd_exponential(hs, range, tau_squared)
  # cov_matern(.1709, .1, hs)
  # cov_matern(c(.1709,.0072),hs)
  # optim(par = c(3,8), fn = cov_matern, d = hs)
}

## This is what we used to create gamma priors
create_gamma_prior <- function(mn, vr){
  alph <- (mn^2)/vr
  bet <- alph/mn
  c(alph, bet)
}


######################################
###### The covariance functions ######
######################################
## Euclidean Distance
###Exponential
Euclid_exponential <- function(d, range, sigma_squared) {
  result <- sigma_squared * exp(-3 * d / range)
  return(result)
}

# Gaussian
Euclid_lin_with_sill <- function(d, range, sigma_squared) {
  result <- sigma_squared * exp(-1*((d/range)**2))
  return(result)
}

# spherical
Euclid_lin_with_sill <- function(d, range, sigma_squared) {
  result <- sigma_squared * (1-3*d/(2*range) + (d**3)/(2*(range**3)))*(d/range <= 1)
  return(result)
}

## Tail Up (h = distance)
###Exponential
tail_up_exponential <- function(tu_flow_con_dist, range, sigma_squared, weight_mat) {
  result <- weight_mat*sigma_squared * exp(-3 * tu_flow_con_dist / range)
  return(result)
}

# linear with sill
tail_up_lin_with_sill <- function(tu_flow_con_dist, range, sigma_squared, weight_mat) {

  result <- weight_mat*sigma_squared * (1-tu_flow_con_dist/range)*(tu_flow_con_dist/range <= 1)
  return(result)
}

# spherical
tail_up_lin_with_sill <- function(tu_flow_con_dist, range, sigma_squared, weight_mat) {
  result <- weight_mat*sigma_squared * (1-3*h/(2*range) + (h**3)/(2*(range**3)))*(h/range <= 1)
  return(result)
}

## Tail Down
###Exponential
tail_down_exponential <- function(tu_flow_con_dist, a, b, range, sigma_squared) {
  # Just get both matrices, and put them through accordingly, then add them up as the covar matrix
  flow_uncon_cov <- sigma_squared * exp(-3 * (a+b) / range)
  flow_con_cov <- sigma_squared * exp(-3 * tu_flow_con_dist / range)
  return(flow_uncon_cov+flow_con_cov)
}


# linear with sill
tail_down_lin_with_sill <- function(tu_flow_con_dist, a, b, range, sigma_squared) {
  ## From what I'm referenceing, the a and b mats are switched, where in the paper
  ## a < b, my matrixes in the exploring code have a > b, so I'm switching them here, 
  ## a is b and b is a. 
  flow_uncon_cov <- sigma_squared * (1-a/range)*(a/range <= 1)
  flow_con_cov <- sigma_sqaured * (1-tu_flow_con_dist/range)*(tu_flow_con_dist/range <= 1)
  result <- flow_uncon_cov+flow_con_cov
  return(result)
}

# spherical
tail_down_spherical <- function(tu_flow_con_dist, a, b, range, sigma_squared) {
  

  flow_uncon_cov <- sigma_squared * (1-3*tu_flow_con_dist/(2*range) + (tu_flow_con_dist**3)/(2*(range**3)))*(tu_flow_con_dist/range <= 1)

  flow_con_cov <- sigma_squared * (1-3*b/(2*range) + (a)/(2*range))*(1-a/range)*(a/range <= 1)

  return(result)
}

####################################################
# The full conditionals for the gibs sampler
####################################################

# See equation 2.18 in paper for where this comes from
beta_FC <- function(y_l, sig_mat_0, sig_mat_1, X_mat, T_mat_l, betas_0, betas_1, thetas_0, thetas_1, l, sig2_time){
  
  
  m_b <- rbind(X_mat[l,]%*%thetas_0+better_solve(t(sig_mat_0[l,-l]),sig_mat_0[-l,-l], betas_0[-l]-X_mat[-l,]%*%thetas_0),
             X_mat[l,]%*%thetas_1+better_solve(t(sig_mat_1[l,-l]),sig_mat_1[-l,-l], betas_1[-l]-X_mat[-l,]%*%thetas_1))
  
  D_b <- rbind(c(sig_mat_0[l,l]+better_solve(t(sig_mat_0[l,-l]),sig_mat_0[-l,-l],sig_mat_0[-l,l]),0),
               c(0,sig_mat_1[l,l]+better_solve(t(sig_mat_1[l,-l]),sig_mat_1[-l,-l],sig_mat_1[-l,l])))
  
  mn <- solve((t(T_mat_l)%*%T_mat_l)/sig2_time[l]+solve(D_b),(t(T_mat_l)%*%y_l)/sig2_time[l]+solve(D_b,m_b))
  
  var <- solve((t(T_mat_l)%*%T_mat_l)/sig2_time[l]+solve(D_b))
  bet_samp <- rmvnorm(1,mn,var)
}

# See equation 2.19  complete conditional distribution
sig2_time_FC <- function(y_l, T_mat_l, alph_sig2_time, bet_sig2_time,beta_0,beta_1){
  Bs <- rbind(beta_0,beta_1)
  n <- nrow(T_mat_l)
  alph <- n/2+alph_sig2_time
  betgam <- (t(y_l-T_mat_l%*%Bs)%*%(y_l-T_mat_l%*%Bs))/2 + bet_sig2_time
  new_sig2 <- rinvgamma(1,alph,betgam)
  new_sig2 
}

# I don't think this is used
onedimdist <- function(x1, x2){
  abs(x2-x1)
}

# The full conditional of theta
theta_FC <- function(sig_mat, X, B, V_inv, m, n_iter, p){ # Sig_inv needs to be calced
  sapply(1:n_iter, function(x) solve(better_solve(t(X),sig_mat,X)+V_inv, better_solve(t(X), sig_mat, B)-V_inv%*%m)+t(chol(solve(better_solve(t(X),sig_mat,X)+V_inv)))%*%rnorm(p))
}

# lines(density(beta_FC(V_inv, m, 1000)),xlim = c(.5,1.5))
# 
log_sigtu_FC <- function(sig2, range_tu,
                                           alpha_sig2_tu, beta_sig2_tu,
                                           coords, data = as.vector(y),
                                           covfunctTU, cov_eu_calcd, cov_td_calcd, cov_calcd,
                                           tu_flow_con_dist, net_dim_list, weight_matrix){
  dgamma(sig2, alpha_sig2_tu, beta_sig2_tu, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_eu_calcd +
              cov_td_calcd +
              cov_tu(net_dim_list, tu_flow_con_dist, range_tu, sig2, weight_matrix, covfunctTU), log = TRUE)
}

log_sigtd_FC <- function(sig2, range_td, a_mat, b_mat,
                            alpha_sig2_td, beta_sig2_td,
                            coords, data = as.vector(y),
                            covfunctTD, cov_eu_calcd, cov_tu_calcd, cov_calcd,
                            tu_flow_con_dist, net_dim_list){
  dgamma(sig2, alpha_sig2_td, beta_sig2_td, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_eu_calcd +
              cov_tu_calcd +
              cov_td(tu_flow_con_dist, a_mat, b_mat, net_dim_list, range_td, sig2, covfunctTD), log = TRUE)
}

log_sig2_eu_FC <- function(sig2, range_eu,
                          alpha_sig2_eu, beta_sig2_eu,
                          coords, data = as.vector(y),
                          covfunctEU, cov_tu_calcd, cov_td_calcd, cov_calcd,
                          eu_dist_mat){
  dgamma(sig2, alpha_sig2_eu, beta_sig2_eu, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_tu_calcd +
              cov_td_calcd +
              covfunctEU(eu_dist_mat, range_eu, sig2), log = TRUE)
}

## Don't calc any matrix that isn't necessary! In this case that's all the matrices beside the sig*identity mat
log_sig2_FC <- function(sig2, alph_sig2, bet_sig2,
                        coords, data,
                        cov_tu_calcd, cov_td_calcd, cov_eu_calcd){
  
  dgamma(sig2, alph_sig2, bet_sig2, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), diag(x=sig2,nrow = length(data)) + 
              cov_eu_calcd +
              cov_tu_calcd +
              cov_td_calcd, log = TRUE)
}



log_tau2_FC <- function(tau2, range_eu,
                        alph_tau2, bet_tau2,
                        coords, data = as.vector(y),
                        covfunctEU, cov_tu_calcd, cov_td_calcd, cov_calcd,
                        eu_dist_mat){
  dgamma(tau2, alph_tau2, bet_tau2, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_tu_calcd +
              cov_td_calcd +
              covfunctEU(eu_dist_mat, range_eu, tau2), log = TRUE)
}

log_rangetu_FC <- function(range_tu, sig2,
                          alpha_range_tu, beta_range_tu,
                          coords, data = as.vector(y),
                          covfunctTU, cov_eu_calcd, cov_td_calcd, cov_calcd,
                          tu_flow_con_dist, net_dim_list, weight_matrix){
  dgamma(range_tu, alpha_range_tu, beta_range_tu, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_eu_calcd +
              cov_td_calcd +
              cov_tu(net_dim_list, tu_flow_con_dist, range_tu, sig2, weight_matrix, covfunctTU), log = TRUE)
}

log_rangetd_FC <- function(range_td, sig2, a_mat, b_mat,
                           alpha_range_td, beta_range_td,
                           coords, data = as.vector(y),
                           covfunctTD, cov_eu_calcd, cov_tu_calcd, cov_calcd,
                           tu_flow_con_dist, net_dim_list){
  dgamma(range_td, alpha_range_td, beta_range_td, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_eu_calcd +
              cov_tu_calcd +
              cov_td(tu_flow_con_dist, a_mat, b_mat, net_dim_list, range_td, sig2, covfunctTD), log = TRUE)
  
}

log_range_eu_FC <- function(range_eu, tau2,
                            alpha_range_eu, beta_range_eu,
                            coords, data = as.vector(y),
                            covfunctEU, cov_tu_calcd, cov_td_calcd, cov_calcd,
                            eu_dist_mat){
  dgamma(range_eu, alpha_range_eu, beta_range_eu, log = TRUE)+
    dmvnorm(t(data), t(mean_funct2(coords)), cov_calcd + 
              cov_tu_calcd +
              cov_td_calcd +
              covfunctEU(eu_dist_mat, range_eu, tau2), log = TRUE)
}

## This puts each indiv network through the covariance function so dont have to calc big mat
## 
cov_eu <- function(eu_dist_mat, range_eu, tau2, cov_funct_eu){
  cov_mat_eu <- cov_funct_eu(eu_dist_mat, range_eu, tau2)
  cov_mat_eu
}

cov_td <- function(tu_flow_con_dist, a_mat, b_mat, dim_list, ran_td, sig2_td, cov_funct_td){
  cov_mat_td <- tu_flow_con_dist
  cov_mat_td[1:dim(tu_flow_con_dist)[2],1:dim(tu_flow_con_dist)[2]] <- 0
  for(i in 1:8){
    # Dim list comes from loading in these data
    ind1 <- dim_list[[i]][1]
    ind2 <- dim_list[[i]][2]
    
    cov_mat_td[ind1:ind2,ind1:ind2] <-  cov_funct_td(tu_flow_con_dist[ind1:ind2,ind1:ind2],
                                                     a_mat[ind1:ind2,ind1:ind2], b_mat[ind1:ind2,ind1:ind2],
                                                     ran_td, sig2_td)
    
  }
  cov_mat_td
}

cov_tu <- function(dim_list, tu_flow_con_dist, ran_tu, sig2_tu, weightrix, cov_funct_tu){
  cov_mat_tu <- weightrix
  cov_mat_tu[1:dim(weightrix)[2],1:dim(weightrix)[2]] <- 0
  for(i in 1:8){
    # Dim list comes from loading in these data
    ind1 <- dim_list[[i]][1]
    ind2 <- dim_list[[i]][2]
    cov_mat_tu[ind1:ind2,ind1:ind2] <- cov_funct_tu(tu_flow_con_dist[ind1:ind2,ind1:ind2],
                                                    range = ran_tu,
                                                    sigma_squared = sig2_tu,
                                                    weightrix[ind1:ind2,ind1:ind2])
  }
  cov_mat_tu
}

# llm <- sapply(1:nrow(e), function(i) {
#   apply(samples, 1, function(x) {
#     log_likelyhood_i(i, x[1],x[2],x[3])
#   })
# })


# S = number of post draws
# T_mat_all contains the 2000+ rows of 1s and time values in 2nd col
# Bets_0 and _1 just the 52XS mats of post samps
# sigs2_time also just 52XS of post samps
log_likelihood_i <- function(y_wit_id,S,T_mat_all,bets_0,bets_1,sigs2_time,site_ids){
  # This is the SXN matrix for loo function we are creating here
  
  res_mat <- matrix(NA,nrow = S, ncol = nrow(y_wit_id))
  
  for(s in 1:S){
    nl=1
    for(l in 1:length(site_ids)){
      bets <- rbind(bets_0[l,s],bets_1[l,s])
      sig2 <- sigs2_time[l,s]
      T_mat <- T_mat_all[y_wit_id[,2]==site_ids[l],]
      y <- y_wit_id[y_wit_id[,2]==site_ids[l],1]
      for(n in 1:length(y)){ # Want y_wit_id to have 2 cols, 1st=y, 2nd=site_id
        res_mat[s,nl] <- dnorm(as.numeric(y[n]),T_mat[n,]%*%bets,sqrt(sig2),log=TRUE)
        nl <- nl+1
      }
    }
  }
  return(res_mat)
}

HierSampler <- function(sig_mat, X_mat, y, current_state, m, V_inv, p, cov_to_calc,
                        cov_tu_calcd,
                        cov_td_calcd,
                        cov_eu_calcd,
                        alph_sig2, bet_sig2, 
                        alpha_range_eu, beta_range_eu,
                        alph_sig2_eu, bet_sig2_eu,
                        alph_sig2_tu, bet_sig2_tu,
                        alph_range_tu, bet_range_tu,
                        alph_range_td, bet_range_td,
                        alph_sig2_td, bet_sig2_td,
                        alph_sig2_time, bet_sig2_time
                        
){
  
  
  
  
  current_state$thetas <- theta_FC(sig_mat, X_mat, y, V_inv, m, 1, p) 
  y_minXB <- y - X_mat%*%current_state$thetas
  # y_minXB <- y - X_mat%*%bet_true
  
  ## ONLY NEED cov matrices relevent to each full conditional, so calc some here
  #
  current_state$sig2 <- slice_sampler(10, c(.00001,10000), log_sig2_FC, current_state$sig2, c(alph_sig2,
                                                                                              bet_sig2,
                                                                                              list(as.data.frame(coords)),
                                                                                              list((y_minXB)),
                                                                                              list(cov_tu_calcd),
                                                                                              list(cov_td_calcd),
                                                                                              list(cov_eu_calcd)), movement_num = 10)
  
  cov_calcd <- diag(current_state$sig2, nrow = nrow(rand_data))
  
  ####################################################################
  
  if('e' %in% cov_to_calc){
  current_state$ran_eu <- slice_sampler(10, c(.000001,10000), log_range_eu_FC, current_state$ran_eu, c(current_state$sig2_eu, alpha_range_eu, beta_range_eu,
                                                                                                       list(as.data.frame(coords)), list(y_minXB),
                                                                                                       Euclid_exponential, list(cov_tu_calcd), list(cov_td_calcd), list(cov_eu_calcd),
                                                                                                       list(euc_dist_mat)), movement_num = 100)
  
  cov_eu_calcd <- cov_eu(euc_dist_mat, current_state$ran_eu, current_state$sig2_eu, Euclid_exponential)

  
  
  current_state$sig2_eu <- slice_sampler(2000, c(.000001,2000000), log_sig2_eu_FC, current_state$sig2_eu, c(current_state$ran_eu, alph_sig2_eu, bet_sig2_eu,
                                                                                                            list(as.data.frame(coords)), list(y_minXB),
                                                                                                            Euclid_exponential, list(cov_tu_calcd), list(cov_td_calcd), list(cov_eu_calcd),
                                                                                                            list(euc_dist_mat)), movement_num = 20000)
  cov_eu_calcd <- cov_eu(euc_dist_mat, current_state$ran_eu, current_state$sig2_eu, Euclid_exponential)
  }else if(sum(cov_eu_calcd^2) != 0){
    stop('You specified no Euc distance, but matrix for it != 0!')
  }
  
  if('u' %in% cov_to_calc){
  current_state$ran_tu <- slice_sampler(600, c(.000001,150000), log_rangetu_FC, current_state$ran_tu, c(current_state$sig2_tu, alph_range_tu, bet_range_tu,
                                                                                                        list(as.data.frame(coords)), list(y_minXB),
                                                                                                        tail_up_exponential, list(cov_eu_calcd), list(cov_td_calcd), list(cov_calcd),
                                                                                                        list(tu_flow_con_dist), list(dim_list), list(weightrix)), movement_num = 1000)
  cov_tu_calcd <- cov_tu(dim_list, tu_flow_con_dist, current_state$ran_tu, current_state$sig2_tu, weightrix, tail_up_exponential)
  
  
  current_state$sig2_tu <- slice_sampler(50, c(.000001,150000), log_sigtu_FC, current_state$sig2_tu, c(current_state$ran_tu, alph_sig2_tu, bet_sig2_tu,
                                                                                                       list(as.data.frame(coords)), list(y_minXB),
                                                                                                       tail_up_exponential, list(cov_eu_calcd), list(cov_td_calcd), list(cov_calcd),
                                                                                                       list(tu_flow_con_dist), list(dim_list), list(weightrix)), movement_num = 50)
  cov_tu_calcd <- cov_tu(dim_list, tu_flow_con_dist, current_state$ran_tu, current_state$sig2_tu, weightrix, tail_up_exponential)
  }else if(sum(cov_tu_calcd^2) != 0){
    stop('You specified no TU distance, but matrix for it != 0!')
  }
  
  if('d' %in% cov_to_calc){
  current_state$ran_td <- slice_sampler(600, c(.000001,150000), log_rangetd_FC, current_state$ran_td, c(current_state$sig2_td, list(a_mat_uncon), list(b_mat_uncon), alph_range_td, bet_range_td,
                                                                                                        list(as.data.frame(coords)), list(y_minXB),
                                                                                                        tail_down_exponential, list(cov_eu_calcd), list(cov_tu_calcd), list(cov_calcd),
                                                                                                        list(tu_flow_con_dist), list(dim_list)), movement_num = 1000)
  
  cov_td_calcd <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, current_state$ran_td, current_state$sig2_td, tail_down_exponential)
  
  current_state$sig2_td <- slice_sampler(600, c(.000001,150000), log_sigtd_FC, current_state$sig2_td, c(current_state$ran_td, list(a_mat_uncon), list(b_mat_uncon), alph_sig2_td, bet_sig2_td,
                                                                                                        list(as.data.frame(coords)), list(y_minXB),
                                                                                                        tail_down_exponential, list(cov_eu_calcd), list(cov_tu_calcd), list(cov_calcd),
                                                                                                        list(tu_flow_con_dist), list(dim_list)), movement_num = 1000)
  
  cov_td_calcd <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, current_state$ran_td, current_state$sig2_td, tail_down_exponential) 
  }else if(sum(cov_td_calcd^2) != 0){
    stop('You specified no TD distance, but matrix for it != 0!')
  }
  
  
  return(current_state)
}

Sampler <- function(n, y, T_mat,site_ids, X_mat, current_state,samp_data,cov_to_calc,
                    # The hyperparameters
                    m,V_inv,
                    alph_sig2_tu,
                    alph_sig2_td,
                    bet_sig2_tu,
                    bet_sig2_td,
                    alph_range_tu,
                    alph_range_td,
                    bet_range_tu,
                    bet_range_td,
                    alph_sig2_eu,
                    bet_sig2_eu,
                    alpha_range_eu,
                    beta_range_eu,
                    alph_sig2,
                    bet_sig2,
                    alph_sig2_time,
                    bet_sig2_time
                    ## Now the 
){
  time_begin <- Sys.time()
  q <- nrow(X_mat)
  p <- ncol(X_mat)
  ## Places to store samples
  bets_0 <- matrix(0, nrow = q, ncol = n)
  bets_1 <- matrix(0, nrow = q, ncol = n)
  sigs2_time <- matrix(0, nrow = q, ncol = n)
  
  current_state <- list(sig2_time = sig2_time, sig2_1= sig2_1, sig2_eu_1 = sig2_eu_1, sig2_tu_1 = sig2_tu_1, sig2_td_1 = sig2_td_1, ran_tu_1 = ran_tu_1, ran_td_1 = ran_td_1, ran_eu_1 = ran_eu_1, betas_1 = betas_1, thetas_1 = thetas_1,
                        sig2_0 = sig2_0, sig2_eu_0 = sig2_eu_0, sig2_tu_0 = sig2_tu_0, sig2_td_0 = sig2_td_0, ran_tu_0 = ran_tu_0, ran_td_0 = ran_td_0, ran_eu_0 = ran_eu_0, betas_0 = betas_0, thetas_0 = thetas_0)
  
  
  ## For hierarchical
  sigs2_1 <- rep(NA,n)
  sigs2_eu_1 <- rep(NA,n)
  sigs2_tu_1 <- rep(NA,n)
  sigs2_td_1 <- rep(NA,n)
  rans_tu_1 <- rep(NA,n)
  rans_td_1 <- rep(NA,n)
  rans_eu_1 <- rep(NA,n)
  thets_1 <- matrix(NA,nrow=p,ncol=n)
  
  sigs2_0 <- rep(NA,n)
  sigs2_eu_0 <- rep(NA,n)
  sigs2_tu_0 <- rep(NA,n)
  sigs2_td_0 <- rep(NA,n)
  rans_tu_0 <- rep(NA,n)
  rans_td_0 <- rep(NA,n)
  rans_eu_0 <- rep(NA,n)
  thets_0 <- matrix(NA,nrow=p,ncol=n)
  
  for(j in 1:n){
    
    
    for(i in 1:length(current_state$betas_1)){
      
      l <- i
      y_l <- as.matrix(y[samp_data$site_id==site_ids[i]])
      T_mat_l <- as.matrix(T_mat[samp_data$site_id==site_ids[i],])
      
      bets_01 <- beta_FC(y_l,sig_mat_0,sig_mat_1,X_mat,T_mat_l,current_state$betas_0,current_state$betas_1,
                         current_state$thetas_0, current_state$thetas_1,l,current_state$sig2_time)
      bets_0[i,j] <- bets_01[1]
      bets_1[i,j] <- bets_01[2]
      
      ## Now to sample the sid2_time
      # bets_01 <- c(B0[i],B1[i])
      sigs2_time[i,j] <- sig2_time_FC(y_l,T_mat_l,alph_sig2_time,bet_sig2_time,bets_01[1],bets_01[2])
    }
    current_state$betas_0 <- bets_0[,j]
    current_state$betas_1 <- bets_1[,j]
    current_state$sig2_time <- sigs2_time[,j]
    
    ## Now for an iteration of the hierarchical model
    # First for the intercept
    
  
    if('u' %in% cov_to_calc){
      cov_tu_calcd_0 <- cov_tu(dim_list, flow_con_mat, current_state$ran_tu_0, current_state$sig2_tu_0, weightrix, tail_up_exponential)
    } else{
      cov_tu_calcd_0 = 0
    }
    if('d' %in% cov_to_calc){
      cov_td_calcd_0 <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, current_state$ran_td_0, current_state$sig2_td_0, tail_down_exponential)
    } else{
      cov_td_calcd_0 = 0
    }
    if('e' %in% cov_to_calc){
      cov_eu_calcd_0 <- cov_eu(euc_dist_mat, current_state$ran_eu_0, current_state$sig2_eu_0, Euclid_exponential)
    } else{
      cov_eu_calcd_0 = 0
    }
    cov_calcd_0 <- diag(current_state$sig2_0, nrow = nrow(rand_data))
    
    ## Make sure matrix only gives what we want
    
    sig_mat_0 <- cov_td_calcd_0 + cov_tu_calcd_0 + cov_eu_calcd_0 + cov_calcd_0
    
    current_state_0 <- current_state[sapply(names(current_state), function(x) grepl('0', x))]
    
    names(current_state_0) <- c('sig2','sig2_eu','sig2_tu','sig2_td','ran_tu','ran_td','ran_eu','betas','thetas')
    current_state_0 <- HierSampler(sig_mat_0,X_mat,B0,current_state_0, m, V_inv,p,cov_to_calc,
                                   cov_tu_calcd_0,
                                   cov_td_calcd_0,
                                   cov_eu_calcd_0,
                                   alph_sig2, bet_sig2, 
                                   alpha_range_eu, beta_range_eu,
                                   alph_sig2_eu, bet_sig2_eu,
                                   alph_sig2_tu, bet_sig2_tu,
                                   alph_range_tu, bet_range_tu,
                                   alph_range_td, bet_range_td,
                                   alph_sig2_td, bet_sig2_td,
                                   alph_sig2_time, bet_sig2_time)
    
    thets_0[,j] <- current_state_0$thetas
    sigs2_0[j] <- current_state_0$sig2
    sigs2_eu_0[j] <- current_state_0$sig2_eu
    rans_eu_0[j] <- current_state_0$ran_eu
    rans_tu_0[j] <- current_state_0$ran_tu
    rans_td_0[j] <- current_state_0$ran_td
    sigs2_tu_0[j] <- current_state_0$sig2_tu
    sigs2_td_0[j] <- current_state_0$sig2_td
    
    temp <- which(sapply(names(current_state), function(x) grepl('0', x)))
    current_state[names(current_state)[temp]] <- current_state_0
    #then for the slopes
    
    current_state_1 <- current_state[sapply(names(current_state), function(x) grepl('1', x))]
    
    if('u' %in% cov_to_calc){
      cov_tu_calcd_1 <- cov_tu(dim_list, flow_con_mat, current_state$ran_tu_1, current_state$sig2_tu_1, weightrix, tail_up_exponential)
    } else{
      cov_tu_calcd_1 = 0
    }
    if('d' %in% cov_to_calc){
      cov_td_calcd_1 <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, current_state$ran_td_1, current_state$sig2_td_1, tail_down_exponential)
    } else{
      cov_td_calcd_1 = 0
    }
    if('e' %in% cov_to_calc){
      cov_eu_calcd_1 <- cov_eu(euc_dist_mat, current_state$ran_eu_1, current_state$sig2_eu_1, Euclid_exponential)
    } else{
      cov_eu_calcd_1 = 0
    }
    cov_calcd_1 <- diag(current_state$sig2_1, nrow = nrow(rand_data))
    
    ## Make sure matrix only gives what we want
    
    sig_mat_1 <- cov_td_calcd_1 + cov_tu_calcd_1 + cov_eu_calcd_1 + cov_calcd_1
    
    
    names(current_state_1) <- c('sig2','sig2_eu','sig2_tu','sig2_td','ran_tu','ran_td','ran_eu','betas','thetas')
    current_state_1 <- HierSampler(sig_mat_1,X_mat,B1,current_state_1, m, V_inv,p,cov_to_calc,
                                   cov_tu_calcd_1,
                                   cov_td_calcd_1,
                                   cov_eu_calcd_1,
                                   alph_sig2, bet_sig2, 
                                   alpha_range_eu, beta_range_eu,
                                   alph_sig2_eu, bet_sig2_eu,
                                   alph_sig2_tu, bet_sig2_tu,
                                   alph_range_tu, bet_range_tu,
                                   alph_range_td, bet_range_td,
                                   alph_sig2_td, bet_sig2_td,
                                   alph_sig2_time, bet_sig2_time)
    
    thets_1[,j] <- current_state_1$thetas
    sigs2_1[j] <- current_state_1$sig2
    sigs2_eu_1[j] <- current_state_1$sig2_eu
    rans_eu_1[j] <- current_state_1$ran_eu
    rans_tu_1[j] <- current_state_1$ran_tu
    rans_td_1[j] <- current_state_1$ran_td
    sigs2_tu_1[j] <- current_state_1$sig2_tu
    sigs2_td_1[j] <- current_state_1$sig2_td
    
    temp <- which(sapply(names(current_state), function(x) grepl('1', x)))
    current_state[names(current_state)[temp]] <- current_state_1
    
  }
  
  time_end <- Sys.time()
  print(time_end - time_begin)
  
  return(list(bets_0 = bets_0, bets_1 = bets_1, sigs2_time = sigs2_time,
              thets_0 = thets_0, thets_1 = thets_1, sigs2_0 = sigs2_0, 
              sigs2_1 = sigs2_1, rans_eu_0 = rans_eu_0, rans_eu_1 = rans_eu_1,
              rans_tu_0 = rans_tu_0, rans_tu_1 = rans_tu_1, rans_td_0 = rans_td_0,
              rans_td_1 = rans_td_1, sigs2_tu_0 = sigs2_tu_0, 
              sigs2_tu_1 = sigs2_tu_1, sigs2_td_0 = sigs2_td_0, 
              sigs2_td_1 = sigs2_td_1, sigs2_eu_0 = sigs2_eu_0, sigs2_eu_1 = sigs2_eu_1))
}



impute_missing_em <- function(data, max_iter = 10000, tol = 1e-7) {
  n <- nrow(data)
  smllnum <- tol
  tries <- max_iter
  
  # Identify rows with missing values
  ind1 <- which(apply(data, 1, function(row) any(is.na(row))))
  new_oli <- data[ind1, ]
  
  # Initialize old predictions (x_olds) with small non-zero values
  x_olds <- list()
  for (i in 1:nrow(new_oli)) {
    rw1 <- new_oli[i, ]
    ind <- which(is.na(rw1))
    x_pred <- rep(1e-5, length(ind))
    x_olds <- append(x_olds, x_pred)
  }
  
  num <- 0
  repeat {
    x_news <- list()
    
    # Step 1: Expectation - Compute column means and covariance matrix
    means <- colMeans(data, na.rm = TRUE)
    covmat <- cov(data, use = 'pairwise.complete.obs')
    
    # Impute missing values row by row
    for (i in 1:nrow(new_oli)) {
      rw1 <- new_oli[i, ]
      ind <- which(is.na(rw1))
      x_pred <- means[ind] + covmat[ind, -ind] %*% 
        solve(covmat[-ind, -ind], rw1[-ind] - means[-ind])
      
      # Store new predictions
      x_news <- append(x_news, x_pred)
      
      # Update the original data with imputed values
      data[ind1[i], ind] <- x_pred
    }
    
    # Check convergence by comparing new and old predictions
    diff_ratio <- max((unlist(x_news) - unlist(x_olds)) / unlist(x_olds))
    
    if (diff_ratio <= smllnum) {
      break  # Converged
    }
    
    if (num >= tries) {
      message("Algorithm did not converge after ", tries, " iterations.")
      break
    }
    
    # Update x_olds for the next iteration
    x_olds <- x_news
    num <- num + 1
  }
  
  return(list(data,ind1))
}
create_diag_table <- function(bets_0, bets_1, sigs2_time, thets_0, thets_1, 
                              sigs2_0, sigs2_1, rans_eu_0=NA, rans_eu_1=NA, rans_tu_0=NA,
                              rans_tu_1=NA, rans_td_0=NA, rans_td_1=NA, sigs2_eu_0=NA, 
                              sigs2_eu_1=NA, sigs2_tu_0=NA, sigs2_tu_1=NA, 
                              sigs2_td_0=NA, sigs2_td_1=NA,cov_to_calc){

df <- data.frame(rbind(t(apply(cbind(rep('beta_0',length(bets_0[,1])),bets_0,B0),1,MCMCdiag)),
                 t(apply(cbind(rep('beta_1',length(bets_1[,1])),bets_1,B1),1,MCMCdiag)),
                 t(apply(cbind(rep('sig2_time',length(sigs2_time[,1])),sigs2_time,sig2_time_true),1,MCMCdiag)),
                 t(apply(cbind(rep('thets',length(thets_0[,1])),thets_0,rep(1,length(thets_0[,1]))),1,MCMCdiag)),
                 t(apply(cbind(rep('thets',length(thets_1[,1])),thets_1,rep(1,length(thets_0[,1]))),1,MCMCdiag)),
                 MCMCdiag(c('sigs2_npspat',sigs2_0,sig2_true_B0)),
                 MCMCdiag(c('sigs2_nospat',sigs2_1,sig2_true_B1)),
                 if('e' %in% cov_to_calc){MCMCdiag(c('rans',rans_eu_0,ran_eu_true_B0))},
                 if('e' %in% cov_to_calc){MCMCdiag(c('rans',rans_eu_1,ran_eu_true_B1))},
                 if('u' %in% cov_to_calc){MCMCdiag(c('rans',rans_tu_0,ran_tu_true_B0))},
                 if('u' %in% cov_to_calc){MCMCdiag(c('rans',rans_tu_1,ran_tu_true_B1))},
                 if('d' %in% cov_to_calc){MCMCdiag(c('rans',rans_td_0,ran_td_true_B0))},
                 if('d' %in% cov_to_calc){MCMCdiag(c('rans',rans_td_1,ran_td_true_B1))},
                 if('e' %in% cov_to_calc){MCMCdiag(c('sigs2',sigs2_eu_0,sig2_eu_true_B0))},
                 if('e' %in% cov_to_calc){MCMCdiag(c('sigs2',sigs2_eu_1,sig2_eu_true_B1))},
                 if('u' %in% cov_to_calc){MCMCdiag(c('sigs2',sigs2_tu_0,sig2_tu_true_B0))},
                 if('u' %in% cov_to_calc){MCMCdiag(c('sigs2',sigs2_tu_1,sig2_tu_true_B1))},
                 if('d' %in% cov_to_calc){MCMCdiag(c('sigs2',sigs2_td_0,sig2_td_true_B0))},
                 if('d' %in% cov_to_calc){MCMCdiag(c('sigs2',sigs2_td_1,sig2_td_true_B1))})

)
  return(df)
}

get_diag_table <- function(res,cov_to_calc){
  bets_0 <- res$bets_0
  # hist(rowMeans(bets_0)[-23])
  bets_1 <- res$bets_1
  # hist(rowMeans(bets_1)[-23])
  sigs2_time <- res$sigs2_time
  
  thets_0 <- res$thets_0
  thets_1 <- res$thets_1
  sigs2_0 <- res$sigs2_0
  sigs2_1 <- res$sigs2_1
  rans_eu_0 <- res$rans_eu_0
  rans_eu_1 <- res$rans_eu_1
  rans_tu_0 <- res$rans_tu_0
  rans_tu_1 <- res$rans_tu_1
  rans_td_0 <- res$rans_td_0
  rans_td_1 <- res$rans_td_1
  
  sigs2_tu_0 <- res$sigs2_tu_0
  sigs2_tu_1 <- res$sigs2_tu_1
  sigs2_td_0 <- res$sigs2_td_0
  sigs2_td_1 <- res$sigs2_td_1
  sigs2_eu_0 <- res$sigs2_eu_0
  sigs2_eu_1 <- res$sigs2_eu_1
  
  diags_table <- create_diag_table(bets_0, bets_1, sigs2_time, thets_0, thets_1,
                                   sigs2_0, sigs2_1, rans_eu_0, rans_eu_1,
                                   rans_tu_0, rans_tu_1, rans_td_0, rans_td_1,
                                   sigs2_eu_0, sigs2_eu_1, sigs2_tu_0, sigs2_tu_1,
                                   sigs2_td_0, sigs2_td_1,cov_to_calc)
  
}

loo_funct <- function(res){
  bets_0 <- res$bets_0[,1001:ncol(res$bets_0)]
  bets_1 <- res$bets_1[,1001:ncol(res$bets_1)]
  sigs2_time <- res$sigs2_time[,1001:ncol(res$sigs2_time)]
  S <- ncol(bets_0)
  sbynmat <- log_likelihood_i(y_wit_id, S, T_mat, bets_0, bets_1, sigs2_time, site_ids)
  return(loo(sbynmat,k_threshold=0.7))
}

## Time to plot and see if models make sense

# make matrix of B0+XB1 for grid of time points, and for every sample for two locations

# the grid


## Function to create visual
lin_vis <- function(loc){
  ys <- time_dat$DOC_mg_L_SCAN[time_dat$site_id==site_ids[loc]]
  tim <- time_dat$datetime[time_dat$site_id==site_ids[loc]]
  timez <- seq(min(time_dat$datetime),max(time_dat$datetime),length.out = 100)
  predictions <- sapply(timez, function(t) {
    bets_0[loc,] + t * bets_1[loc,]  # Compute B0 + time * B1 for each sample
  })
  upr <- apply(predictions,2,function(x) quantile(x,.975))
  lwr <- apply(predictions,2,function(x) quantile(x,.025))
  mean_line <- apply(predictions,2,mean)
  ols <- lm(log(ys)~tim)$coef
  mean_ols <- ols[1]+timez*ols[2]
  
  # now plot all these
  plot(tim,log(ys),main=paste('Site:',site_ids[loc]),xlab='Time', ylab='log(DOC)',ylim=c(min(lwr),max(upr)))
  lines(timez,mean_line)
  lines(timez,upr,lty=2)
  lines(timez,lwr,lty=2)
  lines(timez,mean_ols,col='red',lty=2)
  
}
