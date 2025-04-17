
##########################################################################
# Now... for what we've all been waiting for... RIVER DATAAAAA
##########################################################################
## Load in the data, the distances, and explain functions from SSN that might be handy
# library('devtools')
# devtools::install_local("C:/Users/jlindsey/Documents/Grad WorkClasses/Research/SSN_1.1.17.tar.gz")
# ^^^ How to install package from local file
source("slicy_dicey.R")
source("get_weights.R")
source("functions.R")

library('tidyverse')

# Won't work without data
data <- read.csv('static_data.csv')

## IF YOU DON"T HAVE THE DATA: Just load this in
coords <- read.csv('coords.csv')
# data <- data[-23,]


# # C:/Users/jlindsey/Documents/Grad WorkClasses/Research/official proj data/LSN_Megafire.ssn/LSN_Megafire.ssn/distance/obs/dist.net8.RData
obs_dist_list <- list()
dim_list <- list()
ind <- 0
for(i in 1:8){
  file_path <- file.path('./dist_and_weights_mats/obs',
                         paste("dist.net", i, ".RData", sep = ""))
  
  con <- file(file_path, open = "rb")  # open connection
  obs_dist_list[[i]] <- unserialize(con)  # read from connection
  close(con)  # close connection right after use ✅
  
  if(i == 5){
    obs_dist_list[[i]] <- obs_dist_list[[i]][-12, -12]
  }
  
  # Store matrix placement info
  dim_list[[i]] <- c(ind + 1, ind + ncol(obs_dist_list[[i]]))
  ind <- dim_list[[i]][2]
}


####################################################
########## Cleaning Data NOT REALLY THO
####################################################
# colnm <- colnames(mf04p@obspoints@SSNPoints[[1]]@point.data)
# 
# ssn_object <- mf04p
# 
# for (i in 9:86){
#   mf04p@obspoints@SSNPoints[[1]]@point.data[,paste0(colnm[i])] <- as.numeric(mf04p@obspoints@SSNPoints[[1]]@point.data[,paste0(colnm[i])])
# }
# 
# for(i in 4:5){
#   mf04p@obspoints@SSNPoints[[1]]@point.data[,paste0(colnm[i])] <- as.factor(mf04p@obspoints@SSNPoints[[1]]@point.data[,paste0(colnm[i])])
# }
# data <- mf04p@obspoints@SSNPoints[[1]]@point.data

####################################################
# Creating needed distance matrices and weights
####################################################

# This loads in a lot of the needed stuff for this to work
load('res.RData')

## These next matrixes are essential for computing tail-up and down covariance matrices. They come from the weight_mat function that 
# compiles much of the code from the SSN package, but it still requires the SSN object. I save the matrices so you don't need the SSN package
# to run this analysis-I'll save these as CSV files

euc_dist_mat <- twodimdist(coords)
euc_dist_mat <- euc_dist_mat[-23,-23]

weightrix <- res$weight

net_zero <- res$net.zero
net_inf <- net_zero
net_inf[net_inf == 0] <- Inf
net_inf <- net_inf-net_zero
a_mat <- res$a.mat
b_mat <- res$b.mat

# the point taken out here only has 2 time points, so we take it out
weightrix <- weightrix[-23,-23]
a_mat <- a_mat[-23,-23]
b_mat <- b_mat[-23,-23]
coords <- coords[-23,]


### now for tail up/down matrices
dist_junct <- res$dist_from_junct #+ net_inf## Matrix is only good for tail down models, for its exponential function (a+b)
dist_junct <- dist_junct[-23,-23]
tot_hydro_dist <- dist_junct + t(dist_junct)
## Matrix for tail-up (flow connected only)
flow_con_mat <- 1 - (pmin(dist_junct, t(dist_junct))>0)*1

tot_stream_dist <- dist_junct + t(dist_junct)

tu_flow_con_dist <- tot_stream_dist*(flow_con_mat)
a_mat_uncon <- a_mat*(b_mat > 0)
b_mat_uncon <- b_mat*(b_mat > 0)

# Done.

####################################################
# Setting up Priors
####################################################

## Prior hypparams for gamma for sigma^2 for tail up and down for exp cov function
sig2_tutd_pri_mean <- 40
sig2_tutd_pri_var <- 30^2
res <- create_gamma_prior(sig2_tutd_pri_mean, sig2_tutd_pri_var) 
alph_sig2_tu <- alph_sig2_td <- res[1]
bet_sig2_tu <- bet_sig2_td <- res[2]

## Prior hyppparam for gamma for range for exp cov function tu and td
ran_pri_tu_mean <- 1000
ran_pri_tu_var <- 900^2
res <- create_gamma_prior(ran_pri_tu_mean,ran_pri_tu_var)
# library(invgamma)
plot(seq(1:10000),dgamma(seq(1:10000), res[1], res[2]))
res[1]/(res[2]^2)
alph_range_tu <- alph_range_td <- res[1]
bet_range_tu <- bet_range_td <- res[2]

## Prior hypparam for gamma for tau^2 for exp cov function for euclid
sig2_eu_pri_mean <- 80
sig2_eu_pri_var <- 100^2
res <- create_gamma_prior(sig2_eu_pri_mean, sig2_eu_pri_var) 
alph_sig2_eu <- res[1]
bet_sig2_eu <- res[2]

ran_pri_eu_mean <- .05
ran_pri_eu_var <- 2^2
res <- create_gamma_prior(ran_pri_eu_mean, ran_pri_eu_var) ## Prior hypparam for gamma for range for exp cov function for euclid
alpha_range_eu <- res[1]
beta_range_eu <- res[2]

sig2_pri_mean <- 20
sig2_pri_var <- 10^2
res <- create_gamma_prior(sig2_pri_mean, sig2_pri_var) 
alph_sig2 <- res[1]
bet_sig2 <- res[2]


alph_sig2_time <- 3
bet_sig2_time <- 5

####################################################
# Creating the full conditional functions
####################################################
# It's in the functions.R file

####################################################
# Creating simulated data
####################################################

set.seed(24)
full_df <- read.csv('fullDF.csv')
site_info <- read.csv('site_info.csv')
MegafirePointSamples <- read.csv('MegafirePointSamples.csv')
site_info$site_id %in% data$site_id
static_data <- merge(
  full_df[, c('site_id', colnames(full_df)[!(colnames(full_df) %in% colnames(site_info))])], 
  site_info,
  by = 'site_id',
  all.x = TRUE # This ensures we only keep site_ids from site_info
)
# We're going to try to get Y's like this
plot(density(MegafirePointSamples$DOC_mg_L_SCAN, na.rm = TRUE))
static_data$site
# Make beta0s around 5, perhaps with higher altitude

head(static_data)

data <- data[-23,]

## This is where we make the fake static covariates
## This code uses 'data' a lot to get the dimensions necessary for creating some things: So use this 

############################################################
# to act as a substitute if you don't have data

data <- data.frame(X1 = rep(NA,51),
                   X2 = rep(NA,51))


X_mat <- matrix(c(rep(1,nrow(data)), rbeta(nrow(data),2,5), rnorm(nrow(data),50,10),rbinom(nrow(data),1,.5)), ncol = 4)
colnames(X_mat) <- c('interc', 'fake_perc_burn','fake_var', 'no_tree_cover')

## Making the "true" parameters to create the data with, to then recapture
thet_true_B1 <- c(0,-5,.003,-.5)
thet_true_B0 <- c(5,-3,.01,.5)

sig2_tu_true_B1 <- 9
sig2_td_true_B1 <- 4
sig2_eu_true_B1 <- 8 ## This is he variance of the euclidean distance-I shoulda named it differently
sig2_true_B1 <- 1.5
ran_tu_true_B1 <- 1410
ran_td_true_B1 <- 1510
ran_eu_true_B1 <- .1

sig2_tu_true_B0 <- 16
sig2_td_true_B0 <- 36
sig2_eu_true_B0 <- 12 ## This is he variance of the euclidean distance-I shoulda named it differently
sig2_true_B0 <- 14
ran_tu_true_B0 <- 1600
ran_td_true_B0 <- 1200
ran_eu_true_B0 <- .15

## Some priors for the Theta coefficients
m <- rep(0,ncol(X_mat))
V <- diag(100, nrow = ncol(X_mat))
V[1,1] <- 1000
V_inv <- solve(V)
sig2_time_true <- runif(nrow(X_mat),min=2,max=32)


## Creating true covariance matrixes to create fake beta coefficients
cov_tu_calcd_B1 <- cov_tu(dim_list, flow_con_mat, ran_tu_true_B1, sig2_tu_true_B1, weightrix, tail_up_exponential)
cov_td_calcd_B1 <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, ran_td_true_B1, sig2_td_true_B1, tail_down_exponential)
cov_eu_calcd_B1 <- cov_eu(euc_dist_mat, ran_eu_true_B1, sig2_eu_true_B1, Euclid_exponential)

cov_tu_calcd_B0 <- cov_tu(dim_list, flow_con_mat, ran_tu_true_B0, sig2_tu_true_B0, weightrix, tail_up_exponential)
cov_td_calcd_B0 <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, ran_td_true_B0, sig2_td_true_B0, tail_down_exponential)
cov_eu_calcd_B0 <- cov_eu(euc_dist_mat, ran_eu_true_B0, sig2_eu_true_B0, Euclid_exponential)

## Our simulated (I should stop saying fake) beta coefficients
B1 <- X_mat%*%thet_true_B1 + t(chol((cov_tu_calcd_B1 + cov_td_calcd_B1 + cov_eu_calcd_B1 + diag(sig2_true_B0,nrow=nrow(data)))))%*%rnorm(nrow(data),0,1)
B0 <- X_mat%*%thet_true_B0 + t(chol((cov_tu_calcd_B0 + cov_td_calcd_B0 + cov_eu_calcd_B0 + diag(sig2_true_B1,nrow=nrow(data)))))%*%rnorm(nrow(data),0,1)

## Fake time data creation 
T_sim_total <- seq(0.0,6,by=1/365)


# Now, to make a time variable times for each site by sampling
i=0
fake_data <- data.frame(site_id = character(), time = numeric(), fake_DOC = numeric())

#################### This chunk of code was to make sure each location of the fake data
# Had the same amount of observations as the data used. If you don't have the data, use the next chunk
for(site in data$site_id){
  
  n <- sum(!is.na(MegafirePointSamples[MegafirePointSamples$site_id == site,]$DOC_mg_L_SCAN))
  temp_df <- data.frame(site_id=site,time=sample(T_sim_total, n, replace=FALSE))
  i = i+1
  error <- rnorm(length(temp_df$time),0,sqrt(sig2_time_true[i]))
  # error2 <- chol(diag(sig2_time_true[i],nrow=nrow(temp_df)))%*%rnorm(nrow(temp_df),0,1)
  temp_df$fake_DOC <-  B0[i] + temp_df$time*B1[i] + error
  # temp_df$site_ind <- rep(i,nrow(temp_df))
  fake_data <- rbind(fake_data,temp_df)
}

## For no reference data
site_id <- 1:51

for(site in site_id){
  # This randomly creates the amount of observations a location gets, from 20-90 by 10
  n <- sample(c(20,30,40,50,60,70,80,90),1)
  
  # This samples the corresponding time points
  temp_df <- data.frame(site_id=site,time=sample(T_sim_total, n, replace=FALSE))
  i = i+1
  error <- rnorm(length(temp_df$time),0,sqrt(sig2_time_true[i]))
  # error2 <- chol(diag(sig2_time_true[i],nrow=nrow(temp_df)))%*%rnorm(nrow(temp_df),0,1)
  temp_df$fake_DOC <-  B0[i] + temp_df$time*B1[i] + error
  # temp_df$site_ind <- rep(i,nrow(temp_df))
  fake_data <- rbind(fake_data,temp_df)
}



## Simulated time covariates
T_mat <- cbind(rep(1,nrow(fake_data)),fake_data$time)
y <- fake_data$fake_DOC
## The data are made.
### Finally, time to sample this stuff
###########################################################################
############### THE FULL SAMPLE (on fake data) ############################
###########################################################################

# Initial values
sig2_time <- rep(4,nrow(X_mat))
betas_0 <- rep(0,nrow(X_mat))
betas_1 <- rep(0,nrow(X_mat))
p <- ncol(X_mat)
## For hierarchical
sig2_1 <- 20
sig2_eu_1 <- 50
sig2_tu_1 <- 20
sig2_td_1 <- 20
ran_tu_1 <- 1000
ran_td_1 <- 1000
ran_eu_1 <- 5
p <- ncol(X_mat)
thetas_1 <- rep(0,p)

sig2_0 <- 20
sig2_eu_0 <- 50
sig2_tu_0 <- 20
sig2_td_0 <- 20
ran_tu_0 <- 1000
ran_td_0 <- 1000
ran_eu_0 <- 5
thetas_0 <- rep(0,p)


current_state <- list(sig2_time = sig2_time, sig2_1= sig2_1, sig2_eu_1 = sig2_eu_1, sig2_tu_1 = sig2_tu_1, sig2_td_1 = sig2_td_1, ran_tu_1 = ran_tu_1, ran_td_1 = ran_td_1, ran_eu_1 = ran_eu_1, betas_1 = betas_1, thetas_1 = thetas_1,
                      sig2_0 = sig2_0, sig2_eu_0 = sig2_eu_0, sig2_tu_0 = sig2_tu_0, sig2_td_0 = sig2_td_0, ran_tu_0 = ran_tu_0, ran_td_0 = ran_td_0, ran_eu_0 = ran_eu_0, betas_0 = betas_0, thetas_0 = thetas_0)


################## CURRENT STATE FOR TESTING
# current_state <- list(sig2_time = sig2_time, sig2_1= sig2_true_B1, sig2_eu_1 = sig2_eu_true_B1, sig2_tu_1 = sig2_tu_true_B1, sig2_td_1 = sig2_td_true_B1, ran_tu_1 = ran_tu_true_B1, ran_td_1 = ran_td_true_B1, ran_eu_1 = ran_eu_true_B1, betas_1 = betas_1, thetas_1 = thet_true_B1,
#                       sig2_0 = sig2_true_B0, sig2_eu_0 = sig2_eu_true_B0, sig2_tu_0 = sig2_tu_true_B0, sig2_td_0 = sig2_td_true_B0, ran_tu_0 = ran_tu_true_B0, ran_td_0 = ran_td_true_B0, ran_eu_0 = ran_eu_true_B0, betas_0 = betas_0, thetas_0 = thet_true_B0)

##########################

cov_tu_calcd_0 <- cov_tu(dim_list, flow_con_mat, current_state$ran_tu_0, current_state$sig2_tu_0, weightrix, tail_up_exponential)
cov_td_calcd_0 <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, current_state$ran_td_0, current_state$sig2_td_0, tail_down_exponential)
cov_eu_calcd_0 <- cov_eu(euc_dist_mat, current_state$ran_eu_0, current_state$sig2_eu_0, Euclid_exponential)
cov_calcd_0 <- diag(current_state$sig2_0, nrow = nrow(data))

sig_mat_0 <- cov_td_calcd_0 + cov_tu_calcd_0 + cov_eu_calcd_0 + cov_calcd_0

cov_tu_calcd_1 <- cov_tu(dim_list, flow_con_mat, current_state$ran_tu_1, current_state$sig2_tu_1, weightrix, tail_up_exponential)
cov_td_calcd_1 <- cov_td(tu_flow_con_dist, a_mat_uncon, b_mat_uncon, dim_list, current_state$ran_td_1, current_state$sig2_td_1, tail_down_exponential)
cov_eu_calcd_1 <- cov_eu(euc_dist_mat, current_state$ran_eu_1, current_state$sig2_eu_1, Euclid_exponential)
cov_calcd_1 <- diag(current_state$sig2_1, nrow = nrow(data))

sig_mat_1 <- cov_td_calcd_1 + cov_tu_calcd_1 + cov_eu_calcd_1 + cov_calcd_1


#########################

library(mvtnorm)
library(invgamma)
site_ids <- unique(fake_data$site_id)

n=10000
## First, we need to create another function that samples the betas and sig2_time

res_fake <- Sampler(n, y, T_mat,site_ids, X_mat, current_state,fake_data,c('e','u','d'),
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
)
# save(res_fake, file = "simulation_chain.RData")

# Load the saved result correctly
# load("simulation_chain.RData")
beg <- 1000


bets_0 <- res_fake$bets_0[,beg:n]
bets_1 <- res_fake$bets_1[,beg:n]
sigs2_time <- res_fake$sigs2_time[,beg:n]

thets_0 <- res_fake$thets_0[,beg:n]
thets_1 <- res_fake$thets_1[,beg:n]
sigs2_0 <- res_fake$sigs2_0[beg:n]
sigs2_1 <- res_fake$sigs2_1[beg:n]
rans_eu_0 <- res_fake$rans_eu_0[beg:n]
rans_eu_1 <- res_fake$rans_eu_1[beg:n]
rans_tu_0 <- res_fake$rans_tu_0[beg:n]
rans_tu_1 <- res_fake$rans_tu_1[beg:n]
rans_td_0 <- res_fake$rans_td_0[beg:n]
rans_td_1 <- res_fake$rans_td_1[beg:n]

sigs2_tu_0 <- res_fake$sigs2_tu_0[beg:n]
sigs2_tu_1 <- res_fake$sigs2_tu_1[beg:n]
sigs2_td_0 <- res_fake$sigs2_td_0[beg:n]
sigs2_td_1 <- res_fake$sigs2_td_1[beg:n]
sigs2_eu_0 <- res_fake$sigs2_eu_0[beg:n]
sigs2_eu_1 <- res_fake$sigs2_eu_1[beg:n]

## Plot the Thetaaaas!! 
# Lets take a look ourselves real quick
# Plot for `thets_0` and `thets_1`
clnmes <-  c('interc', 'fake_perc_burn','fake_var', 'no_tree_cover')
df_to_plot_0 <- data.frame(meds = apply(thets_0,1,median),
                           upr = apply(thets_0,1,quantile,probs = .975),
                           lwr = apply(thets_0,1,quantile,probs = .025),
                           tru = c(5, -3, 0.01, 0.5),
                           varname = clnmes)

df_to_plot_1 <- data.frame(meds = apply(thets_1,1,median),
                           upr = apply(thets_1,1,quantile,probs = .975),
                           lwr = apply(thets_1,1,quantile,probs = .025),
                           tru = c(0, -5, 0.003, -0.5),
                           varname = clnmes)

# The intercepts:
## Some help from chatGPT for this visual
ggplot(data = df_to_plot_0, aes(x = meds, y = varname)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +  # Credible intervals
  geom_point(aes(color = "Estimate"), size = 3) +  # Estimated medians
  geom_point(aes(x = tru, y = varname, color = "True Value"), shape = 17, size = 3) +  # True values
  scale_color_manual(values = c("Estimate" = "blue", "True Value" = "red")) +  # Define colors
  labs(
    # title = 'Static Variable Coefficients for Intercept',
    x = 'Value',
    y = 'Variable',
    color = "Legend"
  ) +
  xlim(-20, 20) + 
  theme_minimal()

ggplot(data = df_to_plot_1, aes(x = meds, y = varname)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +  # Credible intervals
  geom_point(aes(color = "Estimate"), size = 3) +  # Estimated medians
  geom_point(aes(x = tru, y = varname, color = "True Value"), shape = 17, size = 3) +  # True values
  scale_color_manual(values = c("Estimate" = "blue", "True Value" = "red")) +  # Define colors
  labs(
    # title = 'Static Variable Coefficients for Slope',
    x = 'Value',
    y = 'Variable',
    color = "Legend"
  ) +
  theme_minimal()




## Make a table for latex that shows whether or not parameters were found in a 95% credible interval

#### Below creates a more granular table of results: Too big for paper

# res_df <- as.data.frame(rbind(t(apply(cbind(rep('beta_0',length(bets_0[,1])),bets_0,B0),1,CI_table_helper)),
#       t(apply(cbind(rep('beta_1',length(bets_1[,1])),bets_1,B1),1,CI_table_helper)),
#       t(apply(cbind(rep('sig2_time',length(sigs2_time[,1])),sigs2_time,sig2_time_true),1,CI_table_helper)),
#       t(apply(cbind(rep('thets_0',length(thets_0[,1])),thets_0,thet_true_B0),1,CI_table_helper)),
#       t(apply(cbind(rep('thets_1',length(thets_1[,1])),thets_1,thet_true_B1),1,CI_table_helper)),
#       CI_table_helper(c('sigs2_0',sigs2_0,sig2_true_B0)),
#       CI_table_helper(c('sigs2_1',sigs2_1,sig2_true_B1)),
#       CI_table_helper(c('rans_eu_0',rans_eu_0,ran_eu_true_B0)),
#       CI_table_helper(c('rans_eu_1',rans_eu_1,ran_eu_true_B1)),
#       CI_table_helper(c('rans_tu_0',rans_tu_0,ran_tu_true_B0)),
#       CI_table_helper(c('rans_tu_1',rans_tu_1,ran_tu_true_B1)),
#       CI_table_helper(c('rans_td_0',rans_td_0,ran_td_true_B0)),
#       CI_table_helper(c('rans_td_1',rans_td_1,ran_td_true_B1)),
#       CI_table_helper(c('sigs2_eu_0',sigs2_eu_0,sig2_eu_true_B0)),
#       CI_table_helper(c('sigs2_eu_1',sigs2_eu_1,sig2_eu_true_B1)),
#       CI_table_helper(c('sigs2_tu_0',sigs2_tu_0,sig2_tu_true_B0)),
#       CI_table_helper(c('sigs2_tu_1',sigs2_tu_1,sig2_tu_true_B1)),
#       CI_table_helper(c('sigs2_td_0',sigs2_td_0,sig2_td_true_B0)),
#       CI_table_helper(c('sigs2_td_1',sigs2_td_1,sig2_td_true_B1))))
# 
# ## This here was produced by the help of LLM chatgpt version 4.o
# res_df <- res_df %>%
#   group_by(param_name) %>%
#   mutate(index = row_number()) %>%
#   ungroup() %>%
#   mutate(param_name = paste0(param_name, "_", str_pad(index, 2, pad = "0"))) %>%
#   select(-index)
# 
# res_df$wthin <- as.numeric(res_df$wthin)
# res_df$truval <- as.numeric(res_df$truval)
# res_df$lwr <- as.numeric(res_df$lwr)
# res_df$upr <- as.numeric(res_df$upr)

# a more simplefied table of the results of parameters recaptured

res_df <- as.data.frame(rbind(t(apply(cbind(rep('beta_0',length(bets_0[,1])),bets_0,B0),1,CI_table_helper)),
                              t(apply(cbind(rep('beta_1',length(bets_1[,1])),bets_1,B1),1,CI_table_helper)),
                              t(apply(cbind(rep('sig2_time',length(sigs2_time[,1])),sigs2_time,sig2_time_true),1,CI_table_helper)),
                              t(apply(cbind(rep('thets',length(thets_0[,1])),thets_0,thet_true_B0),1,CI_table_helper)),
                              t(apply(cbind(rep('thets',length(thets_1[,1])),thets_1,thet_true_B1),1,CI_table_helper)),
                              CI_table_helper(c('sigs2_nospat',sigs2_0,sig2_true_B0)),
                              CI_table_helper(c('sigs2_nospat',sigs2_1,sig2_true_B1)),
                              CI_table_helper(c('rans',rans_eu_0,ran_eu_true_B0)),
                              CI_table_helper(c('rans',rans_eu_1,ran_eu_true_B1)),
                              CI_table_helper(c('rans',rans_tu_0,ran_tu_true_B0)),
                              CI_table_helper(c('rans',rans_tu_1,ran_tu_true_B1)),
                              CI_table_helper(c('rans',rans_td_0,ran_td_true_B0)),
                              CI_table_helper(c('rans',rans_td_1,ran_td_true_B1)),
                              CI_table_helper(c('sigs2',sigs2_eu_0,sig2_eu_true_B0)),
                              CI_table_helper(c('sigs2',sigs2_eu_1,sig2_eu_true_B1)),
                              CI_table_helper(c('sigs2',sigs2_tu_0,sig2_tu_true_B0)),
                              CI_table_helper(c('sigs2',sigs2_tu_1,sig2_tu_true_B1)),
                              CI_table_helper(c('sigs2',sigs2_td_0,sig2_td_true_B0)),
                              CI_table_helper(c('sigs2',sigs2_td_1,sig2_td_true_B1))))

length(as.numeric(res_df$wthin))


# Drop rownames
rownames(res_df) <- NULL
library(tidyverse)


res_df2 <- res_df %>%
  select(-all_of(c('truval','lwr','upr'))) %>%
  mutate(wthin = as.numeric(wthin),param_name=as.factor(param_name)) %>%
  group_by(param_name) %>%
  summarize(params_tot = n(),params_captured = sum(wthin), prop_captured=mean(wthin))

library(xtable)

# Convert dataframe to LaTeX table
latex_code <- xtable(res_df2)

# Print LaTeX code
print(latex_code, type = "latex")

##### Diagnostics raftery-lewis
library(coda)
library(ggplot2)

## effective samp and raftery-lewis

diags_table_fake <- create_diag_table(bets_0, bets_1, sigs2_time, thets_0, thets_1, 
                                 sigs2_0, sigs2_1, rans_eu_0, rans_eu_1, 
                                 rans_tu_0, rans_tu_1, rans_td_0, rans_td_1, 
                                 sigs2_eu_0, sigs2_eu_1, sigs2_tu_0, sigs2_tu_1,
                                 sigs2_td_0, sigs2_td_1,c('e','d','u'))
# rans_eu_1 doesn't seem to be able to get enough ESS

effectiveSize(bets_0[1,])
raftery.diag(bets_0)

################## Lots of different visuals
### Test to see if betas were captured
print(cbind(rowMeans(bets_0[,50:n]),B0,rowMeans(bets_1[,50:n]),B1))

plot(bets_0[1,38000:39000])

for(i in 1:nrow(bets_0)){
  # plot((bets_0[i,beg:n]),main=paste('Number of time points:',sum(fake_data$site_id==site_ids[i]),'     site id:',site_ids[i]),ylab='Posterior Beta_0 Draws')
  plot((bets_0[i,38000:39000]),type='l',main=paste0('Beta 0 for location', i),ylab='Value')
  
  plot(density(bets_0[i,]),main=paste('Number of time points:',sum(fake_data$site_id==site_ids[i]),'site id:',site_ids[i]))
  abline(v=B0[i])
}

for(i in 1:nrow(bets_0)){
  plot((bets_1[i,beg:n]),main=paste('Number of time points:',sum(fake_data$site_id==site_ids[i])))
  plot(density(bets_1[i,beg:n]),main=paste('Number of time points:',sum(fake_data$site_id==site_ids[i],"siteid",site_ids[i])))
  abline(v=B1[i])
}

### Tests to see if sig2_time converged and covered true val
print(cbind(rowMeans(sigs2_time[,beg:n]),sig2_time_true,site_ids))

plot((sigs2_time[23,]))
plot(density(sigs2_time[23,]))

for(i in 1:nrow(sigs2_time)){
  plot(density(sigs2_time[i,beg:n]),main=paste('Number of time points:',sum(fake_data$site_id==site_ids[i])))
  abline(v=sig2_time_true[i])
}

# Define the true values for parameters
true_values_B1 <- list(
  thets = c(0, -5, 0.003, -0.5),
  sigs2_tu = 9,
  sigs2_td = 4,
  sigs2_eu = 8,
  sigs2 = 1.5,
  rans_tu = 1410,
  rans_td = 1510,
  rans_eu = 0.1
)

true_values_B0 <- list(
  thets = c(5, -3, 0.01, 0.5),
  sigs2_tu = 16,
  sigs2_td = 36,
  sigs2_eu = 12,
  sigs2 = 14,
  rans_tu = 1600,
  rans_td = 1200,
  rans_eu = 0.15
)

# Plot for `thets_0` and `thets_1`
for (i in 1:nrow(thets_0)) { 
  plot(thets_0[i, 38000:39000], main = paste("Theta 0 - Parameter", i),
       xlab = "Index", ylab = "Value")
  abline(h = true_values_B0$thets[i], col = 'red')
  plot(density(thets_0[i, beg:n]), main = paste("Density of theta 0 - Parameter", i),
       xlab = paste("Value-True Val",true_values_B0$thets[i]), ylab = "Density")
  abline(v = true_values_B0$thets[i], col = 'red')
}

for (i in 1:nrow(thets_1)) {
  plot(thets_1[i, 38000:39000], type = 'l', main = paste("thets_1 - Parameter", i),
       xlab = "Index", ylab = "Value")
  # abline(h = true_values_B1$thets[i], col = 'red')
  plot(density(thets_1[i, ]), main = paste("Density of thets_1 - Parameter", i),
       xlab = paste("Value-True Val",true_values_B1$thets[i]), ylab = "Density")
  abline(v = true_values_B1$thets[i], col = 'red')
}



plot(1)
# Plot for all other parameters
parameters <- list(
  sigs2_0 = list(data = sigs2_0, true = true_values_B0$sigs2),
  sigs2_1 = list(data = sigs2_1, true = true_values_B1$sigs2),
  rans_eu_0 = list(data = rans_eu_0, true = true_values_B0$rans_eu),
  rans_eu_1 = list(data = rans_eu_1, true = true_values_B1$rans_eu),
  rans_tu_0 = list(data = rans_tu_0, true = true_values_B0$rans_tu),
  rans_tu_1 = list(data = rans_tu_1, true = true_values_B1$rans_tu),
  rans_td_0 = list(data = rans_td_0, true = true_values_B0$rans_td),
  rans_td_1 = list(data = rans_td_1, true = true_values_B1$rans_td),
  sigs2_eu_0 = list(data = sigs2_eu_0, true = true_values_B0$sigs2_eu),
  sigs2_eu_1 = list(data = sigs2_eu_1, true = true_values_B1$sigs2_eu),
  sigs2_tu_0 = list(data = sigs2_tu_0, true = true_values_B0$sigs2_tu),
  sigs2_tu_1 = list(data = sigs2_tu_1, true = true_values_B1$sigs2_tu),
  sigs2_td_0 = list(data = sigs2_td_0, true = true_values_B0$sigs2_td),
  sigs2_td_1 = list(data = sigs2_td_1, true = true_values_B1$sigs2_td)
)

for (param_name in names(parameters)) {
  param_info <- parameters[[param_name]]
  if(is.null(param_info$true)){
    param_info$true = NA
  }
  plot_parameter(param_info$data, param_name, beg, n,param_info$true)
}

### Get loo_cv (More important and simplefied later)
y_wit_id <- cbind(y,site_id = fake_data$site_id)
S <- n-1000
sbynmat <- log_likelihood_i(y_wit_id, S, T_mat, bets_0, bets_1, sigs2_time, site_ids)
dim(sbynmat)

sbynmat[1:10,1:10]
library(loo)
loo_result<-loo(sbynmat,cores=12)
bad_obs <- pareto_k_ids(loo_result, threshold = 0.7)
print(bad_obs)
hist(y)

####################################################
###############################################
# Now its time...
# to run the model on the real data!
###############################################
####################################################


set.seed(24)
full_df <- read.csv('fullDF.csv')
site_info <- read.csv('site_info.csv')
MegafirePointSamples <- read.csv('MegafirePointSamples.csv')
library(DataExplorer)

####################################################
# Data Cleaning
####################################################

static_data <- merge(
  full_df[, c('site_id', colnames(full_df)[!(colnames(full_df) %in% colnames(site_info))])], 
  site_info,
  by = 'site_id',
  all.x = TRUE # This ensures we only keep site_ids from site_info
)
static_data$site_id
data$site_id

## Get the data we have stations for
time_ser_our_data <- MegafirePointSamples[which(MegafirePointSamples$site_id %in% data$site_id),]
nrow(time_ser_our_data)

## What we'll do now is try to predict every variable with a time-series model taking into 
# account the day, and month, starting with the temperature to test things out
time_ser_our_data$datetime <- as.POSIXct(time_ser_our_data$datetime, format = '%m/%d/%Y %H:%M')

## 2016 is an error, get rid of it-it's the first row
time_ser_no_2016 <- time_ser_our_data[2:nrow(time_ser_our_data),]

# Getting rid of all rows with NA in DOC response,and time
time_ser_no_2016 <- time_ser_no_2016[!is.na(time_ser_no_2016$DOC_mg_L_SCAN) & !is.na(time_ser_no_2016$datetime),]
time_ser_no_2016$datetime


# time_ser_no_2016 <- time_ser_no_2016[!is.na(time_ser_no_2016$datetime),]
# percents2 <- sapply(time_ser_no_2016, function(x) sum(is.na(x))/nrow(time_ser_no_2016))
# percents2
time_dat <- time_ser_no_2016
time_dat$DOC_mg_L_SCAN

dim(time_dat)
y <- fix_0s(as.data.frame(time_dat$DOC_mg_L_SCAN))[,1] # turning vals <= 0 to be min of values>0
hist(log(y))

#Make time into the same format as have above, year.(day-1/365)
time_dat$datetime <- as.numeric(time_dat$datetime)
time_dat$datetime <- time_dat$datetime - as.numeric(as.POSIXct('01/01/2018 00:00', format = '%m/%d/%Y %H:%M'))
time_dat$datetime <- ((time_dat$datetime/(60*60*24))-1)/365

T_mat <- cbind(rep(1,nrow(time_dat)),time_dat$datetime)

## Forget static for now, lets use site_id because it has all site locations
####################################################
# Variable selection via correlation
####################################################

static_site_info<-site_info[site_info$site_id %in% data$site_id,]
summary(static_site_info)
vars_to_use <- c('river_length','burn_pct','lu_water_pct','lu_agri_pct','slp_median','Type2','ELEV','PRECIP')

static_site_info_2 <- static_site_info[,colnames(static_site_info)%in%vars_to_use]
static_site_info_2$river_length <- (static_site_info_2$river_length-mean(static_site_info_2$river_length))/sd(static_site_info_2$river_length)
colnames(static_site_info_2) <- c('Type2','river_length','burn_pct','lu_water_pct','lu_agri_pct','slp_median','ELEV','PRECIP')
em_static <- static_site_info_2
em_static[,-1] <- impute_missing_em(static_site_info_2[,-1])[[1]]
em_static[nrow(em_static),ncol(em_static)] <- 30.5
X_mat <- model.matrix(~Type2+river_length+burn_pct+lu_water_pct+lu_agri_pct+slp_median+ELEV+PRECIP,data=em_static)
dim(X_mat)

## Standardizing some cols
sndze <- function(col){
  (mean(col)-col)/sd(col)
}
X_mat[,c(4,8,9,10)] <- apply(X_mat[,c(4,8,9,10)],2,sndze)
X_mat[,c(5,6,7)] <- X_mat[,c(5,6,7)]/100

library(car)
plot_correlation(X_mat)
static_data_em <- static_data
static_data_em[,c('ELEV','PRECIP')] <- em_static[,c('ELEV','PRECIP')]
plot_correlation(static_data_em[,!colnames(static_data_em) %in% c('X','x','y')])
plot(X_mat[,'lu_agri_pct'],X_mat[,'river_length'])
plot((X_mat[,'lu_agri_pct']^(.1)),X_mat[,'river_length'])


temp <- as.data.frame(X_mat)

viftab <- vif(lm(rnorm(nrow(em_static))~Type2Natural_B+Type2Urban+river_length+burn_pct+lu_water_pct+lu_agri_pct+slp_median+ELEV+PRECIP,data=temp))
viftab <- as.matrix(viftab)
colnames(viftab) <- c('VIF')
class(viftab)
tail(static_data)
library(xtable)
print(xtable(as.matrix(viftab)), type = "latex",label='vif')

# percents3 <- sapply(static_site_info, function(x) sum(is.na(x))/nrow(static_site_info))
# static_site_info <- static_site_info[,percents3 < .0001]
# colnames(static_site_info)



# Get correlation plot of variables
library(DataExplorer)

site_ids <- static_site_info$site_id

###########################################################################
############### THE FULL SAMPLE (on REAL data) ############################
###########################################################################

m <- rep(0,ncol(X_mat))
V <- diag(100, nrow = ncol(X_mat))
V[1,1] <- 1000
V_inv <- solve(V)

sig2_time <- rep(4,nrow(X_mat))
betas_0 <- rep(0,nrow(X_mat))
betas_1 <- rep(0,nrow(X_mat))
p <- ncol(X_mat)
## For hierarchical
sig2_1 <- 20
sig2_eu_1 <- 50
sig2_tu_1 <- 20
sig2_td_1 <- 20
ran_tu_1 <- 1000
ran_td_1 <- 1000
ran_eu_1 <- 5
p <- ncol(X_mat)
thetas_1 <- rep(0,p)

sig2_0 <- 20
sig2_eu_0 <- 50
sig2_tu_0 <- 20
sig2_td_0 <- 20
ran_tu_0 <- 1000
ran_td_0 <- 1000
ran_eu_0 <- 5
thetas_0 <- rep(0,p)

n <- 40
# n <- 100

current_state <- list(sig2_time = sig2_time, sig2_1= sig2_1, sig2_eu_1 = sig2_eu_1, sig2_tu_1 = sig2_tu_1, sig2_td_1 = sig2_td_1, ran_tu_1 = ran_tu_1, ran_td_1 = ran_td_1, ran_eu_1 = ran_eu_1, betas_1 = betas_1, thetas_1 = thetas_1,
                      sig2_0 = sig2_0, sig2_eu_0 = sig2_eu_0, sig2_tu_0 = sig2_tu_0, sig2_td_0 = sig2_td_0, ran_tu_0 = ran_tu_0, ran_td_0 = ran_td_0, ran_eu_0 = ran_eu_0, betas_0 = betas_0, thetas_0 = thetas_0)


# Run each sampler in parallel

library(future)
library(loo)

y_wit_id <- cbind(log(y),site_id = time_dat$site_id)


plan(multisession)


mcmc1_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state, time_dat,c('u','d','e'),
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
  )
},seed = TRUE)

########### ########### ########### ########### ########### ########### 
########### Now for Taking out Euclidean
########### ########### ########### ########### ########### ########### 

mcmc2_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state, time_dat, c('u','d'),
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
  )
},seed = TRUE)

## TU

mcmc3_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state, time_dat, c('u'),
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
  )
},seed=TRUE)

## TD

mcmc4_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state, time_dat, c('d'),
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
  )
},seed=TRUE)

## E

mcmc5_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state, time_dat,c('e'),
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
  )
},seed=TRUE)

## Just sig2

mcmc6_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state, time_dat,c(''),
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
  )
},seed=TRUE)

res <- value(mcmc1_future)
resUD <- value(mcmc2_future)
resU <- value(mcmc3_future)
resD <- value(mcmc4_future)
resE <- value(mcmc5_future)
resS <- value(mcmc6_future)

plan(multisession)

loo_futureUDE <- (future({loo_funct(res)},seed=TRUE))
loo_futureUD <- (future({loo_funct(resUD)},seed=TRUE))
loo_futureU <- (future({loo_funct(resUD)},seed=TRUE))
loo_futureD <- (future({loo_funct(resU)},seed=TRUE))
loo_futureE <- (future({loo_funct(resD)},seed=TRUE))
loo_futureS <- (future({loo_funct(resE)},seed=TRUE))

loo_resultUDE <- value(loo_futureUDE)
loo_resultUD <- value(loo_futureUD)
loo_resultU <- value(loo_futureU)
loo_resultD <- value(loo_futureD)
loo_resultE <- value(loo_futureE)
loo_resultS <- value(loo_futureS)

## Make table of looic and SE

loo_table <- data.frame(
  model = c('UDE','UD','U','D','E','S'),
  looic = c(loo_resultUDE$estimates[3,1],loo_resultUD$estimates[3,1],loo_resultU$estimates[3,1],loo_resultD$estimates[3,1],loo_resultE$estimates[3,1],loo_resultS$estimates[3,1]),
  se = c(loo_resultUDE$estimates[3,2],loo_resultUD$estimates[3,2],loo_resultU$estimates[3,2],loo_resultD$estimates[3,2],loo_resultE$estimates[3,2],loo_resultS$estimates[3,2])
)

## Make into latex
print(xtable(loo_table), type = "latex",label='lootable')


# save(res, file = "real_data_chainno23.RData")
# save(resUD, file = "UDsampsno23.RData")
# save(resU, file = "Usampsno23.RData")
# save(resD, file = "Dsampsno23.RData")
# save(resE, file = "Esampsno23.RData")
# save(resS, file = "Ssampsno23.RData")

# load("real_data_chainno23.RData")
# load("UDsampsno23.RData")
# load("Usampsno23.RData")
# load("Dsampsno23.RData")
# load("Esampsno23.RData")
# load("Ssampsno23.RData")

bets_0 <- resD$bets_0
  # hist(rowMeans(bets_0)[-23])
bets_1 <- resD$bets_1
# hist(rowMeans(bets_1)[-23])
sigs2_time <- resD$sigs2_time

thets_0 <- resD$thets_0
thets_1 <- resD$thets_1
sigs2_0 <- resD$sigs2_0
sigs2_1 <- resD$sigs2_1
rans_eu_0 <- resD$rans_eu_0
rans_eu_1 <- resD$rans_eu_1
rans_tu_0 <- resD$rans_tu_0
rans_tu_1 <- resD$rans_tu_1
rans_td_0 <- resD$rans_td_0
rans_td_1 <- resD$rans_td_1

sigs2_tu_0 <- resD$sigs2_tu_0
sigs2_tu_1 <- resD$sigs2_tu_1
sigs2_td_0 <- resD$sigs2_td_0
sigs2_td_1 <- resD$sigs2_td_1
sigs2_eu_0 <- resD$sigs2_eu_0
sigs2_eu_1 <- resD$sigs2_eu_1



### Test to see if betas are similar to the OLS estimates via visuals
naive_coefs_0 <- c()
naive_coefs_1 <- c()
naive_sig2s <- c()
for (site in site_ids){
  y_lm <- log(y[time_dat$site_id==site])
  fit <- lm(y_lm~time_dat[time_dat$site_id==site,]$datetime)
  naive_coefs_0 <- c(naive_coefs_0,fit$coefficients[1])
  naive_coefs_1 <- c(naive_coefs_1,fit$coefficients[2])
  naive_sig2s <- c(naive_sig2s,sum((fit$residuals)^2)/(length(y_lm)-1))
}
#
for(i in 1:nrow(bets_0)){
  # plot((bets_0[i,beg:n]),main=paste('Number of time points:',sum(time_dat$site_id==site_ids[i])))
  plot(density(bets_0[i,]),main=paste('Number of time points:',sum(time_dat$site_id==site_ids[i]),'site_id',site_ids[i]))
  abline(v=naive_coefs_0[i],col='red')
  abline(v=mean(bets_0[i,]))

}

plot(1)
#
for(i in 1:nrow(bets_0)){
  # plot((bets_1[i,beg:n]),main=paste('Number of time points:',sum(time_dat$site_id==site_ids[i])))
  plot(density(bets_1[i,]),main=paste('Number of time points:',sum(time_dat$site_id==site_ids[i]),"site_id",site_ids[i]))
  abline(v=naive_coefs_1[i])
  abline(v=mean(bets_1[i,]))
  
}

plot(1)
### Tests to see if sig2_time converged
print(cbind(rowMeans(sigs2_time[,beg:n]),naive_sig2s,site_ids))
#
plot((sigs2_time[23,beg:n]))
plot(density(sigs2_time[23,beg:n]))

for(i in 1:nrow(sigs2_time)){
  plot(density(sigs2_time[i,beg:n]),main=paste('Number of time points:',sum(time_dat$site_id==site_ids[i])))
  abline(v=naive_sig2s[i])
}
#
#
#
#
#
#
plot(1)
# # Plot for all other parameters
parameters <- list(
  sigs2_0 = list(data = sigs2_0),
  sigs2_1 = list(data = sigs2_1),
  rans_eu_0 = list(data = rans_eu_0),
  rans_eu_1 = list(data = rans_eu_1),
  rans_tu_0 = list(data = rans_tu_0),
  rans_tu_1 = list(data = rans_tu_1),
  rans_td_0 = list(data = rans_td_0),
  rans_td_1 = list(data = rans_td_1),
  sigs2_eu_0 = list(data = sigs2_eu_0),
  sigs2_eu_1 = list(data = sigs2_eu_1),
  sigs2_tu_0 = list(data = sigs2_tu_0),
  sigs2_tu_1 = list(data = sigs2_tu_1),
  sigs2_td_0 = list(data = sigs2_td_0),
  sigs2_td_1 = list(data = sigs2_td_1)
)
# names(parameters)
for (param_name in names(parameters)) {
  param_info <- parameters[[param_name]]
  print(param_name)
  plot_parameter(param_info$data, param_name, beg, n)
}

##### Eff samp size, raftery-lewis, loocv
library(coda)
# Load the saved result correctly

bets_0 <- resD$bets_0[,1001:ncol(resD$bets_0)]
bets_1 <- resD$bets_1[,1001:ncol(resD$bets_0)]
sigs2_time <- resD$sigs2_time[,1001:ncol(resD$bets_0)]

thets_0 <- resD$thets_0[,1001:ncol(resD$bets_0)]
thets_1 <- resD$thets_1[,1001:ncol(resD$bets_0)]
sigs2_0 <- resD$sigs2_0[1001:ncol(resD$bets_0)]
sigs2_1 <- resD$sigs2_1[1001:ncol(resD$bets_0)]

rans_eu_0 <- resD$rans_eu_0[1001:ncol(resD$bets_0)]
rans_eu_1 <- resD$rans_eu_1[1001:ncol(resD$bets_0)]
rans_tu_0 <- resD$rans_tu_0[1001:ncol(resD$bets_0)]
rans_tu_1 <- resD$rans_tu_1[1001:ncol(resD$bets_0)]
rans_td_0 <- resD$rans_td_0[1001:ncol(resD$bets_0)]
rans_td_1 <- resD$rans_td_1[1001:ncol(resD$bets_0)]

sigs2_tu_0 <- resD$sigs2_tu_0[1001:ncol(resD$bets_0)]
sigs2_tu_1 <- resD$sigs2_tu_1[1001:ncol(resD$bets_0)]
sigs2_td_0 <- resD$sigs2_td_0[1001:ncol(resD$bets_0)]
sigs2_td_1 <- resD$sigs2_td_1[1001:ncol(resD$bets_0)]
sigs2_eu_0 <- resD$sigs2_eu_0[1001:ncol(resD$bets_0)]
sigs2_eu_1 <- resD$sigs2_eu_1[1001:ncol(resD$bets_0)]

diags_table <- create_diag_table(bets_0, bets_1, sigs2_time, thets_0, thets_1,
                                  sigs2_0, sigs2_1, rans_eu_0, rans_eu_1,
                                  rans_tu_0, rans_tu_1, rans_td_0, rans_td_1,
                                  sigs2_eu_0, sigs2_eu_1, sigs2_tu_0, sigs2_tu_1,
                                  sigs2_td_0, sigs2_td_1,cov_to_calc = c('d'))

## Plot the Thetaaaas!! 
# Lets take a look ourselves real quick
# Plot for `thets_0` and `thets_1`
df_to_plot_0 <- data.frame(meds = apply(thets_0,1,median),
                           upr = apply(thets_0,1,quantile,probs = .975),
                           lwr = apply(thets_0,1,quantile,probs = .025),
                           varname = colnames(X_mat))

df_to_plot_1 <- data.frame(meds = apply(thets_1,1,median),
                           upr = apply(thets_1,1,quantile,probs = .975),
                           lwr = apply(thets_1,1,quantile,probs = .025),
                           varname = colnames(X_mat))

for (i in 1:nrow(thets_0)) {
  # plot(thets_0[i,], type = 'l', main = paste("thets_0 - Parameter", i),
  #      xlab = "Index", ylab = "Value")
  plot(density(thets_0[i,]), main = paste("Density of thets_0 -", colnames(X_mat)[i]),
       xlab = paste("Value"), ylab = "Density")
  abline(v=df_to_plot_0$meds[i])
  abline(v=df_to_plot_0$upr[i],col='red')
  abline(v=df_to_plot_0$lwr[i],col='red')
  # plot((thets_0[i,38000:38999]), type = 'l',main = paste("Density of thets_0 -", colnames(X_mat)[i]),
  #      xlab = paste("Value"), ylab = "Density")
}

for (i in 1:nrow(thets_1)) {
  # plot(thets_1[i,], type = 'l', main = paste("thets_1 - Parameter", i),
  # xlab = "Index", ylab = "Value")
  plot(density(thets_1[i,]), main = paste("Density of thets_1 -", colnames(X_mat)[i]),
       xlab = paste("Value"), ylab = "Density")
  abline(v=df_to_plot_1$meds[i])
  abline(v=df_to_plot_1$upr[i],col='red')
  abline(v=df_to_plot_1$lwr[i],col='red')
  plot(density(thets_1[i,]), main = paste("Density of thets_1 -", colnames(X_mat)[i]),
       xlab = paste("Value"), ylab = "Density")
}

##############################
# The intercepts:
## Some help from chatGPT for this visual
ggplot(data = df_to_plot_0, aes(x = meds, y = varname)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) + 
  geom_point(size = 3) +
  labs(
    title = 'Static Variable Coefficients for Intercept',
    x = 'Value',
    y = 'Variable'
  ) +
  xlim(-20, 20) + 
  theme_minimal()

ggplot(data = df_to_plot_1, aes(x = meds, y = varname)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) + 
  geom_point(size = 3) +
  labs(
    title = 'Static Variable Coefficients for Slope',
    x = 'Value',
    y = 'Variable'
  ) +
  xlim(-18, 18) + 
  
  theme_minimal()

####################################################
# Intercept vis
####################################################

vals <- rowMeans(bets_0)
response_values <- c(vals[1:22],NA,vals[23:length(vals)])

mf04p@obspoints@SSNPoints[[1]]@point.data$response <- response_values

# Define number of bins (e.g., 7 for clarity)
num_bins <- 6

# Create bins based on response values
bin_breaks <- quantile(response_values, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
binned_responses <- cut(response_values, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)

# Generate color palette
color_palette <- colorRampPalette(c("#F5DEB3", "#D2691E", "#8B0000"))(num_bins+1)
# color_palette[6:7] <- '#468C8E'
# Assign colors based on bins
point_colors <- color_palette[binned_responses]

# Extract coordinates of observation points
point_coords <- mf04p@obspoints@SSNPoints[[1]]@point.data[, c("NEAR_X", "NEAR_Y")]

# Plot the SSN object WITHOUT points
plot(mf04p,
     lwdLineCol = "afvArea",
     lwdLineEx = 5,
     lineCol = "#468C8E",  # Adjusted blue-green for better aesthetics
     addWithLegend = FALSE,  # Manually adding the legend
     xlab = "x-coordinate (m)", 
     ylab = "y-coordinate (m)", 
     xlim = c(400000, 505000),  
     ylim = range(point_coords$NEAR_Y) + c(-16500, 17500),
     main = 'Visual of log(DOC) Intercept by Location'
)

# Overlay points with binned colors
points(point_coords$NEAR_X, point_coords$NEAR_Y, xlim = c(400000,505000),  # Ensure the full range of x is shown
       ylim = y_limits+c(-16500,17500),
       col = point_colors, pch = 19, cex = 1.8)

# Add a legend
legend("topright", 
       legend = round(bin_breaks, 2), 
       col = color_palette, 
       pch = 19, 
       title = "log(DOC)")

#################################################### 
# Vis for SLOPE
####################################################

# Extract response variable
vals <- rowMeans(bets_1)
response_values <- c(vals[1:22],NA,vals[23:length(vals)])

mf04p@obspoints@SSNPoints[[1]]@point.data$response <- response_values

# Define number of bins (e.g., 7 for clarity)
num_bins <- 6

# Create bins based on response values
bin_breaks <- quantile(response_values, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
binned_responses <- cut(response_values, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)

# Generate color palette
color_palette <- colorRampPalette(c("#F5DEB3", "#D2691E", "#8B0000"))(num_bins+1)
# color_palette[6:7] <- '#468C8E'
# Assign colors based on bins
point_colors <- color_palette[binned_responses]

# Extract coordinates of observation points
point_coords <- mf04p@obspoints@SSNPoints[[1]]@point.data[, c("NEAR_X", "NEAR_Y")]

# Plot the SSN object WITHOUT points
plot(mf04p,
     lwdLineCol = "afvArea",
     lwdLineEx = 5,
     lineCol = "#468C8E",  # Adjusted blue-green for better aesthetics
     addWithLegend = FALSE,  # Manually adding the legend
     xlab = "x-coordinate (m)", 
     ylab = "y-coordinate (m)", 
     xlim = c(400000, 505000),  
     ylim = range(point_coords$NEAR_Y) + c(-16500, 17500),
     main = 'Visual of log(DOC) Slope by Location'
)

# Overlay points with binned colors
points(point_coords$NEAR_X, point_coords$NEAR_Y, xlim = c(400000,505000),  # Ensure the full range of x is shown
       ylim = y_limits+c(-16500,17500),
       col = point_colors, pch = 19, cex = 1.8)

# Add a legend
legend("topright", 
       legend = round(bin_breaks, 2), 
       col = color_palette, 
       pch = 19, 
       title = "log(DOC)")

##################
## Graphs (make sure the right bets_0 and bets_1 is in global env)
##################
lin_vis(5)
lin_vis(1)
lin_vis(34)
for(i in 1:51){
  lin_vis(i)
}
length(time_dat$datetime[time_dat$site_id==site_ids[3]])
table(time_dat$site_id)

###################
### Visualization of Natural B land
##################

# Plot the SSN object WITHOUT points
static_data2 <- static_data[match(mf04p@obspoints@SSNPoints[[1]]@point.data$site_id,static_data$site_id),]


## Get the data we have stations for
static_site_info3 <-site_info[site_info$site_id %in% mf04p@obspoints@SSNPoints[[1]]@point.data$site_id,]
static_site_info3 <- static_site_info3[match(mf04p@obspoints@SSNPoints[[1]]@point.data$site_id,static_site_info3$site_id),]
static_site_info3$site_id
# Overlay points with unique colors


plot(mf04p,
     lwdLineCol = "afvArea",
     lwdLineEx = 5,
     lineCol = "#4682B4",
     addWithLegend = FALSE,  # We will manually add the points
     xlab = "x-coordinate (m)", 
     ylab = "y-coordinate (m)", 
     xlim = c(400000,505000),  # Ensure the full range of x is shown
     ylim = y_limits+c(-14500,16500),
     main = "Site Colored by Differing Lines From OLS"
)

points(point_coords$NEAR_X, point_coords$NEAR_Y, xlim=c(400000,505000),ylim=y_limits+c(-14500,16500),
       col = ifelse(static_site_info3[, 'Type2'] == 'Natural_B', 'red',
                    ifelse(static_site_info3[, 'Type2'] == 'Natural', 'black', 'grey')), pch = 19, cex = 1.5)


diffs <- c(1,10,13,21,27+1,45+1,46+1,47+1)
points(point_coords$NEAR_X, point_coords$NEAR_Y, xlim=c(400000,505000),ylim=y_limits+c(-14500,16500),
       col = ifelse(1:52 %in% diffs,'red','black'), pch = 19, cex = 1.5)

legend("topright", 
       legend = unique(static_site_info3[,'Type2']), 
       col = c('black','red','grey'),
       pch = 19, title = "Land Type")

## Visualizing which points differed from OLS
for(i in c(1,10,13,21,27,45,46,47)){
  lin_vis(i)
}
static_data

## Last minute gelman-ruben diag

# lets do four more chains of tail down

library(future)
library(loo)

y_wit_id <- cbind(log(y),site_id = time_dat$site_id)


sig2_time1 <- rep(4,nrow(X_mat))
betas_01 <- rep(0,nrow(X_mat))
betas_11 <- rep(0,nrow(X_mat))
p1 <- ncol(X_mat)
## For hierarchical
sig2_11 <- 20
sig2_eu_11 <- 50
sig2_tu_11 <- 20
sig2_td_11 <- 20
ran_tu_11 <- 1000
ran_td_11 <- 1000
ran_eu_11 <- 5
p1 <- ncol(X_mat)
thetas_11 <- rep(0,p)

sig2_01 <- 20
sig2_eu_01 <- 50
sig2_tu_01 <- 20
sig2_td_01 <- 20
ran_tu_01 <- 1000
ran_td_01 <- 1000
ran_eu_01 <- 5
thetas_01 <- rep(0,p)

n <- 40000
# n <- 100

current_state1 <- list(sig2_time = sig2_time1, sig2_1= sig2_11, sig2_eu_1 = sig2_eu_11, sig2_tu_1 = sig2_tu_11, sig2_td_1 = sig2_td_11, ran_tu_1 = ran_tu_11, ran_td_1 = ran_td_11, ran_eu_1 = ran_eu_11, betas_1 = betas_11, thetas_1 = thetas_11,
                      sig2_0 = sig2_01, sig2_eu_0 = sig2_eu_01, sig2_tu_0 = sig2_tu_01, sig2_td_0 = sig2_td_01, ran_tu_0 = ran_tu_01, ran_td_0 = ran_td_01, ran_eu_0 = ran_eu_01, betas_0 = betas_01, thetas_0 = thetas_01)




sig2_time2 <- rep(2,nrow(X_mat))
betas_02 <- rep(1,nrow(X_mat))
betas_12 <- rep(1,nrow(X_mat))
p1 <- ncol(X_mat)
## For hierarchical
sig2_12 <- 30
sig2_eu_12 <- 60
sig2_tu_12 <- 30
sig2_td_12 <- 30
ran_tu_12 <- 1500
ran_td_12 <- 1500
ran_eu_12 <- 10
p <- ncol(X_mat)
thetas_12 <- rep(1,p)

sig2_02 <- 30
sig2_eu_02 <- 60
sig2_tu_02 <- 30
sig2_td_02 <- 30
ran_tu_02 <- 1500
ran_td_02 <- 1500
ran_eu_02 <- 10
thetas_02 <- rep(0,p)

n <- 40000

current_state2 <- list(sig2_time = sig2_time2, sig2_1= sig2_12, sig2_eu_1 = sig2_eu_12, sig2_tu_1 = sig2_tu_12, sig2_td_1 = sig2_td_12, ran_tu_1 = ran_tu_12, ran_td_1 = ran_td_12, ran_eu_1 = ran_eu_12, betas_1 = betas_12, thetas_1 = thetas_12,
                      sig2_0 = sig2_02, sig2_eu_0 = sig2_eu_02, sig2_tu_0 = sig2_tu_02, sig2_td_0 = sig2_td_02, ran_tu_0 = ran_tu_02, ran_td_0 = ran_td_02, ran_eu_0 = ran_eu_02, betas_0 = betas_02, thetas_0 = thetas_02)


sig2_time3 <- rep(4,nrow(X_mat))
betas_03 <- rep(0,nrow(X_mat))
betas_13 <- rep(0,nrow(X_mat))
## For hierarchical
sig2_13 <- 10
sig2_eu_13 <- 40
sig2_tu_13 <- 10
sig2_td_13 <- 10
ran_tu_13 <- 700
ran_td_13 <- 700
ran_eu_13 <- 2
p <- ncol(X_mat)
thetas_13 <- rep(0,p)

sig2_03 <- 10
sig2_eu_03 <- 40
sig2_tu_03 <- 10
sig2_td_03 <- 10
ran_tu_03 <- 700
ran_td_03 <- 700
ran_eu_03 <- 3
thetas_03 <- rep(0,p)

n <- 40000

current_state3 <- list(sig2_time = sig2_time3, sig2_1= sig2_13, sig2_eu_1 = sig2_eu_13, sig2_tu_1 = sig2_tu_13, sig2_td_1 = sig2_td_13, ran_tu_1 = ran_tu_13, ran_td_1 = ran_td_13, ran_eu_1 = ran_eu_13, betas_1 = betas_13, thetas_1 = thetas_13,
                      sig2_0 = sig2_03, sig2_eu_0 = sig2_eu_03, sig2_tu_0 = sig2_tu_03, sig2_td_0 = sig2_td_03, ran_tu_0 = ran_tu_03, ran_td_0 = ran_td_03, ran_eu_0 = ran_eu_03, betas_0 = betas_03, thetas_0 = thetas_03) 

plan(multisession)


mcmc7_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state1, time_dat,c('d'),
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
  )
},seed=FALSE)

mcmc8_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_stat2, time_dat,c('d'),
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
  )
},seed=FALSE)

mcmc9_future <- future({
  Sampler(n, log(y), T_mat,site_ids, X_mat, current_state3, time_dat,c('d'),
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
  )
},seed=FALSE)



res7 <- value(mcmc7_future)
res8 <- value(mcmc8_future)
res9 <- value(mcmc9_future)




gelmanrubin_table(res7,res8,res9)
