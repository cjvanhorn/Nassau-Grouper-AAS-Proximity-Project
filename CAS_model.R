# December 2, 2022
# Author: Cameron Van Horn
# Email: cjvanhor@ucsd.edu
#############################################

# Below is code to fit a model for assessing spatial and temporal drivers 
# of sound production in Nassau Grouper

# The 'Rethinking' package is the machine of the model, and requires extra
# software to operate. You must:
   # install 'rstan' from mc-stan.org
   # install a C++ compiler (the 'tool chain') also at mc-stan.org
   # within R, install rethinking using the following code:
      # install.packages(c('coda', 'mvtnorm', 'devtools'))
      # library(devtools)
      # devtools::install_github('rmcelreath/rethinking', ref = 'Experimental')
   # If any unexpected issues arise, reference github.com/rmcelreath/rethinking
library(rethinking)
library(RColorBrewer)

# set working directory
setwd("~/Desktop/NG CAS_Model")

# read in the data (found in the github)
data <- read.csv('model_data.csv')
head(data) # 158 hours (data rows) for 5 hydrophones

# for the model to function properly, we must perform all variable 
# transformations outside of the model
# Also, we will need a list of data containing only the required variables
# create list for model from data
data_list <- list(
   tothrs = length(unique(data$hour.of.period)), # cumulative hours monitored
   hr = as.integer(data$hour.of.period), # hour of study period (1 - 158)
   hydro = as.integer(as.factor(data$Hydro)), # 1 = DSG1, 2 = ST2, etc.
   cas = data$cas.count, # counts of CAS 
   m = data$minutes.observed / 60, # minutes annotated per hour conversion
   timeofday = as.integer(data$Hour + 1), # add 1 to allow logs of hour 0
   DAFS = as.integer(as.factor(data$DAFS)), # Days After First Spawn
   fishProx = as.integer(as.factor(data$fish.prox)) # 1 = not observed, 
                                                    # 2 = nearby (20-100m), 
                                                    # 3 = present (0-20m)
)

#######################
##### Build Model #####
#######################

set.seed(472787) # set seed for duplication
n_iter = 1000 # number of iterations per chain
n_chains = 4 # number of chains in MCMC
n_cores = 4 # number of cores used by cpu, adjust accordingly

# use ulam() function to compile lists of formulas for Stan (think of it as
# as a translator)

mod <- { ulam(
   alist(
      # CAS counts are Poisson-distributed with a rate of lambda
      cas ~ dpois(lambda),
      # use a log link function to restrict CAS production rates as positive
      # estimate CAS production rates as a function of
         # time of day (timeofday)
         # days after first spawn (DAFS)
         # fish proximity to each hydrophone (fishProx)
         # random effect of hour of study within each hydrophone (hr,hydro)*
      # and multiply all above terms by minutes annotated per hour (m)
      # * a covariance matrix where each hydrophone has a separate variance,
      # * and covariance of hydrophones is assumed to be fixed across hours
      # * This accounts for correlations in CAS detections across hydrophones
      # * that arise from limited detection ranges and aggregating behaviors
      log(lambda) <- (d[hr,hydro] + bt[timeofday] + bd[DAFS] + 
                         fp[fishProx]) * m, 
      
      # to facilitate mixing among chains, we factor all parameters out of 
      # adaptive priors and insert them into the model
      # First, compose non-centered adaptive prior for covariance matrix d
      # through a Cholesky decomposition
      # Define matrix of varying effects of total hours monitored (tothrs)
      # to be transformed
      # Merge stdev vector (sigma_d) with Cholesky correlation factor (L_Rho_d)
      # multiplied by a z-score matrix (z) defined by dnorm(0,1)
      # This yields a covariance matrix with random effects on correct scale
      # for the model
      transpars> matrix[tothrs, 5]:d <-
         compose_noncentered(rep_vector(sigma_d, 5), L_Rho_d, z), 
      matrix[5, tothrs]:z ~ dnorm(0, 1),
      
      # assign prior for sigma_d as dexp(1)
      # define prior for L_Rho_d with a Cholesky correlation factor of 8
      sigma_d ~ dexp(1),
      cholesky_factor_corr[5]:L_Rho_d ~ lkj_corr_cholesky(8),
      
      # compute correlation matrices from L_Rho_d to interpret correlations
      # at the end of each transition
      gq> matrix[5, 5]:Rho_d <<- Chol_to_Corr(L_Rho_d),
      
      # assign vague priors for parameters timeofday, DAFS, and fishProx
      # normally distributed centered on 0 with stdev of 1
      bt[timeofday] ~ dnorm(0, 1),
      bd[DAFS] ~ dnorm(0, 1),
      fp[fishProx] ~ dnorm(0, 1)
      
      # increase max_treedepth to 14 to allow for adequate simulation steps
      # for model fit
      # adapt_delta = 0.95 as default
   ), data = data_list, chains = n_chains, cores = n_cores, log_lik = T, 
   iter = n_iter, control = list(max_treedepth = 14, adapt_delta = 0.95)) }

# inspect posterior summary
mod_sum <- precis(mod, depth = 3, pars = c("Rho_d", "sigma_d", "bt", "bd", "fp"), 
                  prob = 0.95)
mod_sum

# high n_eff and Rhat4 ~ 1 indicate convergence and mixing

# save posterior summary
write.csv(mod_sum, 'CAS model_summary.csv')

###############################
##### VISUALIZE POSTERIOR #####
###############################

# extract posterior samples for parameters timeofday, DAFS, and fishProx
post <- extract.samples(mod)

post_t <- exp(post$bt)
post_d <- exp(post$bd)
post_fp <- exp(post$fp)

# adjust column names for plotting
# readjust timeofday back to hours 0-23 (from 1-24)
colnames(post_t) <- c(0:23)
# readjust DAFS back to numeric centered on day of first spawn (from 1-7)
# now day 1 is day -4, and the day of first spawn is day 0, and last day is
# day 2
colnames(post_d) <- c(-4:2)
# assign descriptors to categories of fish proximities
colnames(post_fp) <- c('Not Observed', 'Nearby (20-100 m)', 'Present (0-20 m)')

# establish plotting colors
colors <- brewer.pal(n = 11, name = 'BrBG') # 1:5 brown, 6 white, 7:11 bluegreen

# box plot of predictions for time of day
# polygons represent crepusclar hours
boxplot(post_t, outline = F, xlab = 'Hour', 
        ylab = expression('Predicted CAS h'^-1), col = colors[9])
polygon(c(-0.42, -0.42, 7, 7), c(-1.65, 48, 48, -1.65), 
        border = 'dark grey', col = 'dark grey')
polygon(c(20, 20, 25.43, 25.43), c(-1.65, 48, 48, -1.65), 
        border = 'dark grey', col = 'dark grey')
polygon(c(7, 7, 8, 8), c(-1.65, 48, 48, -1.65), 
        border = 'light grey', col = 'light grey')
polygon(c(19, 19, 20, 20), c(-1.65, 48, 48, -1.65), 
        border = 'light grey', col = 'light grey')
par(new = T)
boxplot(post_t, outline = F, xlab = 'Hour', 
        ylab = expression('Predicted CAS h'^-1), col = colors[9])

# box plot of predictions for DAFS,
# with full moon signaling
boxplot(post_d, outline = F, xlab = 'Days Before or After First Spawn (DAFS)', 
        ylab = expression('Predicted CAS h'^-1), col = colors[7]) 
points(2, 100, pch = 19, cex = 5, col = 'Yellow')
points(2, 100, cex = 5)


# box plot of predictions for fish proximity
boxplot(post_fp, outline = F, xlab = '', 
        ylab = expression('Predicted CAS h'^-1), col = colors[4])
