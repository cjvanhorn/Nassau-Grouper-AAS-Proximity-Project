# Plotting and Modeling Nassau Grouper acoustic data from Little Cayman's
# West End during the 2020 winter spawning season

# Contact: Cameron Van Horn
#          (949) 233-8671
#          cameronvanhorn@gmail.com

##############################
##### TERMS AND ACRONYMS #####
##############################

# NG == Nassau Grouper
# AAS == Aggregation-associated sounds
# DAFS == Days After First Spawn
# ST == SoundTrap hydrophone
# LS == Loggerhead Instruments LS1 hydrophone

##########################################
##### LIBRARIES AND LOADING THE DATA #####
##########################################

# The 'Rethinking' package is the machine of the model, and requires extra
# software to operate. You must:
   # Install 'rstan' from mc-stan.org
   # Install a C++ compiler (the 'tool chain') also at mc-stan.org
   # Within R, install rethinking using the following code:
      # install.packages(c('coda', 'mvtnorm', 'devtools'))
      # library(devtools)
      # devtools::install_github('rmcelreath/rethinking', ref = 'Experimental')
# If any unexpected issues arise, reference github.com/rmcelreath/rethinking
library(rethinking)
if(!require('RColorBrewer')) install.packages('RColorBrewer') 
if(!require('ggplot2')) install.packages('ggplot2') 
if(!require('corrplot')) install.packages('corrplot') 
if(!require('tidyverse')) install.packages('tidyverse') 
if(!require('PBSmapping')) install.packages('PBSmapping') 
if(!require('raster')) install.packages('raster')
if(!require('sf')) install.packages('sf')
if(!require('units')) install.packages('units')
if(!require('smoothr')) install.packages('smoothr')
if(!require('marmap')) install.packages('marmap')
if(!require('metR')) install.packages('metR')
if(!require('ggspatial')) install.packages('ggspatial')
if(!require('hrbrthemes')) install.packages('hrbrthemes')
if(!require('dplyr')) install.packages('dplyr')
if(!require('tidyr')) install.packages('tidyr')
if(!require('patchwork')) install.packages('patchwork')
if(!require('gridGraphics')) install.packages('gridGraphics')
if(!require('knitr')) install.packages('knitr')
if(!require('stats')) install.packages('stats')
if(!require('hms')) install.packages('hms')

# Set working directory
# TODO: make a generic working directory call
setwd("/Users/cameronvanhorn/Documents/MASTERS/Important/PAPER RELATED THINGS")

# clear memory
rm(list = ls())

# Read in the data (found in the github)
load('data_for_figures_and_model.RData')

# functions for plotting maps
scale_x_longitude <- function(xmin = -180, xmax = 180, step = 1, ...) {
   xbreaks <- seq(xmin, xmax, step)
   xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text = paste0(sqrt(x^2), ' ^o', '*W')), ifelse(x > 0, parse(text = paste0(x, ' ^o', '*E')), x))))
   return(scale_x_continuous('', breaks = xbreaks, labels = xlabels, expand = c(0.08, 0), ...))
}
scale_y_latitude <- function(ymin = -90, ymax = 90, step = 0.5, ...) {
   ybreaks <- seq(ymin, ymax, step)
   ylabels <- unlist(lapply(ybreaks, function(x) ifelse (x < 0, parse(text = paste0(x, ' ^o', '*S')), ifelse(x > 0, parse(text = paste0(x, ' ^o', '*N')), x))))
   return(scale_y_continuous('', breaks = ybreaks, labels = ylabels, expand = c(0.02, 0), position = 'right'))
}

#####################
##### THE MODEL #####
#####################

# For the model to function properly, we must perform all variable 
# transformations outside of the model
# Also, we will need a list of data containing only the required variables
# Create list for model from data
data_list <- list(
   
   # cumulative hours observed for all hydrophones
   N = nrow(model_formatted_data), 
   
   # total hours monitored during the study period
   tothrs = length(unique(model_formatted_data$HOUR_OF_PERIOD)), 
   
   # hour of study period (1 - 158)
   hr = as.integer(model_formatted_data$HOUR_OF_PERIOD), 
   
   # 1 = LARS1, 2 = ST2, etc.
   hydro = as.integer(as.factor(model_formatted_data$HYDRO_STATION)), 
   # Note there is no 3 because LARS3 malfunctioned and did not record any data
   
   # counts of AAS 
   aas = model_formatted_data$AAS_COUNT, 
   
   # minutes observed per hour conversion
   m = model_formatted_data$MINUTES_OBSERVED / 60, 
   
   # add 1 to allow logs of hour 0
   timeofday = as.integer(model_formatted_data$HOUR_OF_DAY + 1), 
   # (12:00AM == 1 ... 11:00PM == 24)
   
   # Days After First Spawn (4, 3, 2, 1 before spawn; 5, 6, 7 after spawn)
   DayOfSeason = as.integer(as.factor(model_formatted_data$DAFS)), 
   
   # 1 = NG not observed, 2 = NG nearby (20-100m), 3 = NG present (0-20m)
   fishProx = as.integer(as.factor(model_formatted_data$FISH_PROX)) 

                  )

# set seed for duplication
set.seed(234) 

# number of iterations per chain
n_iter = 1000 

# number of chains in MCMC
n_chains = 3 

# before the model can be run, install cmdstan locally
# install_cmdstan()

# Use ulam() function to compile lists of formulas for Stan (think of it as
# as a translator)
mod <- { 
   ulam(
      alist(
         # AAS counts are Poisson-distributed with a rate of lambda.
         aas ~ poisson(lambda),
         # use a log link function to restrict AAS production rates as positive
         # estimate AAS production rates as a function of
            # time of day (timeofday)
            # days after first spawn (DAFS)
            # fish proximity to each hydrophone (fishProx)
            # random effect of hour of study within each hydrophone (hr,hydro)*
         # log(m) serves as an exposure term that accounts for the fact that
         # the number of minutes observed per hour varies.
         # * a covariance matrix where each hydrophone has a separate variance,
         # * and covariance of hydrophones is assumed to be fixed across hours.
         # * This accounts for correlations in AAS detections across hydrophones
         # * that arise from limited detection ranges and aggregating behaviors.
         log(lambda) <- d[hr,hydro] + bt[timeofday] + bd[DayOfSeason] + 
            fp[fishProx] + log(m), 
         
         # to facilitate mixing among chains, we factor all parameters out of 
         # adaptive priors and insert them into the model
         # First, compose non-centered adaptive prior for covariance matrix d
         # through a Cholesky decomposition
         # Define matrix of varying effects of total hours monitored (tothrs)
         # to be transformed
         # Merge stdev vector (sigma_d) with Cholesky correlation factor 
         # (L_Rho_d) multiplied by a z-score matrix (z) defined by dnorm(0,1)
         # This yields a covariance matrix with random effects on correct scale
         # for the model
         transpars> matrix[tothrs, 5]:d <-
            compose_noncentered(rep_vector(sigma_d, 5), L_Rho_d, z), 
         matrix[5, tothrs]:z ~ normal(0, 1),
      
         # assign prior for sigma_d as dexp(1)
         # define prior for L_Rho_d with a Cholesky correlation factor of 1
         sigma_d ~ exponential(1),
         cholesky_factor_corr[5]:L_Rho_d ~ lkj_corr_cholesky(1),
         
         # compute correlation matrices from L_Rho_d to interpret correlations
         # at the end of each transition
         gq> matrix[5, 5]:Rho_d <<- Chol_to_Corr(L_Rho_d),
         
         # assign vague priors for parameters timeofday, DAFS, and fishProx
         # normally distributed centered on 0 with stdev of 1
         bt[timeofday] ~ normal(0, 1),
         bd[DayOfSeason] ~ normal(0, 1),
         fp[fishProx] ~ normal(0, 1)
         
         # increase max_treedepth to 15 to allow for adequate simulation steps
         # for model fit
         # adapt_delta = 0.95 as default
            ), data = data_list, chains = n_chains, log_lik = T, 
               iter = n_iter, control = list(max_treedepth = 15, 
                                             adapt_delta = 0.95)
      ) 
   }

# inspect posterior summary
mod_sum <- precis(mod, depth = 3, pars = c("Rho_d", "sigma_d", 
                                           "bt", "bd", "fp"), 
                  prob = 0.95)

mod_sum
# high n_eff and Rhat4 ~ 1 indicate convergence and mixing

# save posterior summary
write.csv(mod_sum, 'AAS_model_summary.csv')

##########################################
##### FIGURE 1: Map of Little Cayman #####
##########################################
# Credit and thanks are given to Brian Stock for his help in writing the code
# responsible for this figure
# Map Section A (Caribbean) ----------------------------------------------------
lons <- c(-89, -74)
lats <- c(15, 27)

carib_bathy[carib_bathy > -1] <- -1

worldLLfull$X <- worldLLfull$X - 360

labs <- data.frame(x = c(-81, -80.5, -88.15),
                   y = c(26.5, 22.4, 20.8),
                   region = c('Florida', 'Cuba', 'Mexico'))

xbreaks <- seq(lons[1], lons[2], 5)
xlabels <- paste0(xbreaks, 'º', 'W')

ybreaks <- seq(lats[1], lats[2], 4)
ylabels <- paste0(ybreaks, 'º', 'N')

# Make the plot
panelA <- {
   ggplot(carib_bathy, aes(x = x, y = y)) + 
   coord_quickmap() +
   geom_polygon(data = worldLLfull, 
                aes(x = X, 
                    y = Y, 
                    group = PID), 
                fill = "grey", 
                color = "black") +
   scale_size_manual(values = c(4,2.5), 
                     guide = 'none') +
   scale_shape_manual(values = c(17,19), 
                      guide = 'none') +
   geom_text(data = labs, 
             aes(label = region), 
             size = 10) +
   xlab('') +
   ylab('') +
   scale_x_continuous(limits = lons, 
                      position = 'top', 
                      breaks = xbreaks, 
                      labels = xlabels) +
   scale_y_continuous(limits = lats, 
                      position = "left", 
                      breaks = ybreaks, 
                      labels = ylabels) +  
   theme_bw() +
   theme(panel.border = element_rect(colour = "black", 
                                     fill = NA, 
                                     linewidth = 5), 
         panel.grid.minor = element_blank(),
         axis.title=ggplot2::element_text(size = 16), 
         axis.text=ggplot2::element_text(size = 25),
         legend.text=ggplot2::element_text(size = 14))
}

# View the plot
panelA

# Save the plot
# IMPORTANT: change the file designation to match your machine
ggsave(filename = '~/Documents/MASTERS/Important/map_panelA.png',
       plot = panelA,
       width = 15,
       height = 15)

# Map Section B (Cayman Islands) -----------------------------------------------
lons <- c(-80.55, -79.35)
lats <- c(19.5, 19.9)

spags <- data.frame(x = c(-80.12083),
                    y = c(19.64998),
                    site = c("LCWE"),
                    shp = c("A"))

labs <- data.frame(x = c(-80.27, -79.7),
                   y = c(19.7, 19.8),
                   region = c('Little Cayman', 'Cayman Brac'))

xbreaks <- seq(lons[1], lons[2], 0.6)

xlabels <- paste0(xbreaks, 'º', 'W')

ybreaks <- seq(lats[1], lats[2], 0.2)

ylabels <- paste0(ybreaks, 'º', 'N')

# Make the plot
panelB <- {
   ggplot(carib_bathy, aes(x = x, y = y)) + 
   coord_quickmap() +
   geom_polygon(data = worldLLfull, 
                aes(x = X, 
                    y = Y, 
                    group = PID), 
                fill = "grey", 
                color = "black") +
   geom_point(data = spags, 
              shape = 18, 
              size = 6, 
              fill = 'black', 
              color = 'black') +
   scale_size_manual(values = c(4, 2.5), 
                     guide = 'none') +
   scale_shape_manual(values = c(17, 19), 
                      guide = 'none') +
   geom_text(data = labs, 
             aes(label = region), 
             size = 12) +
   xlab("") +
   ylab("") +
   scale_x_continuous(limits = lons, 
                      position = 'bottom', 
                      breaks = xbreaks, 
                      labels = xlabels) +
   scale_y_continuous(limits = lats, 
                      position = "left", 
                      breaks = ybreaks, 
                      labels = ylabels) +  
   theme_bw()+
   theme(panel.border = element_rect(colour = "black", 
                                     fill = NA, 
                                     linewidth = 5), 
         panel.grid.minor = element_blank(),
         axis.title=ggplot2::element_text(size = 16), 
         axis.text=ggplot2::element_text(size = 25),
         legend.text=ggplot2::element_text(size = 14))
}

# View the plot
panelB

# Save the plot
# IMPORTANT: change the file designation to match your machine
ggsave(filename = '~/Documents/MASTERS/Important/map_panelB.png',
       plot = panelB,
       width = 16,
       height = 8)

# Map Section C (Little Cayman west end shelf edge, array) ---------------------

lcb.ll <- projectRaster(lcbathy, 
                        crs = CRS('+proj=longlat +ellps=WGS84'))
depths <- c(-40, -30, -20, -10, -0.001, Inf)
tmp <- polys <- polys.ll <- list()

for(d in 1:(length(depths) - 1)) {
   tmp[[d]] <- cut(lcbathy, 
                   breaks = c(depths[d], depths[d + 1]))
   polys[[d]] <- rasterToPolygons(tmp[[d]], function (x) {x == 1}, 
                                  dissolve = T) %>% st_as_sf()
}

polys.ll <- lapply(polys, 
                   st_transform, 
                   crs = CRS('+proj=longlat +ellps=WGS84')@projargs)
polys.ll <- lapply(polys.ll, 
                   smoothr::smooth, 
                   method = 'ksmooth', 
                   smoothness = 3)

polys.ll[[1]] <- polys.ll[[1]] %>% mutate(aux = c('1'))
polys.ll[[2]] <- polys.ll[[2]] %>% mutate(aux = c('2'))
polys.ll[[3]] <- polys.ll[[3]] %>% mutate(aux = c('3'))
polys.ll[[4]] <- polys.ll[[4]] %>% mutate(aux = c('4'))
polys.ll[[5]] <- polys.ll[[5]] %>% mutate(aux = c('5'))

df <- do.call(rbind, polys.ll)

hydro_coordinates$NUMBER <- c(1:6)
hydro_coordinates$TYPE <- c('DSG', 'ST', 'DSG', 'ST', 'DSG', 'ST')
hydro_coordinates$LABEL <- c('LS1', 'ST2', 'LS3', 'ST4', 'LS5', 'ST6')

spag_obs <-  c(-80.122850, 19.652769)
spag_hist <- c(-80.12120, 19.65047)
hydro_coordinates[7,1] <- spag_obs[1]
hydro_coordinates[7,2] <- spag_obs[2]
hydro_coordinates[8,1] <- spag_hist[1]
hydro_coordinates[8,2] <- spag_hist[2]

blues <- brewer.pal(n = 9, name = 'Blues')
blues <- blues[c(2, 4, 6, 8)]
colors <- c(blues, '#000000')

xlims <- c(-80.125, -80.109)
ylims <- c(19.648, 19.663)

# Make the plot
panelC <- {
   ggplot(df) +
   geom_sf(aes(geometry = geometry, 
               fill = aux), 
           color = 'black') +
   scale_fill_manual(values = colors, 
                     labels = c('30-40m', '20-30m', '10-20m', 
                                '0-10m', 'Land')) +
   scale_shape(name = '') +
   guides(fill = guide_legend(reverse = T)) +
   geom_point(data = hydro_coordinates[c(1, 3, 5), ], 
              aes(x = LONGITUDE, 
                  y = LATITUDE), 
              shape = 24, 
              size = 12, 
              fill = 'white', 
              color = 'black') +
   geom_point(data = hydro_coordinates[c(2, 4, 6), ], 
              aes(x = LONGITUDE, y = LATITUDE), 
              shape = 24, 
              size = 12, 
              fill = 'black', 
              color = 'black') +
   geom_point(data = hydro_coordinates[7, ], 
              aes(x = LONGITUDE, y = LATITUDE), 
              shape = 4, 
              size = 12, 
              col = 'darkred',
              stroke = 2.5) +
   geom_point(data = hydro_coordinates[8, ], 
              aes(x = LONGITUDE, y = LATITUDE), 
              shape = 8, 
              size = 12, 
              col = 'darkred',
              stroke = 2.5) +
   geom_label(data = hydro_coordinates, 
              aes(x = LONGITUDE, y = LATITUDE, label = LABEL), 
              hjust = -0.6, 
              vjust = -0.3, 
              alpha = 0.7,
              size = 7.8) +
   geom_text(data = data.frame(x = -80.109, 
                               y = 19.661, 
                               label = 'Little\nCayman'), 
             aes(x = x, y = y, label = label), 
             color = 'white', 
             size = 10) +
   coord_sf(xlim = xlims, 
            ylim = ylims) +
   scale_x_longitude(xmin = -80.125, 
                     xmax = -80.110, 
                     step = 0.005) +
   theme_bw() +
   scale_y_latitude(ymin = 19.65, 
                    ymax = 19.66, 
                    step = 0.01) +
   annotation_scale(width_hint = 0.4,
                    height = unit(1, 'cm'),
                    text_cex = 2,
                    pad_x = unit(1, 'cm'),
                    pad_y = unit(1, 'cm')) +
   annotation_north_arrow(location = 'tl', 
                          which_north = 'true',
                          style = 
                             north_arrow_fancy_orienteering(text_size = 30),
                          height = unit(5, 'cm'),
                          width = unit(5, 'cm'),
                          pad_x = unit(1, 'cm'),
                          pad_y = unit(1, 'cm')) +
   theme(legend.position = c(.14, .66),
         legend.title = element_blank(),
         panel.border = element_rect(colour = "black", 
                                     fill = NA, 
                                     linewidth = 4), 
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 16), 
         axis.text = element_text(size = 25),
         legend.text = element_text(size = 20)) +
   guides(fill = guide_legend(keywidth = 0.5,
                              keyheight = 0.5,
                              default.unit = 'inch',
                              reverse = T))
}

# View the plot
panelC

# Save the plot
ggsave(filename = 'map_panelC.png',
       plot = panelC,
       width = 16,
       height = 15)

# use photoshop to add panel labels (a, b, c) & lines

###########################################
##### FIGURE 2: Spatiotemporal Mosaic #####
###########################################
# Data formatting --------------------------------------------------------------
array_data$H4_AAS_COUNT[160] <- NA
array_data$H4_AAS_COUNT[160] <- NA


array_data$HOUR_OF_DAY <- as.factor(array_data$HOUR_OF_DAY)
array_data$DAFS <- as.factor(array_data$DAFS)

# normalized vocalizations per effort
array_data$H1_NORMALIZED_AAS <- array_data$H1_AAS_PER_MIN * 60
array_data$H2_NORMALIZED_AAS <- array_data$H2_AAS_PER_MIN * 60
array_data$H4_NORMALIZED_AAS <- array_data$H4_AAS_PER_MIN * 60
array_data$H5_NORMALIZED_AAS <- array_data$H5_AAS_PER_MIN * 60
array_data$H6_NORMALIZED_AAS <- array_data$H6_AAS_PER_MIN * 60

array_data$H4_NORMALIZED_AAS[160] <- NA

array_data$DATE_TIME <- as.POSIXct(array_data$DATE_TIME, 
                                   tz = 'EST',
                                   format = '%m/%d/%y %H:%M')

# Panel A (Box plots of hours by hydrophone) -----------------------------------
colors <- brewer.pal(n = 11, name = 'BrBG') # 1:5 brown, 6 white, 7:11 bluegreen
# LS1
# Make the plot
h1_boxplot_hr <- {
   ggplot(array_data, 
          aes(x = HOUR_OF_DAY, 
              y = H1_NORMALIZED_AAS)) + 
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) + 
   labs(title = 'A) Hour', 
        x = element_blank(), 
        y = '') +
   theme_linedraw() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank()) +
   theme(axis.text.x = element_blank(),
         axis.ticks = element_blank()) +
   geom_rect(data = array_data, 
             aes(xmin = 0,
                 xmax = 7,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_rect(data = array_data, 
             aes(xmin = 7,
                 xmax = 8,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 19,
                 xmax = 20,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 20,
                 xmax = 25,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) +
   annotate(geom = 'text', 
            x = 3, 
            y = 80, 
            label = 'LS1', 
            color = 'white',
            cex = 8)
}

# View the plot
h1_boxplot_hr


# ST2
# Make the plot
h2_boxplot_hr <- {
   ggplot(array_data, 
          aes(x = HOUR_OF_DAY, 
              y = H2_NORMALIZED_AAS)) + 
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = '') +
   theme_linedraw() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank()) +
   theme(axis.text.x = element_blank(),
         axis.ticks = element_blank()) +
   geom_rect(data = array_data, 
             aes(xmin = 0,
                 xmax = 7,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_rect(data = array_data, 
             aes(xmin = 7,
                 xmax = 8,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 19,
                 xmax = 20,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 20,
                 xmax = 25,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) +
   annotate(geom = 'text', 
            x = 3, 
            y = 295, 
            label = 'ST2', 
            color = 'white',
            cex = 8) +
   ylim(0, 300)
}

# View the plot
h2_boxplot_hr


# ST4
# Make the plot
h4_boxplot_hr <- {
   ggplot(array_data, 
          aes(x = HOUR_OF_DAY, 
              y = H4_NORMALIZED_AAS)) + 
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = '') +
   theme_linedraw() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank()) +
   theme(axis.text.x = element_blank(),
         axis.ticks = element_blank()) +
   geom_rect(data = array_data, 
             aes(xmin = 0,
                 xmax = 7,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_rect(data = array_data, 
             aes(xmin = 7,
                 xmax = 8,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 19,
                 xmax = 20,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 20,
                 xmax = 25,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) +
   annotate(geom = 'text', 
            x = 3, 
            y = 320, 
            label = 'ST4', 
            color = 'white',
            cex = 8) 
}

# View the plot
h4_boxplot_hr


# LS5
h5_boxplot_hr <- {
   ggplot(array_data, 
          aes(x = HOUR_OF_DAY, 
              y = H5_NORMALIZED_AAS)) + 
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = '') +
   theme_linedraw() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank()) +
   theme(axis.text.x = element_blank(),
         axis.ticks = element_blank()) +
   geom_rect(data = array_data, 
             aes(xmin = 0,
                 xmax = 7,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_rect(data = array_data, 
             aes(xmin = 7,
                 xmax = 8,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 19,
                 xmax = 20,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 20,
                 xmax = 25,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) +
   annotate(geom = 'text', 
            x = 3, 
            y = 315, 
            label = 'LS5', 
            color = 'white', 
            cex = 8) 
}

# View the plot
h5_boxplot_hr


# ST6
h6_boxplot_hr <- {
   ggplot(array_data, 
          aes(x = HOUR_OF_DAY, 
              y = H6_NORMALIZED_AAS)) + 
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = 'Hour of Day', 
        y = '') +
   theme_linedraw() + 
   geom_rect(data = array_data, 
             aes(xmin = 0,
                 xmax = 7,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_rect(data = array_data, 
             aes(xmin = 7,
                 xmax = 8,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 19,
                 xmax = 20,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey75',
             fill = 'grey75') +
   geom_rect(data = array_data, 
             aes(xmin = 20,
                 xmax = 25,
                 ymin = -Inf,
                 ymax = Inf),
             alpha = 0.5,
             color = 'grey50',
             fill = 'grey50') +
   geom_boxplot(na.rm = T, 
                outlier.shape = NA) +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.text.x = element_text(angle = 45)) +
   annotate(geom = 'text', 
            x = 3, 
            y = 325, 
            label = 'ST6', 
            color = 'white', 
            cex = 8) +
   scale_x_discrete(breaks = seq(0, 24, 2))
}

# View the plot
h6_boxplot_hr


# Panel B (Box plots of DAFS by hydrophone) ------------------------------------
# LS1
# Make the plot
h1_boxplot_DAFS <- {
   ggplot(array_data, 
          aes(x = DAFS, 
              y = H1_NORMALIZED_AAS)) + 
   geom_boxplot(fill = colors[5], 
                na.rm = T, 
                outlier.shape = NA) + 
   labs(title = 'B) Day', 
        x = element_blank(), 
        y = element_blank()) +
   theme_minimal() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank()) +
   geom_point(data = data.frame(x = 2, 
                                y = 70), 
              aes(x, y), 
              size = 5, 
              color = 'black') +
   geom_point(data = data.frame(x = 2,
                                y = 70), 
              aes(x, y), 
              size = 4, 
              color = 'yellow')
                    }

# View the plot
h1_boxplot_DAFS


# ST2
# Make the plot
h2_boxplot_DAFS <- {
   ggplot(array_data, 
          aes(x = DAFS, 
              y = H2_NORMALIZED_AAS, 
              ymax = 30)) + 
   geom_boxplot(fill = colors[5], 
                na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = element_blank()) +
   theme_minimal() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank()) +
   geom_point(data = data.frame(x = 2, 
                                y = 150), 
              aes(x, y), 
              size = 5, 
              color = 'black') +
   geom_point(data = data.frame(x = 2, 
                                y = 150), 
              aes(x, y), 
              size = 4, 
              color = 'yellow') +
   ylim(0, 175)
                    }

# View the plot
h2_boxplot_DAFS


# ST4
# Make the plot
h4_boxplot_DAFS <- {
   ggplot(array_data, 
          aes(x = DAFS, 
              y = H4_NORMALIZED_AAS)) + 
      geom_boxplot(fill = colors[5], 
                   na.rm = T, 
                   outlier.shape = NA) + 
      labs(title = element_blank(), 
           x = element_blank(), 
           y = element_blank()) + 
      theme_minimal() + 
      geom_point(data = data.frame(x = 2, y = 210), 
                 aes(x, y), 
                 size = 5, 
                 color = 'black') +
      geom_point(data = data.frame(x = 2, y = 210), 
                 aes(x, y), 
                 size = 4, 
                 color = 'yellow') +
      geom_rect(data = array_data, 
                aes(xmin = 8,
                    xmax = 11,
                    ymin = -Inf,
                    ymax = Inf),
                alpha = 0.5,
                color = 'black',
                fill = 'black') +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", 
                                            size = 0.5),
            axis.text.x = element_blank()) +
      ylim(0, 250)
                    }

# View the plot
h4_boxplot_DAFS


# LS5
# Make the plot
h5_boxplot_DAFS <- {
   ggplot(array_data, 
          aes(x = DAFS, 
              y = H5_NORMALIZED_AAS)) + 
   geom_boxplot(fill = colors[5], 
                na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = element_blank()) +
   theme_minimal() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank()) +
   geom_point(data = data.frame(x = 2, y = 190), 
              aes(x, y), 
              size = 5, 
              color = 'black') +
   geom_point(data = data.frame(x = 2, y = 190), 
              aes(x, y), 
              size = 4, 
              color = 'yellow') +
   ylim(0, 225)
                    }

# View the plot
h5_boxplot_DAFS


# ST6
# Make the plot
h6_boxplot_DAFS <- {
   ggplot(array_data, 
          aes(x = DAFS, 
              y = H6_NORMALIZED_AAS)) + 
   geom_boxplot(fill = colors[5], 
                na.rm = T, 
                outlier.shape = NA) + 
   labs(title = element_blank(), 
        x = 'Days Before or After First Spawn', 
        y = element_blank()) +
   theme_minimal() + 
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_text(angle = 45)) +
   geom_point(data = data.frame(x = 2, y = 255), 
              aes(x, y), 
              size = 5, 
              color = 'black') +
   geom_point(data = data.frame(x = 2, y = 255), 
              aes(x, y), 
              size = 4, 
              color = 'yellow') +
   ylim(0, 300)
                    }

# View the plot
h6_boxplot_DAFS

# Panel C (Time series with splines by hydrophone) -----------------------------
# LS1
# Make the plot
h1_time_series <- {
   ggplot(array_data, 
          aes(x = DATE_TIME, 
              y = H1_NORMALIZED_AAS)) +
   geom_point(colour = colors[2], 
              size = 0.5) +
   geom_smooth(n = 1000, 
               level = 0.89, 
               method = 'loess', 
               span = 0.1,
               formula = 'y ~ x', 
               colour = colors[10], 
               size = 0.5) + 
   labs(title = 'C) Continuous Time', 
        x = element_blank(), 
        y = element_blank()) +
   geom_vline(
      xintercept = 
         array_data$DATE_TIME[as.numeric(which(array_data$HOUR_OF_DAY == 0))], 
      linetype = 3) +
   theme_minimal() +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank())
                   }

# View the plot
h1_time_series


# ST2
# Make the plot
h2_time_series <- {
   ggplot(array_data, 
          aes(x = DATE_TIME, 
              y = H2_NORMALIZED_AAS)) +
   geom_point(colour = colors[2], 
              size = 0.5) +
   geom_smooth(n = 1000, 
               level = 0.89, 
               method = 'loess', 
               span = 0.1,
               formula = 'y ~ x', 
               colour = colors[10], 
               size = 0.5) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = element_blank ()) + 
   geom_vline(
      xintercept = 
         array_data$DATE_TIME[as.numeric(which(array_data$HOUR_OF_DAY == 0))], 
      linetype = 3) +
   theme_minimal() +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank())
}

# View the plot
h2_time_series


# ST4
# Make the plot
h4_time_series <- {
   ggplot(array_data[2:159, ], 
          aes(x = 2:159, 
              y = H4_NORMALIZED_AAS)) +
   xlim(1, 218) +
   geom_point(colour = colors[2], 
              size = 0.5) +
   geom_smooth(n = 1000, 
               level = 0.89, 
               method = 'loess', 
               span = 0.1,
               formula = 'y ~ x', 
               colour = colors[10], 
               size = 0.5) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = element_blank()) + 
   geom_vline(
      xintercept = as.numeric(which(array_data$HOUR_OF_DAY == 0)), 
      linetype = 3) +
   theme_minimal() + 
   geom_rect(aes(ymin = -Inf,
                 ymax = Inf,
                 xmin = 159,
                 xmax = Inf),
             colour = 'black', 
             fill = 'black') +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank())
                   }

# View the plot
h4_time_series


# LS5
# Make the plot
h5_time_series <- {
   ggplot(array_data, 
          aes(x = DATE_TIME, 
              y = H5_NORMALIZED_AAS)) +
   geom_point(colour = colors[2], 
              size = 0.5) +
   geom_smooth(n = 1000, 
               level = 0.89, 
               method = 'loess', 
               span = 0.1,
               formula = 'y ~ x', 
               colour = colors[10], 
               size = 0.5) + 
   labs(title = element_blank(), 
        x = element_blank(), 
        y = element_blank()) + 
   geom_vline(
      xintercept = 
         array_data$DATE_TIME[as.numeric(which(array_data$HOUR_OF_DAY == 0))], 
              linetype = 3) +
   theme_minimal() +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_blank())
                   }

# View the plot
h5_time_series


# ST6
# Make the plot
h6_time_series <- {
   ggplot(array_data, 
          aes(x = DATE_TIME, 
              y = H6_NORMALIZED_AAS)) +
   geom_point(colour = colors[2], 
              size = 0.5) +
   geom_smooth(n = 1000, 
               level = 0.89, 
               method = 'loess', 
               span = 0.1,
               formula = 'y ~ x', 
               colour = colors[10], 
               size = 0.5) + 
   labs(title = element_blank(), 
        x = 'Time', 
        y = element_blank()) + 
   geom_vline(
      xintercept = 
         array_data$DATE_TIME[as.numeric(which(array_data$HOUR_OF_DAY == 0))], 
              linetype = 3) +
   theme_minimal() +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", 
                                         size = 0.5),
         axis.text.x = element_text(angle = 45, vjust = 0.7))
                   }

# View the plot
h6_time_series

# Final Mosaic -----------------------------------------------------------------
# Make the mosaic
mosaic <- (h1_boxplot_hr | h1_boxplot_DAFS | h1_time_series) / 
          (h2_boxplot_hr | h2_boxplot_DAFS | h2_time_series) / 
          (h4_boxplot_hr | h4_boxplot_DAFS | h4_time_series) /
          (h5_boxplot_hr | h5_boxplot_DAFS | h5_time_series) /
          (h6_boxplot_hr | h6_boxplot_DAFS | h6_time_series)

# Add y axis label
mosaic <- wrap_elements(mosaic) +
            labs(tag = expression('Effort Normalized AAS Counts h'^-1)) +
            theme(plot.tag = element_text(size = rel(1.5),
                                          angle = 90),
                  plot.tag.position = 'left'
                  )

# View the mosaic
mosaic

# Save the mosaic
ggsave(filename = 'spatiotemporal_mosaic.png',
       plot = mosaic,
       width = 12,
       height = 15)

############################################
##### FIGURE 3: Correlation Pairs Plot #####
############################################
# extract posterior from the model
post <- extract.samples(mod)

# color
colors <- brewer.pal(n = 11, name = 'BrBG') 
# 1:5 brown, 6 white, 7:11 bluegreen

# find the median posterior of the covariance estimate
hydro_cor <- apply(post$Rho_d, c(2, 3), median) 

# attach hydrophone stations
colnames(hydro_cor) <- c("LS1", "ST2", "ST4", "LS5", "ST6")
rownames(hydro_cor) <- c("LS1", "ST2", "ST4", "LS5", "ST6")

# Set parameters to save the plot
pdf('cor_matrix_plot.pdf',
    width = 8,
    height = 7)

# Make / View the plot
corrplot(hydro_cor, 
         method = c('ellipse'),
         tl.srt = 0, 
         tl.col = 'black',
         tl.offset = 0.9,
         tl.cex = 1,
         col = COL2('BrBG', n = 100),
         outline = 'black',
         mar = c(0, 0, 2, 1) # bottom, left, top, right
) 
text(x = 6.15, 
     y = 3, 
     label = expression(paste('Correlation (',rho,')')), 
     srt = 270,
     cex = 1)

# Close the plot save window
dev.off()

#############################################
##### FIGURE 4: Correlation v. Distance #####
#############################################
# Building the dataset ---------------------------------------------------------
# Extract 25% and 75% quantiles for correlation posterior median
cor_quantile <- apply(post$Rho_d, c(2,3), quantile)

# 25% quantiles
hold <- vector()
quantile_25 <- vector()
n <- 1
for (i in 1:5) {
   m <- n + 4
   for (j in 1:5) {
      hold[j] <- cor_quantile[2, i, j]
   }
   
   quantile_25[n:m] <- hold
   n <- n + 5
}

# 75% quantiles
hold <- vector()
quantile_75 <- vector()
n <- 1
for(i in 1:5) {
   m <- n + 4
   for(j in 1:5) {
      hold[j] <- cor_quantile[4, i, j]
   }
   quantile_75[n:m] <- hold
   n <- n + 5
}

# build the dataset
cor_dist_data <- 
   data.frame(CORRELATION_COEFFICIENT = as.vector(hydro_cor),
              DISTANCE_BETWEEN_PAIRS = as.vector(hydro_distance_matrix))
cor_dist_data$COR_QUANTILE_25 <- quantile_25
cor_dist_data$COR_QUANTILE_75 <- quantile_75
cor_dist_data$HYDRO_PAIRS <- c('LS1-LS1', 'LS1-ST2', 'LS1-ST4', 'LS1-LS5', 
                               'LS1-ST6', 'ST2-LS1', 'ST2-ST2', 'ST2-ST4', 
                               'ST2-LS5', 'ST2-ST6', 'ST4-LS1', 'ST4-ST2', 
                               'ST4-ST4', 'ST4-LS5', 'ST4-ST6', 'LS5-LS1', 
                               'LS5-ST2', 'LS5-ST4', 'LS5-LS5', 'LS5-ST6', 
                               'ST6-LS1', 'ST6-ST2', 'ST6-ST4', 'ST6-LS5', 
                               'ST6-ST6')

# remove direct pairwise comparisons (e.g. LS1-LS1)
cor_dist_data <- 
   cor_dist_data[!duplicated(cor_dist_data$CORRELATION_COEFFICIENT),]
cor_dist_data <- cor_dist_data[-1, ]

# Modeling correlation v distance ----------------------------------------------
# format data as a list for rethinking model
cor_dist_list <- list(
   Rho = cor_dist_data$CORRELATION_COEFFICIENT,
   dist = cor_dist_data$DISTANCE_BETWEEN_PAIRS / 
          max(cor_dist_data$DISTANCE_BETWEEN_PAIRS)
                      )

# set seed for duplication
set.seed(472787)

# run the model
cor_dist_mod <- quap(
   alist(
      Rho ~ dnorm(mu, sigma),
      mu <- a + bd * dist,
      a ~ dnorm(0.9, 0.05),
      bd ~ dnorm(0, 0.5),
      sigma ~ dexp(1)
   ), data = cor_dist_list
                     )

# extract posterior estimates, attaching an arbitrary distance seq
distance_seq <- seq(from = 0, to = 1, length.out = 100)
mu <- link(cor_dist_mod, 
           data = list(dist = distance_seq))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI, prob = 0.95)

# Plotting correlation v distance ----------------------------------------------
# colors for plot
colors <- brewer.pal(n = 11, name = 'BrBG') # 1:5 brown, 6 white, 7:11 bluegreen

# Set parameters to save the plot
pdf('cor_dist_plot.pdf',
    width = 11,
    height = 7)

# create plot frame
plot(NULL, 
     xlab = '', 
     xaxt = 'n', 
     ylab = '', 
     yaxt = 'n', 
     ylim = c(-0.6, 1), 
     xlim = c(0, 1))

# add model predicted line
lines(distance_seq, mu_mean, lwd = 2)

# add model predicted shade (95% confidence)
shade(mu_PI, distance_seq)

# overlay new plot of descriptive axes
par(new = T)
plot(NULL,  
     xlab = 'Distance between pairs (m)', 
     ylab = expression(paste('Correlation (',rho,')')), 
     xlim = c(0, 650), 
     ylim = c(-0.6, 1))

# add dashed line at y = 0
abline(h = 0, lty = 2)

# add lines of correlation 25% and 75% quantiles
for (i in 1:10) {
   lines(c(cor_dist_data$DISTANCE_BETWEEN_PAIRS[i], 
           cor_dist_data$DISTANCE_BETWEEN_PAIRS[i]),
         c(cor_dist_data$COR_QUANTILE_25[i], 
           cor_dist_data$COR_QUANTILE_75[i]), 
         col = 'black')
}

# add points
points(cor_dist_data$DISTANCE_BETWEEN_PAIRS, 
       cor_dist_data$CORRELATION_COEFFICIENT, 
       col = colors[8], 
       pch = 19, 
       cex = 1.5)
points(cor_dist_data$DISTANCE_BETWEEN_PAIRS, 
       cor_dist_data$CORRELATION_COEFFICIENT, 
       col = 'black', 
       cex = 1.5)

# add labels to points of interest (click on plot)
   # the labels will adjust to side of point you click on 
for (i in 1:nrow(cor_dist_data)) {
   if (cor_dist_data$DISTANCE_BETWEEN_PAIRS[i] < 200) {
      
      text(x = cor_dist_data$DISTANCE_BETWEEN_PAIRS[i] + 35,
           y = cor_dist_data$CORRELATION_COEFFICIENT[i],
           labels = cor_dist_data$HYDRO_PAIRS[i])
   
   } else if (cor_dist_data$DISTANCE_BETWEEN_PAIRS[i] > 200 &
              cor_dist_data$DISTANCE_BETWEEN_PAIRS[i] < 400 &
              cor_dist_data$CORRELATION_COEFFICIENT[i] > 0) {
      
      text(x = cor_dist_data$DISTANCE_BETWEEN_PAIRS[i] + 35,
           y = cor_dist_data$CORRELATION_COEFFICIENT[i],
           labels = cor_dist_data$HYDRO_PAIRS[i])
      
   } else {
      
      text(x = cor_dist_data$DISTANCE_BETWEEN_PAIRS[i] - 35,
           y = cor_dist_data$CORRELATION_COEFFICIENT[i],
           labels = cor_dist_data$HYDRO_PAIRS[i])
      
   }
}

# close the plot save window
dev.off()




###########################################
##### FIGURE 5: Model Posterior Plots #####
###########################################
# Building the datasets --------------------------------------------------------
# extract posterior estimates
post <- extract.samples(mod)

# create posterior prediction datasets by parameter
post_hr <- data.frame(exp(post$bt))
post_dafs <- data.frame(exp(post$bd))
post_fp <- data.frame(exp(post$fp))

# fix column names for pivoting
colnames(post_hr) <- as.character(c(0:23))
colnames(post_dafs) <- as.character(c(-4:2))
colnames(post_fp) <- c('Not \nObserved', 
                       'Nearby \n(20-100 m)', 
                       'Present \n(0-20 m)')

# pivot datasets for ggplot functionality
post_hr <- pivot_longer(post_hr, 
                        cols = everything(), 
                        names_to = 'HOUR_OF_DAY',
                        values_to = 'POSTERIOR')
post_dafs <- pivot_longer(post_dafs,
                          cols = everything(),
                          names_to = 'DAFS',
                          values_to = 'POSTERIOR')
post_fp <- pivot_longer(post_fp,
                        cols = everything(),
                        names_to = 'FISH_PROX',
                        values_to = 'POSTERIOR')

# factor names_to columns
post_hr$HOUR_OF_DAY <- factor(post_hr$HOUR_OF_DAY,
                              levels = c(0:23))
post_dafs$DAFS <- factor(post_dafs$DAFS,
                         levels = c(-4:2))
post_fp$FISH_PROX <- factor(post_fp$FISH_PROX,
                            levels = c('Not \nObserved',
                                       'Nearby \n(20-100 m)',
                                       'Present \n(0-20 m)'))

# Panel A (Box plots of posterior predictions by hour of day) ------------------
# Make the plot
fig5_panelA <- {
   ggplot(data = post_hr,
          aes(x = HOUR_OF_DAY,
              y = POSTERIOR)) +
      geom_rect(aes(ymin = -Inf,
                    ymax = Inf,
                    xmin = -Inf,
                    xmax = 6),
                colour = 'darkgrey', 
                fill = 'darkgrey') +
      geom_rect(aes(ymin = -Inf,
                    ymax = Inf,
                    xmin = 19,
                    xmax = Inf),
                color = 'darkgrey',
                fill = 'darkgrey') +
      geom_rect(aes(ymin = -Inf,
                    ymax = Inf,
                    xmin = 6,
                    xmax = 7),
                color = 'lightgrey',
                fill = 'lightgrey') +
      geom_rect(aes(ymin = -Inf,
                    ymax = Inf,
                    xmin = 18,
                    xmax = 19),
                color = 'lightgrey',
                fill = 'lightgrey') +
      geom_boxplot(outliers = F,
                   fill = colors[9]) +
      scale_y_discrete(breaks = c(2, 4, 6),
                       limits = c(2, 4, 6)) +
      labs(x = 'Hour of Day',
           y = element_blank()) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 15)) 
}

# View the plot
fig5_panelA


# Panel B (Box plots of posterior predictions by DAFS) -------------------------
# Make the plot
fig5_panelB <- {
   ggplot(data = post_dafs,
          aes(x = DAFS,
              y = POSTERIOR)) + 
      geom_boxplot(outliers = F,
                   fill = colors[7]) +
      labs(x = 'Days Before or After First Spawn (DAFS)',
           y = element_blank()) +
      scale_y_discrete(breaks = c(2, 4, 6, 8),
                       limits = c(2, 4, 6, 8)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 15)) +
      geom_point(data = data.frame(x = 2,
                                   y = 7.8),
                 aes(x = x,
                     y = y),
                 color = 'black',
                 size = 11) +
      geom_point(data = data.frame(x = 2,
                                   y = 7.8),
                 aes(x = x,
                     y = y),
                 color = 'yellow',
                 size = 10) 
}

# View the plot
fig5_panelB


# Panel C (Box plots of posterior predictions by fish proximity) ---------------
# Make the plot
fig5_panelC <- {
   ggplot(data = post_fp,
          aes(x = FISH_PROX,
              y = POSTERIOR)) +
      geom_boxplot(outliers = F,
                   fill = colors[4]) +
      labs(x = 'Bulk of Fish Proximity to Hydrophone',
           y = element_blank()) +
      scale_y_discrete(breaks = c(5, 10, 15),
                       labels = c(5, 10, 15),
                       limits = c(5, 10, 15)) +
      expand_limits(y = 23) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 15)) 
}

# View the plot
fig5_panelC
# removed rows are outliers that are not displayed

# Final Mosaic -----------------------------------------------------------------
# Make the shared y-axis label
ylab <- 
   ggplot() +
   annotate(geom = 'text',
            x = 1, 
            y = 1,
            label = expression('Predicted AAS h'^-1),
            angle = 90,
            size = 9) +
   coord_cartesian(clip = 'off') +
   theme_void()

# Make the mosaic
fig5_mosaic <- {
   (ylab |
   (wrap_plots(A = fig5_panelA,
              B = fig5_panelB,
              C = fig5_panelC,
              design = '
                       AAAAA
                       BBBCC
                       ') +
       plot_layout(heights = c(1.5, 2))
    )) +
   plot_layout(widths = c(.05, 1))
}

# View the mosaic
fig5_mosaic 

# Save the mosaic
ggsave(filename = 'posterior_mosaic.png',
       plot = fig5_mosaic,
       width = 12,
       height = 12)

# use photoshop to add panel labels (a, b, c) & fishprox cartoons


##############################################
##### FIGURE S1: Map of Aggregation Site #####
##############################################
# Panel A (Little Cayman bathymetry) -------------------------------------------
# set colors
blues <- brewer.pal(n = 9, name = 'Blues')
blues <- blues[c(2, 4, 6, 8)]
colors <- c(blues, '#000000')

# new limits for coordinates
xlims_A <- c(-80.115, -79.956)
ylims_A <- c(19.640, 19.730)

s1_panelA <- { 
   
   ggplot(df) +
      geom_sf(aes(geometry = geometry, 
                  fill = aux), 
              color = 'black') +
      scale_fill_manual(values = colors, 
                        labels = c('30-40m', '20-30m', 
                                   '10-20m', '0-10m', 
                                   'Land')) +
      scale_shape(name = '') +
      guides(fill = guide_legend(reverse = T)) +
      geom_text(data = data.frame(x = -80.045, 
                                  y = 19.688, 
                                  label = 'Little\nCayman'), 
                aes(x = x, 
                    y = y, 
                    label = label), 
                color = 'white', 
                size = 10) +
      coord_sf(xlim = xlims_A, 
               ylim = ylims_A) +
      scale_x_longitude(xmin = -80.120, 
                        xmax = -79.960, 
                        step = 0.080) + 
      scale_y_latitude(ymin = 19.645, 
                       ymax = 19.725, 
                       step = 0.040) +
      annotation_scale(width_hint = 0.4,
                       text_cex = 2,
                       pad_x = unit(1, 'cm'),
                       pad_y = unit(0.5, 'cm'),
                       location = 'br') + 
      theme_bw() +
      theme(legend.position = c(.74, .3),
            legend.title = element_blank(),
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        linewidth = 4), 
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 16), 
            axis.text = element_text(size = 25),
            legend.text = element_text(size = 20)) +
      guides(fill = guide_legend(keywidth = 0.5,
                                 keyheight = 0.5,
                                 default.unit = 'inch',
                                 reverse = T)) +
      geom_rect(xmin = -80.124, 
                xmax = -80.118, 
                ymin = 19.649, 
                ymax = 19.654, 
                fill = NA, 
                color = 'black', 
                size = 0.75
      ) +
      annotation_north_arrow(location = 'tl', 
                             which_north = 'true',
                             style = 
                                north_arrow_fancy_orienteering(text_size = 30),
                             height = unit(5, 'cm'),
                             width = unit(5, 'cm'),
                             pad_x = unit(1, 'cm'),
                             pad_y = unit(1, 'cm'))
   
}

# View the plot
s1_panelA

# Save the plot
ggsave(filename = 's1_panelA.png',
       plot = s1_panelA,
       width = 16,
       height = 10)

# Panel B (Little Cayman West End Site) ----------------------------------------
xlims_B <- c(-80.124, -80.118)
ylims_B <- c(19.649, 19.654)

s1_panelB <- {
   
   ggplot(df) +
      geom_sf(aes(geometry = geometry, 
                  fill = aux), 
              color = NA, 
              show.legend = F) +
      scale_fill_manual(values = colors, 
                        labels = c('30-40m', '20-30m', 
                                   '10-20m', '0-10m', 'Land')) +
      geom_point(data = landmark_coordinates, 
                 aes(x = LONGITUDE, 
                     y = LATITUDE), 
                 size = 12, 
                 shape = landmark_coordinates$SHAPE, 
                 fill = landmark_coordinates$COLOR, 
                 color = landmark_coordinates$OUTLINE) +
      geom_point(data = landmark_coordinates[14:15, ],
                 aes(x = LONGITUDE,
                     y = LATITUDE),
                 size = 12, 
                 shape = landmark_coordinates$SHAPE[14:15],
                 fill = landmark_coordinates$COLOR[14:15],
                 color = landmark_coordinates$OUTLINE[14:15],
                 stroke = 2.5) + 
      
      # Landmarks
      geom_label(data = landmark_coordinates[7, ], 
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = 1.1, 
                 vjust = 1, 
                 alpha = 0.7,
                 size = 7.8) +
      geom_label(data = landmark_coordinates[8:9, ], 
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = 1.15, 
                 vjust = 0.1, 
                 alpha = 0.7,
                 size = 7.8) +
      
      # LS Hydrophones
      geom_label(data = landmark_coordinates[c(1,5), ], 
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = -0.3, 
                 vjust = -0.3, 
                 alpha = 0.7,
                 size = 7.8) +
      geom_label(data = landmark_coordinates[3, ], 
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = 1.25, 
                 vjust = -0.25, 
                 alpha = 0.7,
                 size = 7.8) +
      
      # Moorings
      geom_label(data = landmark_coordinates[10:13, ],
                 inherit.aes = F,
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = -0.05, 
                 vjust = -0.5, 
                 alpha = 0.7,
                 size = 7.8) +
      
      # ST Hydrophones
      geom_label(data = landmark_coordinates[c(2,4,6), ], 
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = -0.1, 
                 vjust = -0.5, 
                 alpha = 0.7,
                 size = 7.8) +
      
      # FSA
      geom_label(data = landmark_coordinates[c(14,15), ], 
                 aes(x = LONGITUDE, 
                     y = LATITUDE, 
                     label = LABEL), 
                 hjust = -0.15, 
                 vjust = 0.65, 
                 alpha = 0.7,
                 size = 7.8) +
      
      coord_sf(xlim = xlims_B, 
               ylim = ylims_B) +
      scale_x_continuous('') +
      scale_y_continuous('') +
      
      annotation_scale(width_hint = 0.4,
                       text_cex = 2,
                       pad_x = unit(1, 'cm'),
                       pad_y = unit(0.5, 'cm')) + 
      
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        linewidth = 4), 
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 16), 
            axis.text = element_text(size = 25),
            legend.text = element_text(size = 20))
   
}

# View the plot
s1_panelB

# make legend to overlay
table.legend <- data.frame()
table.legend[1:6, 'Labels'] <- c('LS Hydrophones', 'ST Hydrophones', 
                                 'Moorings', 'Landmarks', 
                                 'Historic FSA Center', 'Observed FSA Center')
table.legend[1:6, 'Shapes'] <- c(24, 24, 23, 21, 8, 4)
table.legend[1:6, 'Colors'] <- c('white', 'black', 'green', 'blue', 
                                 'red', 'red')
table.legend[1:6, 'Outline'] <- c(rep('black', 4), rep('red', 2))
table.legend[1:6, 'Stroke'] <- c(rep(0.5, 4), rep(2.5, 2))

level_order <- rev(table.legend$Labels)

legend <- ggplot(table.legend) + 
   geom_point(data = table.legend, aes(x = factor(Labels, 
                                                  level = level_order), 
                                       y = 1.99), 
              shape = table.legend$Shapes, 
              fill = table.legend$Colors,
              color = table.legend$Outline, 
              stroke = table.legend$Stroke,
              size = 7.8) +
   theme_void() + 
   coord_flip() +
   ylim(1.9, 2.1) +
   theme(axis.text.y = element_text(size = 20),
         plot.background = element_rect(fill = 'white',
                                        color = 'black'),
         plot.margin = unit(c(0.3, 0.7, 0.3, 0.9), 'cm')) +
   aes(x = fct_inorder(table.legend$Labels))

legend

# final plot with legend
s1_panelB_final <- s1_panelB + 
   inset_element(p = legend,
                 left = 0.65,
                 bottom = 0.6,
                 right = 0.95,
                 top = 0.9) 

# View the plot
s1_panelB_final

# save the plot
ggsave(filename = 's1_panelB.png',
       plot = s1_panelB_final,
       width = 16,
       height = 16)

####################################################
##### NEW FIGURE: fish proximity to hydrophone #####
####################################################
# Plotting panel A (fish prox) and B (aas rates) in a wrap --------------------- 
# colors for plot
colors <- brewer.pal(n = 11, name = 'BrBG') # 1:5 brown, 6 white, 7:11 bluegreen

# Make panel A
fish_prox_tile <- {
   ggplot(data = visual_acoustic_data, 
          aes(x = DATE_TIME,
              y = HYDROPHONE,
              fill = FISH_PROX)) +
   geom_tile(width = 4000) +
   labs(x = '',
        y = '') +
   scale_fill_gradientn('Fish\nProximity',
                        limits = c(0, 2),
                        breaks = c(0, 1, 2),
                        colors = colors[c(11, 10, 8)]) +
   theme_classic() + 
   annotate(geom = 'rect',
            color = 'black',
            fill = 'black',
            xmin = as.POSIXct('2020-02-12 00:00:00', tz = 'EST'), 
            xmax = as.POSIXct('2020-02-15 00:00:00', tz = 'EST'),
            ymin = 6.5,
            ymax = 6.6) +
   # geom_hline(yintercept = 
   #               fish_prox_data$DATE_TIME[
   #                  which(fish_prox_data$HOUR_OF_DAY == 0)]) +
   # geom_hline(yintercept = 
   #               fish_prox_data$DATE_TIME[
   #                  which(fish_prox_data$HOUR_OF_DAY == 12)],
   #            linetype = 'dashed') +
   theme(axis.line = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_text(size = 15),
         axis.text = element_text(size = 15),
         axis.text.y = element_text(margin = margin(r = -34)),
         axis.text.x = element_blank(),
         legend.position = c(.12, .8),
         legend.spacing = unit(0.5, 'cm'),
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15),
         plot.margin = unit(c(1, 1, -1, 1), 'cm')) +
   guides(fill = guide_legend(byrow = T,
                              reverse = T)) +
      annotate('text',
               x = visual_acoustic_data$DATE_TIME[1240],
               y = 6.2,
               label = '(a)',
               color = 'white',
               size = 12)
}

# View panel A
fish_prox_tile

# Make panel B
aas_rate_tile <- {
   ggplot(data = visual_acoustic_data,
          aes(x = DATE_TIME,
              y = HYDROPHONE,
              fill = AAS_NORMALIZED_RATE)) +
      geom_tile(width = 4000) +
      labs(x = 'Time',
           y = '') +
      scale_fill_gradientn('AAS Detection\nRate',
                           colors = colors[11:6], 
                           breaks = c(100, 200, 300)) +
      theme_classic() +
      annotate(geom = 'rect',
               color = 'black',
               fill = 'black',
               xmin = as.POSIXct('2020-02-12 00:00:00', tz = 'EST'), 
               xmax = as.POSIXct('2020-02-15 00:00:00', tz = 'EST'),
               ymin = 6.5,
               ymax = 6.6) +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 15),
            axis.text.y = element_text(margin = margin(r = -34)),
            axis.text.x = element_text(margin = margin(t = -5)),
            legend.position = c(.83, .5),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            plot.margin = unit(c(-0.2, 1, 1, 1), 'cm')) +
      annotate('text',
               x = visual_acoustic_data$DATE_TIME[1240],
               y = 6.2,
               label = '(b)',
               color = 'white',
               size = 12)
}

# View panel B
aas_rate_tile

# Wrap the plots
heat_map <- (fish_prox_tile / plot_spacer() / aas_rate_tile)  +
   plot_layout(heights = c(1, -0.05, 1))

# View the wrap
heat_map

# Save the wrap
ggsave(filename = 'heat_map.png',
       plot = heat_map,
       width = 12,
       height = 18)


#############################################
##### FIGURE S2: Anthropogenic Activity #####
#############################################
# Make the time series
anthro_plot <- {
   ggplot(data = anthro_data,
          aes(x = DATE_TIME,
              y = INTERFERENCE_PRESENT)) + 
   labs(x = 'Time',
        y = 'Background Interference Present') +
   geom_col(na.rm = T) +
   # geom_vline(xintercept = anthro_data$DATE_TIME[which(anthro_data$TIME == '09:00')],
   #            linewidth = 5,
   #            alpha = 0.5) +
   facet_wrap(vars(HYDROPHONE),
              nrow = 5) + 
   scale_y_continuous(breaks = c(0, 1),
                      limits = c(0, 1)) +
   scale_x_datetime(date_breaks = '1 day') +
   theme_bw() +
   theme(panel.background = element_rect(fill = 'white'), 
         panel.grid = element_blank(),
         strip.background = element_rect(fill = 'black'),
         strip.text = element_text(color = 'white',
                                   size = 15),
         axis.title = element_text(size = 15),
         axis.text = element_text(size = 15))
}

# View the time series
anthro_plot

# Save the time series
ggsave(filename = 'anthro_plot.png',
       plot = anthro_plot,
       width = 15,
       height = 10)


# Make the histogram
anthro_hist <- {
   ggplot(data = anthro_data %>%
             filter(INTERFERENCE_PRESENT == 1),
          aes(x = HOUR)) +
   geom_histogram(stat = 'count',
                  color = 'black',
                  fill = 'black') +
   labs(x = 'Hour of Day',
        y = 'Total Minutes Removed') +
   facet_wrap(vars(HYDROPHONE),
              nrow = 5) +
   coord_flip() +
   scale_x_discrete(drop = F,
                    position = 'top',
                    breaks = c('04', '08', '12', '16', '20')) +
   scale_y_reverse(expand = c(0, 0),
                   breaks = c(50, 100, 150, 200),
                   limits = c(225, 0)) +
   theme_bw() +
   theme(panel.grid.major.y = element_blank(),
         panel.background = element_rect(fill = 'white'),
         strip.background = element_rect(fill = 'black'),
         strip.text = element_text(color = 'white',
                                   size = 15),
         axis.title = element_text(size = 15),
         axis.text = element_text(size = 15))
}

# View the histogram
anthro_hist

# Make a wrap
anthro_wrap <- (anthro_plot | anthro_hist) + 
   plot_layout(widths = c(1, 0.2))

# View the wrap
anthro_wrap

# Save the wrap
ggsave(filename = 'anthro_wrap.png',
       plot = anthro_wrap,
       width = 17,
       height = 10)

