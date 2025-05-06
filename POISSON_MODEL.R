###########################################################################
######## Time regression analysis for count of  dolphin click trains ##########
###########################################################################

# Clean:
rm(list = ls(all = TRUE)) #clear all
graphics.off() # close all
gc() # Clear memory

# Load data:
data_clx <- read.csv("count_dol_trns.csv")


# #Reduce the number of observations to check if the model works
# index_to_take_out <- runif(n = round(nrow(data_clx)*0.9), min = 1, max = nrow(data_clx))
# data_clx <- data_clx[-index_to_take_out, ] 

# Check the response variable:
x11()
hist(data_clx$Noftrns)

#Plot observed counts for daytime
require(ggplot2)
count_trns <- ggplot() +
  geom_point(data = data_clx,
             aes(x = HORA,
                 y = Noftrns,
                 fill = SITIO,
                 color = SITIO),
             shape = 21)
x11()
count_trns



# Bayesian analysis: MCMC sampling of posterior with JAGS software
# model in JAGS code
sink("model_dolph.txt")
cat("model
{
  
  # Likelihood:
  for (i in 1:n_count) {
    # Model:
    count[i] ~ dpois(lambda[i])  # Random effect of population on the mean counts
    
    # Log link function:
    log(lambda[i]) <- log_lambda[i]
    
    # As a function of time:
    log_lambda[i] <- alpha[sitio[i]] + beta[sitio[i]]*hour_sc[i] + beta2[sitio[i]]*pow(hour_sc[i], 2) + beta3[sitio[i]]*pow(hour_sc[i], 3)
  }
  
   # Hyper-distributions:
  for (j in 1:n_sitio) {
   alpha[j] ~ dnorm(hypmu_alpha, 1/pow(hypsd_alpha, 2))
   beta[j] ~ dnorm(hypmu_beta, 1/pow(hypsd_beta, 2))
   beta2[j] ~ dnorm(hypmu_beta2, 1/pow(hypsd_beta2, 2))
   beta3[j] ~ dnorm(hypmu_beta3, 1/pow(hypsd_beta3, 3))
   
  }
  
  
   # Hyper-priors:
  hypmu_alpha ~ dunif(0, 10)
  hypmu_beta ~ dunif(-10, 5)
  hypmu_beta2 ~ dunif(-10, 10)
  hypmu_beta3 ~ dunif(-5, 5)
  
  hypsd_alpha ~ dunif(0, 10)
  hypsd_beta ~ dunif(0, 10)
  hypsd_beta2 ~ dunif(0, 10)
  hypsd_beta3 ~ dunif(0, 5)
  

  
  # Other predictions:
  for (k in 1:n_predict) {
  
    # Hyper-preditions:
    hyper_lambda[k] <- exp(hypmu_alpha + hypmu_beta*hour_sc_pred[k] + hypmu_beta2*pow(hour_sc_pred[k], 2) + hypmu_beta3*pow(hour_sc_pred[k], 3))
    
    # Estimate the residuals:
    residuals[k] <- count[k] - hyper_lambda[k]
    
    # Group predictions:
    for (l in 1:n_sitio) {
      lambda_sitio[k, l] <- exp(alpha[l] + beta[l]*hour_sc_pred[k] + beta2[l]*pow(hour_sc_pred[k], 2) + beta3[l]*pow(hour_sc_pred[k], 3))
    }
  }
  
  # Bayesian R-squared (a proportion of variance explained by the predictions):
  bay_r_sq <- pow(sd(hyper_lambda[]), 2) / (pow(sd(hyper_lambda[]), 2) + pow(sd(residuals[]), 2))
  
 
  }

",fill = TRUE)
sink()
# end of JAGS model code


#Data for the model as a list:
  data.list <- list(
    # Sample sizes:
    n_count = nrow(data_clx),
    n_sitio = length(unique(data_clx$SITIO)),
    n_predict = length(unique(data_clx$HORA)),
    
    # Counts:
    count = data_clx$Noftrns,
    sitio = as.numeric(as.factor(data_clx$SITIO)),
    
    # Predictor:
    hour_sc = (data_clx$HORA - mean(data_clx$HORA))/sd(data_clx$HORA),
    
    # For predictions:
    hour_sc_pred = unique((data_clx$HORA - mean(data_clx$HORA))/sd(data_clx$HORA))
  )
# specify which stochastic quantities ("parameters") to monitor
  monitor <- c(
    "hypmu_alpha",
    "hypmu_beta",
    "hypmu_beta2",
    "hypmu_beta3",
    
    "hypsd_alpha",
    "hypsd_beta",
    "hypsd_beta2",
    "hypsd_beta3",
    
    
    "hyper_lambda",
    "lambda_sitio",
    
    "bay_r_sq")

# MCMC control
n.chains <- 3
n.iter <- 10000
n.burnin <- n.iter*0.2
n.thin <- 20

# Run analysis:
require(jagsUI)
out <- jags(data = data.list,
            # inits = inits.list,
            parameters.to.save = monitor,
            model.file = "model_dolph.txt",
            n.chains = n.chains,
            n.thin = n.thin,
            n.iter = n.iter,
            n.burnin = n.burnin,
            DIC = TRUE,
            parallel = TRUE
)

#### SUMMARIZE RESULTS: ####
require(MCMCvis)
model_summary <- MCMCsummary(object = out, 
                             params = "all",
                             probs = c(0.025, 0.125, 0.25, 
                                       0.5,
                                       0.75, 0.875, 0.975),
                             round = 3,
                             func = function(x) ecdf(x)(0)) # the fraction of observations less or equal to (t)


# Add the n.eff in %:
model_summary$n.eff <- do.call("c", out$n.eff)
model_summary$n.eff.perc <- round((do.call("c", out$n.eff)/out$mcmc.info$n.samples)*100)


# #### VISUALIZE POSTERIORS: ####
MCMCtrace(out, 
          params = 'all',
          iter = n.iter,
          Rhat = TRUE,
          n.eff = FALSE,
          ind = TRUE,
          type = 'density',
          # priors = priors_matrix,
          open_pdf = TRUE)


###########################################################################
###########################################################################
######## START THE DIAGNOSTICS FOR THE MCMC  ##############################
###########################################################################
###########################################################################

# Libraries for diagnostics:
require(mcmcplots)
require(ggmcmc)  # For MCMC diagnostics
require(coda)    # For MCMC analysis
require(lattice) # For quick posterior plotting and diagnostics

# Extract MCMC-type list: 
mcmc_list <- as.mcmc(do.call("cbind", out$sims.list)) # Column-bind all parameters

# Autocorrelation plots:
x11()
autocorr.plot(mcmc_list)

x11()
acfplot(mcmc_list)

# Gelman plot (convergence):
x11()
gelman.plot(out$samples) # shrink factor by iteration bins

#Posterior correlations:
chains_res <- ggs(mcmc_list)

# Correlation plot:
x11()
ggs_pairs(chains_res,
          lower = list(continuous = "density"))



# Alternative: bivariate density plots of parameters
library(fields)                        # for functions image.smooth and tim.colors

mcmc_df <- as.data.frame(mcmc_list)

x11()
ncut <- 50                             # number of bins for variables
ndx.par <- 1:ncol(mcmc_df)             # all parameters
ndx <- combn(ndx.par, 2)               # matrix of all combinations of indices to be plotted
par(mfrow = c(floor(sqrt(ncol(ndx))),
              ceiling(ncol(ndx)/floor(sqrt(ncol(ndx))))),
    mar = c(5,5,2,1))


for (k in 1:ncol(ndx)) {
  
  x <- seq(min(mcmc_df[,ndx[1,k]]),
           max(mcmc_df[,ndx[1,k]]),
           length = ncut)
  
  y <- seq(min(mcmc_df[,ndx[2,k]]),
           max(mcmc_df[,ndx[2,k]]),
           length = ncut)
  
  z <- as.matrix(table(cut(mcmc_df[,ndx[1,k]],ncut),
                       cut(mcmc_df[,ndx[2,k]],ncut)))
  
  Pcor <- round(cor(mcmc_df[,ndx[1,k]],
                    mcmc_df[,ndx[2,k]]), 3)
  
  image(x,y,z = image.smooth(z)$z,
        col = tim.colors(64),
        xlab = names(mcmc_df)[ndx[1,k]],
        ylab = names(mcmc_df)[ndx[2,k]],
        main = paste("Pearson cor =", 
                     Pcor))  # smoothed color plot
}


###########################################################################
###########################################################################
######## END OF THE DIAGNOSTICS FOR THE MCMC  #############################
###########################################################################
###########################################################################

# Build the dataframe for predictions:

# hypermean predicitons:
hy_lambda_dframe <- model_summary[grep("hyper_lambda", 
                                          rownames(model_summary)), ]
names(hy_lambda_dframe) <- paste("pred_", 
                              names(hy_lambda_dframe),
                              sep = "")
require(stringr)
names(hy_lambda_dframe) <- str_sub(names(hy_lambda_dframe), 
                                      start = 1, 
                                    end = -2)
# Add hour:
hy_lambda_dframe$hour <- unique(data_clx$HORA)

#Predictions of the regions
sitio_lambda_dframe <- model_summary[grep("lambda_sitio", 
                                         rownames(model_summary)), ]
names(sitio_lambda_dframe) <- paste("pred_", 
                                   names(sitio_lambda_dframe),
                                   sep = "")
names(sitio_lambda_dframe) <- str_sub(names(sitio_lambda_dframe), 
                                     start = 1, 
                                     end = -2)
# Add a identifier for each region:
for (i in 1:length(unique(data_clx$SITIO))) {
  sitio_lambda_dframe[grep(paste(",", i, "]", sep = ""), 
                          rownames(sitio_lambda_dframe)),
                     "Sitio"] <- paste("Sitio ", i, sep = "")
}

# Add hour:
sitio_lambda_dframe$hour <- unique(data_clx$HORA)

#DATA FRAME FOR THE REGIONS
predict_boca <- sitio_lambda_dframe[sitio_lambda_dframe$Sitio == "Sitio 1", ]
predict_cnl <- sitio_lambda_dframe[sitio_lambda_dframe$Sitio == "Sitio 2", ]
predict_ens <- sitio_lambda_dframe[sitio_lambda_dframe$Sitio == "Sitio 3", ] 


# Data frame for the R-Squared:
rsqrd_dframe <- model_summary[grep("bay_r_sq", 
                                   rownames(model_summary)), ]


# Plot
require(ggplot2)
post_posteriors <- ggplot() +
  #Predictions:
  geom_ribbon(data = hy_lambda_dframe,
              aes(x = hour,
                  ymin = pred_25,
                  ymax = pred_75),
              fill = "black",
              color = "black",
              linetype = "dashed",
              alpha = 0.3) +
  geom_path(data = hy_lambda_dframe,
            aes(x = hour,
                y = pred_50),
            color = "black") +
  # Observed counts:
  geom_point(data = data_clx,
             aes(x = HORA,
                 y = Noftrns),
                 color = "black") + 
  # Predictions for lambda sitio 1 = BOCA:
  geom_path(data = predict_boca,
              aes(x = hour,
                  y = pred_50),
              color = "red") + 
  # Predictions for lambda sitio 2 = CANAL:
  geom_path(data = predict_cnl,
            aes(x = hour,
                y = pred_50),
            color = "blue") +
  # Predictions for lambda sitio 3 = ENSENADA:
  geom_path(data = predict_ens,
            aes(x = hour,
                y = pred_50),
            color = "yellow") +
  # R-squared:
  geom_text(aes(x = 3,
                y = 2000),
            color = "red",
            label = paste("BR2 = ",
                          round(rsqrd_dframe$`50%`,
                                digits = 2),
                          " (", round(rsqrd_dframe$`2.5%`,
                                      digits = 2), 
                          " - ", round(rsqrd_dframe$`97.5%`,
                                       digits = 2), ")",
                          sep = "")) +
  theme_bw() +
  labs(x = "Hora", y = "Counts of echolocation trains")

x11()
post_posteriors


#PLOT FOR SITIO 1 = BOCA
boca_df <- data_clx[data_clx$SITIO == 'BOCA', ]

count_trns_boca <- ggplot() +
  geom_point(data = boca_df,
             aes(x = HORA,
                 y = Noftrns),
                 color = "red") + 
  geom_path(data = predict_boca,
            aes(x = hour,
                y = pred_50),
            color = "darkred") 
  
  
x11()
count_trns_boca

#PLOT FOR SITIO 2 = CANAL
cnl_df <- data_clx[data_clx$SITIO == 'CANAL', ]

count_trns_cnl <- ggplot() +
  geom_point(data = cnl_df,
             aes(x = HORA,
                 y = Noftrns),
             color = "blue") + 
  geom_path(data = predict_cnl,
            aes(x = hour,
                y = pred_50),
            color = "darkblue") 


x11()
count_trns_cnl


#PLOT FOR SITIO 3 = ENSENADA
ens_df <- data_clx[data_clx$SITIO == 'ENSENADA', ]

count_trns_ens <- ggplot() +
  geom_point(data = ens_df,
             aes(x = HORA,
                 y = Noftrns),
             color = "green") + 
  geom_path(data = predict_ens,
            aes(x = hour,
                y = pred_50),
            color = "darkgreen") 


x11()
count_trns_ens