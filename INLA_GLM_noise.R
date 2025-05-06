setwd("~/R/GLM")
# Clean:
rm(list = ls(all = TRUE)) #clear all
graphics.off() # close all
gc() # Clear memory

#load libraries 
library(dplyr)
library(reshape2)
library (ggplot2)
library(ggbreak)
library(lubridate)
library(viridis)
library(scales)
library(suncalc)
library(rcompanion)
library(psych)

#read database
data <- read.csv("base_glm_final.csv")
data <- data %>% 
  filter(!is.na(SEL_mid)) %>%
  filter(!SITE == "ENSENADA")
data$SITE <- as.factor(data$SITE)


# Scale covariates:
data$julian_sc <- (data$julian_day - 
                        mean(data$julian_day,
                             na.rm = T)) /  sd(data$julian_day, 
                                               na.rm = T)
data$julian_sc_copy <- data$julian_sc

data$SEL_mid_sc <- (data$SEL_mid -
                      mean(data$SEL_mid,
                           na.rm = TRUE))/ sd(data$SEL_mid,
                                              na.rm = TRUE)
data$SHIP_sc <- (data$SHIP -
                      mean(data$SHIP,
                           na.rm = TRUE))/ sd(data$SHIP,
                                              na.rm = TRUE)
data$tide_sc <- (data$altura_mm -
                   mean(data$altura_mm,
                        na.rm = TRUE))/ sd(data$altura_mm,
                                           na.rm = TRUE)



# Add type of row:
data$type <- "Observations"
# x11()
# hist(data$SEL_mid)
# plotNormalDensity(data$SEL_mid) 

# INLA MODELS:
require(INLA)
models_list <- list()

models_list[[1]] <- SEL_mid ~ julian_sc
models_list[[2]] <- SEL_mid ~ julian_sc + I(SHIP_sc)
models_list[[3]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + f(timeofday,
                                                          model = "seasonal",
                                                          season.length = 4)


#best model for short term
models_list[[4]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + f(HOUR,
                                                          model = "seasonal",
                                                          season.length = 24) + f(julian_sc_copy,
                                                                                  model = "rw1")
models_list[[5]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + f(HOUR,
                                                          model = "seasonal",
                                                          season.length = 24) + f(julian_sc_copy,
                                                                                  model = "ar1")
# models_list[[1]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + f(HOUR,
#                                                           model = "seasonal",
#                                                           season.length = 24) + f(julian_sc_copy,
#                                                                                   model = "ar",
#                                                                                   order = 3)
models_list[[6]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + I(tide_sc) +f(HOUR,
                                                                      model = "seasonal",
                                                                      season.length = 24) + f(julian_sc_copy,
                                                                                              model = "ar1")
models_list[[7]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + I(tide_sc) +f(HOUR,
                                                                      model = "seasonal",
                                                                      season.length = 24) + f(julian_sc_copy,
                                                                                              model = "rw1")
models_list[[8]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + I(tide_sc) +f(HOUR,
                                                                      model = "seasonal",
                                                                      season.length = 24) + f(julian_sc_copy,
                                                                                              model = "ar",
                                                                                              order = 3)
models_list[[9]]  <- SEL_mid ~ julian_sc + I(SHIP_sc) + I(tide_sc) +f(SITE,
                                                                      model = "iid")
#  


# Likelihoods:
families_list <- list()
names(inla.models()$likelihood)

# families_list[[1]] <- "gaussian"
families_list[[1]] <- "lognormal"



#### BUILD DATAFRAMES FOR PREDICTIONS:

# Hypermean predictions only:
n_predictions <- 100
hyp_predict_df <- as.data.frame(matrix(NA,
                                       nrow = n_predictions,
                                       ncol = ncol(data)))
names(hyp_predict_df) <- names(data)
hyp_predict_df$julian_sc <- seq(min(data$julian_sc, na.rm = T),
                              max(data$julian_sc, na.rm = T),
                              length = nrow(hyp_predict_df))
hyp_predict_df$julian_day <- seq(min(data$julian_day, na.rm = T),
                           max(data$julian_day, na.rm = T),
                           length = nrow(hyp_predict_df))
hyp_predict_df$HOUR <- seq(min(data$HOUR, na.rm = T),
                                 max(data$HOUR, na.rm = T),
                                 length = nrow(hyp_predict_df))
hyp_predict_df$type <- "Hyper_mean_predictions"

# Partial seasonal predictions:
seas_df <- as.data.frame(matrix(NA,
                                ncol = ncol(data),
                                nrow = 24))
names(seas_df) <- names(data)
seas_df$HOUR <- 0:23
seas_df$type <- "seasonal_predictions"

# Partial long-term fixed effects:
longt_df <- as.data.frame(matrix(NA,
                                 ncol = ncol(data),
                                 nrow = 100))
names(longt_df) <- names(data)
longt_df$julian_day <- seq(min(data$julian_day),
                   max(data$julian_day),
                   length = nrow(longt_df))
longt_df$type <- "long_term_predictions"

# Predictions for each site:
uniq_site <- unique(data$SITE)
site_df_list <- list()
for (i in 1:length(uniq_site)) {
  # i <- 1
  site_df_pred <- hyp_predict_df
  site_df_pred$SITE <- uniq_site[i]
  site_df_pred$type <- "site_predictions"
  site_df_list[[i]] <- site_df_pred

}

site_df_pred <- do.call("rbind",
                       site_df_list)

# Join data frames for INLA inference:
data <- rbind(data,
              seas_df,
              site_df_pred,
              hyp_predict_df,
              longt_df)

## RUN ALL MODELS:

# Prelocate:
model_objs <- list()
summ_results <- list()
list_pit <- list()
models_df <- data.frame()
count <- 0

require(svMisc) # for progress bar

# Loop by alternative models:
for (i in 1:length(models_list)) {
  # i <- 1
  # Loop on alternative likelihoods:
  for (j in 1:length(families_list)) {
    # j <- 1
    count <- count + 1
    
    gc()
    progress(count,
             max.value = length(models_list)*length(families_list)) # set progress bar
    
    
    # Run model
    model_objs[[count]] <- inla(data = data,
                                formula = models_list[[i]],
                                family = families_list[[j]],
                                control.predictor = list(compute = TRUE,
                                                         link = 1),
                                
                                quantiles = c(0.025, 0.125, 0.25,
                                              0.5,
                                              0.75, 0.875, 0.975),   # CI intervals for results
                                
                                control.compute = list(
                                  return.marginals.predictor = TRUE,  # WARNING!! Too slow (calculates marginals of each prediction)
                                  cpo = TRUE,   # Conditional predictive ordinate (CPO): # Similar to the Bayesian p-value
                                  # Cross-validation criterion for model assessment
                                  # Large values indicate a better fit of the model to the data,
                                  # small values indicate a bad fitting of the model to that observation and,
                                  # perhaps, that it is an outlier.
                                  # Since the sum of all CPO of all values is converted to negative,
                                  # smaller values point to a better model fit.
                                  dic = TRUE,   # Information criterion for model comparison
                                  mlik = TRUE,  # Marginal likelihood:
                                  # Probability of the observed
                                  # data under a given model
                                  # higher values are better.
                                  waic = TRUE,  # Information criterion for model comparison
                                  config = TRUE) # To be able to sample from the posteriors
    ) 
    
    # Summarized results:
    summ_results[[count]] <- summary(model_objs[[count]])
    
    # Save diagnostics into data.frame:
    models_df[count, "model_no"] <- count
    models_df[count, "waic"] <- summ_results[[count]]$waic$waic
    models_df[count, "dic"] <- summ_results[[count]]$dic$dic
    models_df[count, "cpo"] <- sum(summ_results[[count]]$cpo$cpo)
    models_df[count, "mlik_integ"] <- summ_results[[count]]$mlik[1,]
    models_df[count, "mlik_gaussian"] <- summ_results[[count]]$mlik[2,]
    models_df[count, "pit"] <- mean(summ_results[[count]]$cpo$pit, na.rm = TRUE)   # Predictive integral transform (PIT):
    # For each observation, it is the probability 
    # of a new value to be lower
    # than the actual observed value.
    
    # The mean is fine, but, it should be evaluated for each value of the observations. 
    # It is a measure of how far an observation is from the model's predictions.
    # Save PIT in list.
    list_pit[[count]] <- c(summ_results[[count]]$cpo$pit)
    
    models_df[count, "formula"] <- paste("No.", i, as.character(c(models_list[[i]])), sep = " ")
    models_df[count, "family"] <- as.character(families_list[[j]])
    
    # Close progress bar:
    if (count == length(models_list)*length(families_list)) message("Done!")
    gc()
  }
}

# Organize the data frame by lowest WAIC:
models_df <- models_df[order(models_df$waic),]
models_df$delta_waic <- cumsum(c(0, diff(models_df$waic)))


# Choose model to plot:
model_objs <- model_objs[[as.numeric(row.names(models_df)[1])]]

# Summary mode:
summ_results <- summ_results[[as.numeric(row.names(models_df)[1])]]

# Predictive integral transform (PIT):
pit_vector <- c(list_pit[[as.numeric(row.names(models_df)[1])]])

#extract the values of the marginals posterior
marg.fixed.intercept <- inla.tmarginal(fun = "exp",
                                       marginal = model_objs$marginals.fixed[[1]])
marg.fixed.jul <- inla.tmarginal(fun = "exp",
                                       marginal = model_objs$marginals.fixed[[2]])
marg.fixed.ship <- inla.tmarginal(fun = "exp",
                                 marginal = model_objs$marginals.fixed[[3]])
# marg.rndm.hour <- inla.tmarginal(fun = "exp",
#                                   marginal = model_objs[["marginals.random"]][["HOUR"]][["index.1"]])

#plot marginal for the fixed effects
x11()
marg.fixed.intercept_df <- as.data.frame(marg.fixed.intercept)
ggplot(data = marg.fixed.intercept_df,
       aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_path(color = "blue", size = 0.3)
x11()
marg.fixed.jul_df <- as.data.frame(marg.fixed.jul)
ggplot(data = marg.fixed.jul_df,
       aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_path(color = "blue", size = 0.3)
x11()
marg.fixed.ship_df <- as.data.frame(marg.fixed.ship)
ggplot(data = marg.fixed.ship_df,
       aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_path(color = "blue", size = 0.3)

#plot marginal for the random effects


# Add fitted values to original observations:
names(model_objs$summary.fitted.values) <- paste("fit_", 
                                                 names(model_objs$summary.fitted.values), 
                                                 sep = "")
data <- cbind(data, model_objs$summary.fitted.values)



# Split data frames:
obs_df <- data[data$type == "Observations", ]
hyp_mean_pred <- data[data$type == "Hyper_mean_predictions", ]
site_pred <- data[data$type == "site_predictions", ]
seas_pred <- data[data$type == "seasonal_predictions",]
longt_df<- data[data$type == "long_term_predictions",]

# Bayesian R-squared:
mean_r_squared <- var(exp(obs_df$fit_0.5quant))/(var(exp(obs_df$fit_0.5quant)) + 
                                              var(obs_df$SEL_mid - exp(obs_df$fit_0.5quant), na.rm = TRUE))

low_r_squared <- var(exp(obs_df$fit_0.025quant))/(var(exp(obs_df$fit_0.025quant)) + 
                                                  var(obs_df$SEL_mid - exp(obs_df$fit_0.025quant), na.rm = TRUE))

hi_r_squared <- var(exp(obs_df$fit_0.975quant))/(var(exp(obs_df$fit_0.975quant)) + 
                                               var(obs_df$SEL_mid - exp(obs_df$fit_0.975quant), na.rm = TRUE))


# For the R-squared label:
r_sq_text = as.character(paste("italic(R[B]^2==",
                               round(mean_r_squared, 
                                     digits = 2),")",
                               sep = ""))

data_plot <- obs_df %>%
  filter(SEL_mid > 175) %>%
  filter(SEL_mid < 190)
seas_pred$hour_sc <- 0:23

spline_sel <- as.data.frame(c(spline(seas_pred$hour_sc, seas_pred$fit_0.5quant),
                              spline(seas_pred$hour_sc, seas_pred$fit_0.025quant),
                              spline(seas_pred$hour_sc, seas_pred$fit_0.975quant)))
colnames(spline_sel)<- c("hora",
                         "y",
                         "horax",
                         "low",
                         "horay",
                         "high")





#Plot:
  require(ggplot2)
model_plot <- ggplot() +
  # # Hypermean predictions:
  geom_ribbon(data = spline_sel,
              aes(x = horax,
                  ymin = exp(low),
                  ymax = exp(high)),
              color = NA,
              fill = "gray",
              alpha = 0.5) +
  #Season predictions - hour
  geom_path(data = spline_sel,
            aes(x = hora,
                y = exp(y)),
            color = "black")+
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))+
  scale_y_continuous(name = expression(paste("SEL (1kHz - 20 kHz) ( dB re 1 " , mu, Pa^2, s, ")")),
                     limits = c(182,187))+
  scale_x_continuous(name = "Hora del día",
                     limits = c(0,23),
                     breaks = seq(0,23,2))
# x11()
# model_plot

# ggsave("SELmid_short_pattern.jpg", plot = model_plot, dpi = 120,
#        width = 1536, height = 869, units = "px")
#Plot:
box_data <- ggplot() +
  # # Hypermean predictions:
  geom_boxplot(data = data_plot,
              aes(x = factor(HOUR),
                  y = SEL_mid,
                  color = SITE),
              lwd = 0.5)+
  scale_color_manual(values = c("#FF0000",
                                "#330099"),
                     labs(fill = "Sitio"))+
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = expression(paste("SEL (1kHz - 20 kHz) ( dB re 1 " , mu, Pa^2, s, ")")),
                     limits = c(175,190))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))+
  scale_x_discrete(name = "Hora del día",
                   breaks = seq(0,23,2))
# x11()
# box_data
# ggsave("SELmid_short_pattern_box.jpg", plot = box_data, dpi = 120,
#        width = 1536, height = 869, units = "px")
library(cowplot)
x11()
x<-plot_grid(box_data, model_plot,
          # labels = "AUTO",
          # label_x = 0, label_y = 0,
          # hjust = -0.5, vjust = -0.5,
          align = "v",
          ncol = 1,
          rel_widths = c(2,1))
x11()
x
ggsave("SELmid_short_mix.jpg", plot = x, dpi = 120,
       width = 1536, height = 869, units = "px")
