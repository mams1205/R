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
data$SITE <- as.factor(data$SITE)
data$DATE <-as.Date(data$DATE)

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
data$tide_sc <- (data$altura_mm -
                   mean(data$altura_mm,
                        na.rm = TRUE))/ sd(data$altura_mm,
                                           na.rm = TRUE)
data$ship_sc <- (data$SHIP -
                   mean(data$SHIP,
                        na.rm = TRUE))/ sd(data$SHIP,
                                           na.rm = TRUE)
data$Dolclx_sc <- abs((data$Dolclx - 
                     mean(data$Dolclx,
                          na.rm = T)) /  sd(data$Dolclx, 
                                            na.rm = T))
 
# Add type of row:
data$type <- "Observations"
data <- data %>% 
  filter(!SITE == "ENSENADA")
# INLA MODELS:
require(INLA)
models_list <- list()

# models_list[[1]] <- Dolclx_01 ~ julian_sc
# models_list[[2]] <- Dolclx_01 ~ julian_sc + I(ship_sc) + f(SITE, model = "iid")
# models_list[[3]] <- Dolclx_01 ~ julian_sc + I(ship_sc)
# models_list[[4]] <- Dolclx_01 ~ julian_sc + f(HOUR,
#                                              model = "seasonal",
#                                              season.length = 24) + f(julian_sc_copy,
#                                                                      model = "rw1")
models_list[[1]] <- Dolclx_01 ~ julian_sc + I(ship_sc) + f(HOUR,
                                                        model = "seasonal",
                                                        season.length = 24) + f(julian_sc_copy,
                                                                                model = "ar",
                                                                                order = 3)
# models_list[[6]] <- Dolclx_01 ~ julian_sc + I(SEL_mid_sc) + f(HOUR,
#                                                         model = "seasonal",
#                                                         season.length = 24) + f(julian_sc_copy,
#                                                                                 model = "rw1")
# 
# models_list[[7]] <- Dolclx_01 ~ julian_sc + f(timeofday,
#                                               model = "seasonal",
#                                               season.length = 4) + f(julian_sc_copy,
#                                                                       model = "rw1")
# 
# models_list[[8]] <- Dolclx_01 ~ julian_sc + I(ship_sc) + f(timeofday,
#                                                         model = "seasonal",
#                                                         season.length = 4) + f(julian_sc_copy,
#                                                                                 model = "rw1")
# 
# models_list[[9]] <- Dolclx_01 ~ julian_sc + I(tide_sc) + I(ship_sc) + I(SEL_mid_sc) +  f(HOUR,
#                                                                                        model = "seasonal",
#                                                                                        season.length = 24) + f(julian_sc_copy,
#                                                                                                              model = "rw1")
# models_list[[10]] <- Dolclx_01 ~ julian_sc + I(tide_sc) + f(HOUR,
#                                                            model = "seasonal",
#                                                            season.length = 24) + f(julian_sc_copy,
#                                                                                       model = "rw1")
# models_list[[11]] <- Dolclx_01 ~ julian_sc + I(tide_sc) + I(ship_sc) + I(SEL_mid_sc) +  f(SITE,
#                                                                                          model = "iid")


# Likelihoods:
families_list <- list()
names(inla.models()$likelihood)

# families_list[[1]] <- "binomial"
families_list[[1]] <- "betabinomial"



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
hyp_predict_df$HOUR <- rep(seq(0, 23),
                     length = nrow(hyp_predict_df))

hyp_predict_df$type <- "Hyper_mean_predictions"

# Partial seasonal predictions:
seas_df <- as.data.frame(matrix(NA,
                                ncol = ncol(data),
                                nrow = 24))
names(seas_df) <- names(data)
seas_df$HOUR <- 0:23
seas_df$type <- "seasonal_predictions"

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
#predictions for number of ship
predict_ship <- as.data.frame(matrix(NA,
                                       nrow = 100,
                                       ncol = ncol(data)))
names(predict_ship) <- names(data)

predict_ship$ship_sc <- seq(min(data$ship_sc, na.rm = T),
                                 max(data$ship_sc, na.rm = T),
                                 length = nrow(predict_ship))

predict_ship$SHIP <- seq(min(data$SHIP, na.rm = T),
                            max(data$SHIP, na.rm = T),
                            length = nrow(predict_ship))

predict_ship$type <- "ship_predictions"


# Join data frames for INLA inference:
data <- rbind(data,
              hyp_predict_df,
              site_df_pred,
              seas_df,
              predict_ship)

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
                                control.link (model = "logit"),
                                
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

marg.fixed.intercept <- inla.tmarginal(fun = "exp",
                                       marginal = model_objs$marginals.fixed[[1]])

#plot marginal for the intercept
x11()
marg.fixed.intercept_df <- as.data.frame(marg.fixed.intercept)
ggplot(data = marg.fixed.intercept_df,
       aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_path(color = "blue", size = 0.3)


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
pred_df<- data[data$type == "Predictions",]
pred_ship<- data[data$type == "ship_predictions",]


# Bayesian R-squared:
mean_r_squared <- var(obs_df$fit_0.5quant)/(var(obs_df$fit_0.5quant) + 
                                              var(obs_df$Dolclx_01 - obs_df$fit_0.5quant, na.rm = TRUE))
low_r_squared <- var(obs_df$fit_0.025quant)/(var(obs_df$fit_0.025quant) + 
                                             var(obs_df$Dolclx_01 - obs_df$fit_0.025quant, na.rm = TRUE))
hi_r_squared <- var(obs_df$fit_0.975quant)/(var(obs_df$fit_0.975quant) + 
                                               var(obs_df$Dolclx_01 - obs_df$fit_0.975quant, na.rm = TRUE))

# For the R-squared label:
r_sq_text = as.character(paste("italic(R[B]^2==",
                               round(mean_r_squared, 
                                     digits = 2),")",
                               sep = ""))

prob_data <- data %>%
  group_by(HOUR) %>%
  summarise(sum(Dolclx_01, na.rm = TRUE)/length(Dolclx_01))
colnames(prob_data) <- c("HOUR", "prob")

sum_data <- data %>%
  group_by(HOUR) %>%
  summarise(sum(Dolclx_01, na.rm = TRUE))
colnames(sum_data) <- c("HOUR", "suma")

prob_data_ship <- data %>%
  group_by(HOUR) %>%
  summarise(sum(SHIP_01, na.rm = TRUE)/length(SHIP_01))
colnames(prob_data_ship) <- c("HOUR", "prob")

spline_clx <- as.data.frame(c(spline(seas_pred$HOUR, seas_pred$fit_0.5quant),
                              spline(seas_pred$HOUR, seas_pred$fit_0.025quant),
                              spline(seas_pred$HOUR, seas_pred$fit_0.975quant)))
colnames(spline_clx)<- c("hora",
                         "y",
                         "horax",
                         "low",
                         "horay",
                         "high")

#Plot:
data$rug <- c(0)
# data$DR_clx <- data$DR_clx/1000
require(ggplot2)
model_plot <- ggplot() +
# #Hypermean predictions:
# geom_ribbon(data = hyp_mean_pred,
#             aes(x = julian_day,
#                 ymin = exp(fit_0.025quant),
#                 ymax = exp(fit_0.975quant)),
#             color = NA,
#             fill = "gray",
#             alpha = 0.5) +
# # Population predictions:
# geom_path(data = site_pred,
#           aes(x = julian_day,
#               y = fit_0.5quant,
#               color = SITE)) +
#Observations:
# geom_point(data = prob_data,
#             aes(x = HOUR,
#                 y = prob),
#             color = "red",
#            size = 3) +
# geom_path(data = prob_data_ship,
#             aes(x = HOUR,
#                 y = prob),
#             color = "dodgerblue",
#             linewidth = 2)+
#Season predictions - hour
  geom_ribbon(data = pred_ship,
            aes(x = SHIP,
                ymin = fit_0.025quant,
                ymax = fit_0.975quant),
            color = "gray",
            alpha = 0.2)+
  geom_path(data = pred_ship,
          aes(y = fit_0.5quant,
              x = SHIP),
            color = "black",
            linewidth = 1)+
  # scale_x_continuous(breaks = seq(0,20,2))+
  ylab("Probabilidad de detección de clics") +
  geom_rug(data = data,
           aes(x = Dolclx_sc))+
  # scale_y_continuous(sec.axis = sec_axis(trans = ~.+1000, name = "clx"))+
  xlab("Número de dembarcaciones por hora")+
  theme_bw()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

x11()
model_plot

# ggsave("dolclx_short_pattern_final.jpg", plot = model_plot, dpi = 120,
#        width = 1536, height = 869, units = "px")

# data <- data %>%
#   filter(Dolclx > 0)
# data <- data %>%
#   filter(Dolclx < 5000)

box_plot <- ggplot()+
  geom_boxplot(data = data,
               aes(x = factor(HOUR),
                   y = DR_clx),
               color = "black",
               fill = NA,
               lwd = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))+
  scale_x_discrete(breaks = seq(0,23,2))+
  xlab("Hora del día")+
  ylab("Tasa de detección de clics/hora")
x11()
box_plot
# ggsave("dolclx_box_final.jpg", plot = box_plot, dpi = 120,
#        width = 1536, height = 869, units = "px")
library(cowplot)
x11()
x <- plot_grid(box_plot, model_plot,
          labels = "AUTO",
          align = "hv")

ggsave("dolclx_boxymodel.jpg", plot = x, dpi = 120,
       width = 1536, height = 869, units = "px")

# Sample data
set.seed(123)
df <- data.frame(
  x = 1:10,
  y1 = rnorm(10, mean = 10, sd = 2),
  y2 = rnorm(10, mean = 50, sd = 10)
)

# Plotting
plot <- ggplot(df, aes(x)) +
  geom_line(aes(y = y1, color = "Variable 1")) +
  geom_line(aes(y = y2 * 0.1, color = "Variable 2")) +
  scale_y_continuous(
    sec.axis = sec_axis(~.*10, name = "Variable 2")
  ) +
  labs(
    y = "Variable 1",
    color = "Variables"
  ) +
  theme_minimal()
x11()
plot

x <- data%>%
  group_by(SHIP) %>%
  summarise(mean(DR_clx, na.rm = TRUE))

