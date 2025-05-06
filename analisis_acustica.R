setwd("~/R/Eco_cuantitativa/Data")

# Clean:
rm(list = ls(all = TRUE)) #clear all
graphics.off() # close all
gc() # Clear memory


#read the database
data <- read.csv("base_proyecto.csv")


library(ggplot2)
#explore the response variable

#histogram to explore one of the possible response variable
#detection rate of whistles per hour, there is an excess of 0s in the data
hist_drw <- hist(data$DR_W)




#explore the covariates for outliers
# my_plots<-list()
# 
# for(i in 8:18){
#   p1<-eval(substitute(
#     ggplot(data,aes(x=data[,i],y=c(1:nrow(data)), color = factor(SITE)))+
#       geom_rect(aes(xmin=mean(data[,i],na.rm=T)-3*sd(data[,i],na.rm=T),
#                     xmax=mean(data[,i],na.rm=T)+3*sd(data[,i],na.rm=T),
#                     ymin=0,ymax=nrow(data)),
#                 fill="palegreen")+
#       geom_vline(xintercept=mean(data[,i],na.rm=T),size=1)+
#       geom_vline(xintercept=mean(data[,i],na.rm=T)-3*sd(data[,i],na.rm=T),color="black",size=1)+
#       geom_vline(xintercept=mean(data[,i],na.rm=T)+3*sd(data[,i],na.rm=T),color="black",size=1)+
#       geom_point(aes(x=data[,i],y=c(1:nrow(data))))+
#       ylab("Observations")+
#       xlab(names(data[i]))+
#       theme_bw()+
#       theme(panel.grid=element_blank(),
#             axis.text.y=element_blank(),
#             axis.ticks.y=element_blank()))
#     ,list(i = i))
#   print(i)
#   # print(p1)
#   my_plots[[i]] <- p1
# 
# }
# 
# x11()
# gridExtra::grid.arrange(my_plots[[8]],my_plots[[9]],
#                         my_plots[[10]],my_plots[[11]],
#                         my_plots[[12]],my_plots[[13]],
#                         my_plots[[14]],my_plots[[15]],
#                         my_plots[[16]],my_plots[[17]],
#                         ncol=5)


# #correlation of the variables
# source("panel_cor.R")
# x11()
# pairs(data[,c(9:15)], upper.panel = panel.cor)

#explore the binary data W_01 against the predictors
library(dplyr)
data_bin <- data %>%
  select("DATE", "HOUR","SITE","W_01","SHIP","SST",
                   "SEL_low","SEL_mid","SEL_hi","solar_angle","altura_mm", 
                   "moon_illumination")

data_bin$period <- as.factor(data$period)
data_bin$HOUR <- as.factor(data$HOUR)
data_bin$DATE <- as.factor(data$DATE)
data_bin$SITE <- as.factor(data$SITE)

# data.canal <- data.2 %>%
#   filter(SITE == "CANAL")
# 
# data.boca <- data.2 %>%
#   filter(SITE == "BOCA")

# my_plots.2<-list()
# 
# for(i in 3:10){
#   
#   p2<-eval(substitute(
#     ggplot(data = data.2,
#            aes(x=data.2[,i],
#                y=W_01))+
#       geom_point(data = data.canal,
#                  aes(x = data.canal[,i],
#                      y = W_01),
#                  color = "red", alpha=.5)+
#       geom_point(data = data.boca,
#                  aes(x = data.boca[,i],
#                      y = W_01),
#                  color = "gray60", alpha=.5)+
#       geom_smooth(data = data.canal,
#                   aes(x = data.canal[,i],
#                       y = W_01),
#                   color="yellow",fill="green")+
#       geom_smooth(data = data.boca,
#                   aes(x = data.boca[,i],
#                       y = W_01),
#                   color="black",fill="dodgerblue4")+
#       xlab(names(data.2[i]))+
#       theme_bw()+
#       theme(panel.grid=element_blank(),
#             panel.background=element_rect(size=2,color="black"),
#             axis.text=element_blank(),
#             axis.title=element_text(size=20))),
#     list(i=i))
#   print(i)
#   print(p2)
#   my_plots.2[[i]]<-p2
# }
# 
# x11()
# gridExtra::grid.arrange(my_plots.2[[3]],my_plots.2[[4]],
#                         my_plots.2[[5]],my_plots.2[[6]],
#                         my_plots.2[[7]],my_plots.2[[8]],
#                         my_plots.2[[9]],my_plots.2[[10]],
#                         ncol=4)


# my_plots.2<-list()
# 
# for(i in 2:9){
#   
#   p2<-eval(substitute(
#     ggplot(data = data_bin,
#            aes(x=data_bin[,i],
#                y=W_01))+
#       geom_point(color = "red", alpha=.5)+
#       geom_smooth(color="black",fill="dodgerblue4")+
#       xlab(names(data_bin[i]))+
#       theme_bw()+
#       theme(panel.grid=element_blank(),
#             panel.background=element_rect(size=2,color="black"),
#             axis.text=element_blank(),
#             axis.title=element_text(size=20))),
#     list(i=i))
#   print(i)
#   print(p2)
#   my_plots.2[[i]]<-p2
# }

# x11()
# gridExtra::grid.arrange(my_plots.2[[3]],my_plots.2[[4]],
#                         my_plots.2[[5]],my_plots.2[[6]],
#                         my_plots.2[[7]],my_plots.2[[8]],
#                         my_plots.2[[9]],my_plots.2[[2]],
#                         ncol=4)


#VIF
library(car)
model_vif<-lm(W_01~SHIP+SST+SEL_low+SEL_mid+SEL_hi+solar_angle+altura_mm+
                moon_illumination+period,
              data=data_bin)
x <- as.data.frame(vif(model_vif))

#GLM binary model
library(lme4)
#run the complete model
model_complete_bin <- glmer(W_01 ~ scale(SHIP) + scale(SST) + scale(SEL_low) + 
                              scale(SEL_mid)+ scale(solar_angle) + 
                              scale(altura_mm) + scale(moon_illumination) + SITE + 
                              (1 | HOUR) + (1 | DATE), 
                            family = binomial, 
                            data = data_bin)
library(jtools)
#check the summary of the data
summ(model_complete_bin)
acf(resid(model_complete_bin))
library(arm)
x <- predict(model_complete_bin)
y <- resid(model_complete_bin)
binnedplot(x,y)
# 
# x11()
# coefplot(model_complete_bin)

coef.df_complete<-data.frame(summary(model_complete_bin)$coefficients)
coef.df_complete$variable <- c("Intercept","# de embarcaciones","SST",
                               "SEL (0-1kHz)","SEL (1-20kHz)",           
                               "Ángulo del sol","Altura de la marea",
                               "Iluminación de la luna",
                               "SITIO (CANAL)") 
coef.df_complete$significancy <- ifelse(abs(coef.df_complete$Estimate)>1.96*coef.df_complete$Std..Error,
                               "significant",
                               "no-significant")
coef.df_complete$sign <- ifelse(coef.df_complete$Estimate>0,
                       "positive",
                       "negative")
coef.df_complete$sign<- paste(coef.df_complete$significancy,coef.df_complete$sign)
coef.df_complete <- subset(coef.df_complete, variable != "Intercept")

                                      
# x11()
effects_complete <-ggplot()+
  geom_point(data = coef.df_complete,
             aes(x = variable,
                 y = Estimate,
                 fill = sign),
             size = 5,
             shape = 21,
             stroke = 1)+
  geom_hline(yintercept=0,linetype=2)+
  geom_linerange(data = coef.df_complete,
                  aes(x = variable,
                      ymin = Estimate-1.96*Std..Error,
                      ymax = Estimate+1.96*Std..Error),
                 color = "black",
                 size = .7,
                 linewidth=1,
                 position=position_dodge(.3))+
  scale_fill_manual(values=c("gray60","gray60", "dodgerblue"))+
  coord_flip()+
  ylab("Tamaño de los efectos")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ggtitle("Efectos para el modelo completo")

  


#model for only night AND DAY
data_night <- data_bin %>%
  filter(period =="NIGHT")
data_day <-data_bin %>%
  filter(period =="DAY")

#run the model for night
model_bin_night <- glmer(W_01 ~ scale(SHIP) + scale(SST) + scale(SEL_low) + 
                         scale(SEL_mid)+ scale(solar_angle) + 
                         scale(altura_mm) + scale(moon_illumination) + SITE + 
                         (1|HOUR) + (1|DATE), 
                       family = binomial, 
                       data = data_night)

#check the summary of the data
summ(model_bin_night)

# x <- predict(model_bin_night)
# y <- resid(model_bin_night)
# binnedplot(x,y)


# #plot the coefficients
# x11()
# coefplot(model_bin_night)
# acf(resid(model_bin_night))


coef.df_night<-data.frame(summary(model_bin_night)$coefficients)
coef.df_night$variable <- c("Intercept","# de embarcaciones","SST",
                            "SEL (0-1kHz)","SEL (1-20kHz)",           
                            "Ángulo del sol","Altura de la marea",
                            "Iluminación de la luna",
                            "SITIO (CANAL)")
coef.df_night$significancy <- ifelse(abs(coef.df_night$Estimate)>1.96*coef.df_night$Std..Error,
                                        "significant",
                                        "no-significant")
coef.df_night$sign <- ifelse(coef.df_night$Estimate>0,
                                "positive",
                                "negative")
coef.df_night$sign<- paste(coef.df_night$significancy,coef.df_night$sign)
coef.df_night <- subset(coef.df_night, variable != "Intercept")


# x11()
effects_night<-ggplot()+
  geom_point(data = coef.df_night,
             aes(x = variable,
                 y = Estimate,
                 fill = sign),
             size = 5,
             shape = 21,
             stroke = 1)+
  geom_hline(yintercept=0,linetype=2)+
  geom_linerange(data = coef.df_night,
                 aes(x = variable,
                     ymin = Estimate-1.96*Std..Error,
                     ymax = Estimate+1.96*Std..Error),
                 color = "black",
                 size = .7,
                 linewidth=1,
                 position=position_dodge(.3))+
  scale_fill_manual(values=c("gray60","gray60", "dodgerblue"))+
  coord_flip()+
  ylab("Tamaño de los efectos")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ggtitle("Efectos para la noche")



#model for day data
model_bin_day<- glmer(W_01 ~ scale(SHIP) + scale(SST) + scale(SEL_low) + 
                        scale(SEL_mid)+ scale(solar_angle) + 
                        scale(altura_mm) + scale(moon_illumination) + SITE + 
                        (1|HOUR) + (1|DATE), 
                      family = binomial, 
                      data = data_day)

#check the summary of the data
summ(model_bin_day)

# x <- predict(model_bin_day)
# y <- resid(model_bin_day)
# binnedplot(x,y)

# x11()
# coefplot(model_bin_day)
# acf(resid(model_bin_day))

coef.df_day<-data.frame(summary(model_bin_day)$coefficients)
coef.df_day$variable <- c("Intercept","# de embarcaciones","SST",
                          "SEL (0-1kHz)","SEL (1-20kHz)",           
                          "Ángulo del sol","Altura de la marea",
                          "Iluminación de la luna",
                          "SITIO (CANAL)")
coef.df_day$significancy <- ifelse(abs(coef.df_day$Estimate)>1.96*coef.df_day$Std..Error,
                                     "significant",
                                     "no-significant")
coef.df_day$sign <- ifelse(coef.df_day$Estimate>0,
                             "positive",
                             "negative")
coef.df_day$sign<- paste(coef.df_day$significancy,coef.df_day$sign)
coef.df_day <- subset(coef.df_day, variable != "Intercept")


# x11()
effects_day <- ggplot()+
  geom_point(data = coef.df_day,
             aes(x = variable,
                 y = Estimate,
                 fill = sign),
             size = 5,
             shape = 21,
             stroke = 1)+
  geom_hline(yintercept=0,linetype=2)+
  geom_linerange(data = coef.df_day,
                 aes(x = variable,
                     ymin = Estimate-1.96*Std..Error,
                     ymax = Estimate+1.96*Std..Error),
                 color = "black",
                 size = .7,
                 linewidth=1,
                 position=position_dodge(.3))+
  scale_fill_manual(values=c("gray60","gray60", "dodgerblue"))+
  coord_flip()+
  ylab("Tamaño de los efectos")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  ggtitle("Efectos para el día")

###################################################
#predictions for day data

predictions_day <- data_day  
predicted_fit <- as.data.frame(predict(model_bin_day, newdata = predictions_day, 
                                          se = T,
                                          type = "response"))
predictions_day <- cbind(predictions_day,predicted_fit)
predictions_day$model <- "DÍA"

predictions_night <- data_night 
predicted_fit_night <- as.data.frame(predict(model_bin_night, newdata = predictions_night, 
                                       se = T,
                                       type = "response"))
predictions_night <- cbind(predictions_night,predicted_fit_night)
predictions_night$model <- "NOCHE"

predictions.df<-rbind(predictions_day, predictions_night)


#predictions
# x11()
ship <- ggplot(predictions.df, aes(x = SHIP, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
                                                 color = "black", fill = "dodgerblue") +
  labs(x = "Num. de embarcaciones", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))

moon <- ggplot(predictions.df, aes(x = moon_illumination, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = "black", fill = "dodgerblue") +
  labs(x = "Iluminación de la luna", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))

sel_low <- ggplot(predictions.df, aes(x = SEL_low, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = "black", fill = "dodgerblue") +
  labs(x = "SEL (0-1kHz)", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))

sel_mid <- ggplot(predictions.df, aes(x = SEL_mid, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = "black", fill = "dodgerblue") +
  labs(x = "SEL (1-20kHz)", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))

solar_angle <- ggplot(predictions.df, aes(x = solar_angle, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = "black", fill = "dodgerblue") +
  labs(x = "Ángulo del sol", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))

sst <- ggplot(predictions.df, aes(x = SST, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = "black", fill = "dodgerblue") +
  labs(x = "SST", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))
tide <- ggplot(predictions.df, aes(x = altura_mm, y = fit)) +
  # geom_point(color = "gray60", alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              color = "black", fill = "dodgerblue") +
  labs(x = "Altura de la marea", y = "Probabilidad Predecida", color = "Observed W_01") +
  theme_classic()+
  facet_wrap(~model)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.line=element_line(linewidth=1))

# x11()
gridExtra::grid.arrange(ship, moon, sel_low, 
                        sel_mid,solar_angle, 
                        sst,tide,
                        ncol=3)

gridExtra::grid.arrange(effects_complete, effects_day, effects_night,
                        ncol = 2)


data_hour <- data %>% 
  group_by(HOUR, SITE) %>%
  summarise(W_01 = sum(W_01))

data_covariates <- data %>% 
  summarise(W_01 = sum(W_01))

ggplot(data = data_hour, aes(x = HOUR, y = W_01, fill = SITE))+
  geom_col()+
  scale_color_manual("red", "dodgerblue")+
  facet_wrap(~SITE)+
  ylab("Suma de horas con presencia de delfines")+
  xlab("Hora")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))





























