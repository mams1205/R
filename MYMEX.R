setwd("~/R/PSAW_FOSSA_tutorial")

# Clean:
rm(list = ls(all = TRUE)) #clear all
graphics.off() # close all
gc() # Clear memory

library(PAMpal)

# avoid popups by specifying all arguments
#ADD THE DB
db <- list.files('./BW43_daytime_review/',
                 pattern='sqlite', 
                 full.names=TRUE)
db<- db[1:3]

#ADD THE BINARIES
bin <- './BW43_daytime_review/MYMEX_Leg1_RT_BinaryStorage'
#
#ADD DB AND BINARIES TO PPS
pps <- PAMpalSettings(db=db,
                      bin = bin,
                      sr_hz='auto', filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025,
                      settings = 'XMLSettings.xml')



#
# #read excel with the events
# library(readxl)
# events <- read_excel('./BW43_daytime_review/MYMEX_leg1_BW43_summary')


# give the correct colnames
# library(dplyr)
# events <- events %>%
#   rename(start = "CT start date/time",
#          end = "CT end date/time",
#          id = "Event")
# events$db <- if_else(events$id > 8,
#                      './BW43_daytime_review/MYMEX_Leg1_RT_Tracking_20240609.sqlite3',
#                      './BW43_daytime_review/MYMEX_Leg1_BW43_Review.sqlite3')


#save the data file
eventdata <- processPgDetections(pps, mode='db', id='MYMEX')

saveRDS(eventdata, file='./MYMEX_eventdata.rds')
# eventdata <- readRDS('MYMEX_eventdata.rds')
#select only the events of BW43 
#create an empty vector for the length of the events
comments <- vector("numeric", 185)
for (i in 1:185){
  comments[i]<- eventdata@events[[i]]@ancillary[["eventComment"]]
}

#select the events that you want to drop
indices_to_drop <- which(comments == "Manual Click Train Detection"|is.na(comments))

#empty vector with all the events
indices_events <- seq(1:185)
#drop the events of the indices
indices_events <- indices_events[-indices_to_drop]
#calculate average spectrum and concatenated spectrogram

avspec <- calculateAverageSpectra(eventdata, 
                                  evNum=indices_events,plot = TRUE, 
                                  ylim= c(-25,0), filterfrom_khz = 15, showBreaks = F)
avspec_plot<- recordPlot(calculateAverageSpectra(eventdata, 
                                   evNum=indices_events,plot = TRUE, 
                                   ylim= c(-10,0), filterfrom_khz = 15, showBreaks = F))


avspec_df <- as.data.frame(avspec$freq)
avspec_df$spec <- avspec$avgSpec
colnames(avspec_df) <- c("freq", "spec")

avspec_df <- avspec_df %>%
  # filter(spec >= -10) %>%
  filter(freq > 15000)

peaks <- PAMmisc::peakTrough(cbind(avspec_df$freq/1e3, avspec_df$spec),
                                        plot=TRUE,
                                        smooth = 3,
                                        freqBounds=c(0,40))
peaks <- unlist(peaks[1,1:4])
peaks <- unname(peaks*1000)

peaks_plot <- avspec_df %>%
  filter(freq %in% peaks) %>%
  arrange(desc(spec))
peaks_plot$peaks <- c("Highest Peak", "Second Peak", "Third Peak", "Trough/Notch")



library(ggplot2)
colnames(avspec_df) <- c("freq", "spec")
peak_plot <- ggplot()+
  geom_path(data = avspec_df,
            aes(x = freq/1e3,
                y = spec)) +
  geom_point(data = peaks_plot,
             aes(x = freq/1e3,
                 y = spec,
                 color = peaks),
             size = 3.5)+
  scale_color_manual(values = c("darkblue", "dodgerblue", "skyblue", "red"))+
  ylab("Relative Amplitude (dB)")+
  xlab("Frequency (kHz)")+
  ggtitle("Power Spectrum Average")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
  
x11()
peak_plot


#create new event data with just the events you want
event_data_new<- eventdata[-indices_to_drop]

#calculate the ICI
ICI <-calculateICI(event_data_new)
ICI <-getICI(ICI, type = c("data"))

#click data info
clicks <- getClickData(event_data_new)

#select only the icis higher than 0
ICI <- ICI %>%
  filter(ici > 0)

#select only the ones of Bw43
clicks_f <- clicks 
  # filter(eventLabel == "Bw43")

#select all ici lower than 1
ICI_all <- ICI %>%
  filter(eventId %in% clicks_f$eventId) %>%
  filter(ici < 1)
#histogram plot for ici
ICI_plot<- ggplot(data = ICI_all,
                  aes(x = ici))+
  geom_histogram(color = "black", fill = "gray60")+
  theme_bw()+
  xlab("ICI (ms)")+
  ylab("")+
  ggtitle("ICI histogram")+
  theme(plot.title = element_text(hjust = 0.5))
  
ICI_plot

#plot waveform
plot_waveform <- recordPlot(plotWaveform(event_data_new,
                                         UID = 1707002208,
                                         # UID  =  clicks_f$UID[which(clicks_f$noiseLevel == min(clicks_f$noiseLevel))],
                                         length= 150))
wave_data <- plot_waveform[1]
x<- wave_data[[1]][[4]][[2]][[2]][["x"]]
y <- wave_data[[1]][[4]][[2]][[2]][["y"]]

wave_data <- data.frame(time = x, amplitude = y)
wave_data <- as.data.frame(c(spline(wave_data$time, wave_data$amplitude)))


wave_plot <- ggplot()+
  geom_path(data = wave_data,
            aes(x = x,
                y = y),
            color = "dodgerblue")+
  xlab("Time (s)")+
  ylab("Amplitude")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),       
        axis.text.y = element_blank())
wave_plot


plot_wigner <- recordPlot(plotWigner(event_data_new,
                                     UID = 1707002208,
                                     # UID  =  clicks_f$UID[which(clicks_f$noiseLevel == min(clicks_f$noiseLevel))],
                                     length = 75))
wigner_plot <- ggplot() +
  xlim(2.8, 3.4) +
  ylim(0, 100) +
  labs( x = "Time (ms)", y = "Frequency (kHz)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )
x11()
wigner_plot

x11()
cowplot::plot_grid(wave_plot,wigner_plot,
                   peak_plot, ICI_plot, 
                   ncol = 2,
                   labels = c("A", "B", "C", "D"))








