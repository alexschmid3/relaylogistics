#install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(dplyr)

guides(
  fill = guide_legend(
    override.aes = list(color = NA, alpha = 0.8)
  ),
  color = "none"
)


IBMblue = "#648FFF"
IBMpurple = "#785EDC"
IBMpink = "#DC267F"
IBMorange = "#FE6100"
IBMyellow = "#FFB000"

#---------------------------------------------------------------#

#Read in data
alldata <- read_csv('driverhours_combined.csv')

#---------------------------------------------------------------#

maxhours <- alldata %>%
  group_by(legendname) %>%
  summarize(max=max(hours))
avghours <- alldata %>%
  group_by(legendname) %>%
  summarize(mean=mean(hours))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="driverequityhistogram.png",
    width=2400, height=1800)

alldata %>%
  ggplot(aes(x=hours, color=legendname, fill=legendname)) +
  geom_density(alpha=0.3, size=3) + 
  geom_vline(data = maxhours, aes(xintercept = max, color = legendname), size=5)+
  #geom_vline(data = avghours, aes(xintercept = mean, color = legendname), linetype="dotted", size=5)+
  scale_colour_manual("", breaks = c("With driver equity", "No driver equity"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("With driver equity", "No driver equity"),
                    values = c("#DA853E", "#649BCB")) +
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  labs(x="Delivery time (hours)", y="Density") +
  scale_x_continuous(limits=c(0,48)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 70)) +
  theme(axis.text = element_text(size = 60)) +
  theme(axis.title = element_text(size = 75))

dev.off()



