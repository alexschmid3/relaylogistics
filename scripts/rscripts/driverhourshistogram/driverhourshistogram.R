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

#---------------------------------------------------------------#

#Read in data
alldata <- read_csv('driverhours_combined.csv')

#---------------------------------------------------------------#

maxhours <- alldata %>%
  group_by(legendname) %>%
  summarize(max=max(hours))
maxhours[1,2] = 36.8
maxhours[2,2] = 35.9

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="driverequityhistogram.png",
    width=2400, height=1800)

alldata %>%
  ggplot(aes(x=hours, color=legendname, fill=legendname)) +
  geom_density(alpha=0.3, size=3) + 
  geom_vline(data = maxhours, aes(xintercept = max, color = legendname), size=5)+
  scale_colour_manual("", breaks = c("With driver equity", "No driver equity"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("With driver equity", "No driver equity"),
                    values = c("#DA853E", "#649BCB")) +
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  labs(x="Delivery time (hours)", y="Density") +
  scale_x_continuous(limits=c(0,40)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 70)) +
  theme(axis.text = element_text(size = 60)) +
  theme(axis.title = element_text(size = 75))

dev.off()



