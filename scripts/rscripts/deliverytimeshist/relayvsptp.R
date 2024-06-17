#install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#

#Read in data
ptpdata <- read_csv('ex4_exp64_ptp_rundate2024-06-13_deliverytimes.csv')
relaydata <- read_csv('ex4_exp69_relay_rundate2024-06-13_deliverytimes.csv')

#---------------------------------------------------------------#

#Filter to orders completed in both runs


#---------------------------------------------------------------#

relaydata1 <- relaydata %>% select(id, order, actual)
vec1 <- rep("Actual delivery time", nrow(relaydata1))
relaydata1 <- relaydata1 %>% mutate(datalabel = vec1)
relaydata1 <- rename(relaydata1, deliverytime = actual)

relaydata2 <- relaydata %>% select(id, order, best)
vec2 <- rep("Best delivery time", nrow(relaydata2))
relaydata2 <- relaydata2 %>% mutate(datalabel = vec2)
relaydata2 <- rename(relaydata2, deliverytime = best)

relaydata <- rbind(relaydata1, relaydata2) 

#---------------------------------------------------------------#

relaymeandelivtime <- relaydata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(deliverytime))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="relaydelivtimes_raw.png",
    width=2400, height=1800)

relaydata %>%
  ggplot(aes(x=deliverytime, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3, size=3)+ 
  geom_vline(data = relaymeandelivtime, aes(xintercept = mean, 
                                     color = datalabel), size=5)+
  scale_colour_manual("", breaks = c("Actual delivery time", "Best delivery time"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Actual delivery time", "Best delivery time"),
                    values = c("#DA853E", "#649BCB")) +
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  labs(x="Delivery time (hours)", y="Density") +
  scale_x_continuous(limits=c(0,150)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 70)) +
  theme(axis.text = element_text(size = 60)) +
  theme(axis.title = element_text(size = 75))

dev.off()


#---------------------------------------------------------------#

ptpdata <- read_csv('ex4_exp64_ptp_rundate2024-06-13_deliverytimes.csv')


ptpdata1 <- ptpdata %>% select(id, order, actual)
vec1 <- rep("Actual delivery time", nrow(ptpdata1))
ptpdata1 <- ptpdata1 %>% mutate(datalabel = vec1)
ptpdata1 <- rename(ptpdata1, deliverytime = actual)

ptpdata2 <- ptpdata %>% select(id, order, best)
vec2 <- rep("Best delivery time", nrow(ptpdata2))
ptpdata2 <- ptpdata2 %>% mutate(datalabel = vec2)
ptpdata2 <- rename(ptpdata2, deliverytime = best)

ptpdata <- rbind(ptpdata1, ptpdata2) 

#---------------------------------------------------------------#

ptpmeandelivtime <- ptpdata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(deliverytime))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="ptpdelivtimes_raw.png",
    width=2400, height=1800)

ptpdata %>%
  ggplot(aes(x=deliverytime, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3, size=3)+ 
  geom_vline(data = ptpmeandelivtime, aes(xintercept = mean, 
                                            color = datalabel), size=5)+
  scale_colour_manual("", breaks = c("Actual delivery time", "Best delivery time"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Actual delivery time", "Best delivery time"),
                    values = c("#DA853E", "#649BCB")) +
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  labs(x="Delivery time (hours)", y="Density") +
  scale_x_continuous(limits=c(0,150)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 70)) +
  theme(axis.text = element_text(size = 60)) +
  theme(axis.title = element_text(size = 70))

dev.off()


