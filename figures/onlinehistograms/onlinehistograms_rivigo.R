install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
#library(leaflet)
#library(lubridate)
library(dplyr)
#library(network)
#library(outliers)
#library(geosphere)
#library(fields)
#library(hash)
#library(naniar)

#---------------------------------------------------------------#
#---------------------------------------------------------------#

#Read in detour and delay distribution data
abcgdata <- read_csv('orderhist_abcg.csv')
otddata <- read_csv('orderhist_otd.csv')
rrdata <- read_csv('orderhist_rr.csv')

#---------------------------------------------------------------#
#---------------------------------------------------------------#

joineddata <- full_join(abcgdata, otddata, by = "orderid") %>%
  select(orderid, orderdetour.x, orderdetour.y)

joineddata <- rename(joineddata, detour_abcg = orderdetour.x)
joineddata <- rename(joineddata, detour_otd = orderdetour.y)

joineddata <- full_join(joineddata, rrdata, by = "orderid") %>%
  select(orderid, detour_abcg, detour_otd, orderdetour)

joineddata <- rename(joineddata, detour_rr = orderdetour)

vec_mvg <- rep("RPDP/MAG", 831)
vec_otd <- rep("Orders-then-drivers", 831)
vec_rr <- rep("Orders-then-drivers (historical) ", 831)

mvghist <- joineddata %>% select(orderid, detour_abcg)
otdhist <- joineddata %>% select(orderid, detour_otd)
rrhist <- joineddata %>% select(orderid, detour_rr)

mvghist <- rename(mvghist, orderdetour = detour_abcg)
otdhist <- rename(otdhist, orderdetour = detour_otd)
rrhist <- rename(rrhist, orderdetour = detour_rr)

mvghist <- mvghist %>%
  mutate(datalabel = vec_mvg)
otdhist <- otdhist %>%
  mutate(datalabel = vec_otd)
rrhist <- rrhist %>%
  mutate(datalabel = vec_rr)

alldata = rbind(mvghist, rrhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("RPDP/MAG", "Orders-then-drivers", "Orders-then-drivers (historical) "), 
                            labels=c("RPDP/MAG", "Orders-then-drivers", "Orders-then-drivers (historical) "))

alldata <- subset(alldata, !is.na(orderdetour))

mean_detour <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(orderdetour))

#---------------------------------------------------------------#

#Detour
update_geom_defaults("text", list(size = 48))

dev.new(width=8, height=4)

png(file="detourdensity_leg.png", width=900, height=1200)
   # width=450, height=600)

alldata %>%
  ggplot(aes(x=orderdetour, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3,size=1)+ 
  geom_vline(data = mean_detour, aes(xintercept = mean, 
                                     color = datalabel), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("RPDP/MAG", "Orders-then-drivers (historical) "),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("RPDP/MAG", "Orders-then-drivers (historical) "),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Detour (as % of shortest path)", y="Density")+
  scale_x_continuous(limits=c(0,0.15), labels = scales::percent_format(accuracy = 1)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 30)) +
  theme(axis.title = element_text(size = 30))

dev.off()

#---------------------------------------------------------------#

joineddata <- full_join(abcgdata, otddata, by = "orderid") %>%
  select(orderid, orderdelay.x, orderdelay.y)

joineddata <- rename(joineddata, delay_abcg = orderdelay.x)
joineddata <- rename(joineddata, delay_otd = orderdelay.y)

joineddata <- full_join(joineddata, rrdata, by = "orderid") %>%
  select(orderid, delay_abcg, delay_otd, orderdelay)

joineddata <- rename(joineddata, delay_rr = orderdelay)

vec_mvg <- rep("Optimized solution", 831)
vec_otd <- rep("Orders-then-drivers", 831)
vec_rr <- rep("Orders-then-drivers (historical) ", 831)

mvghist <- joineddata %>% select(orderid, delay_abcg)
otdhist <- joineddata %>% select(orderid, delay_otd)
rrhist <- joineddata %>% select(orderid, delay_rr)

mvghist <- rename(mvghist, orderdelay = delay_abcg)
otdhist <- rename(otdhist, orderdelay = delay_otd)
rrhist <- rename(rrhist, orderdelay = delay_rr)

mvghist <- mvghist %>%
  mutate(datalabel = vec_mvg)
otdhist <- otdhist %>%
  mutate(datalabel = vec_otd)
rrhist <- rrhist %>%
  mutate(datalabel = vec_rr)

alldata = rbind(mvghist, rrhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("Optimized solution", "Orders-then-drivers", "Orders-then-drivers (historical) "), 
                            labels=c("Optimized solution", "Orders-then-drivers", "Orders-then-drivers (historical) "))

alldata <- subset(alldata, !is.na(orderdelay))

mean_delay <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(orderdelay))

#---------------------------------------------------------------#

#Delay
update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="delaydensity.png",
    width=450, height=600)

alldata %>%
  ggplot(aes(x=orderdelay, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3, size=1)+ 
  geom_vline(data = mean_delay, aes(xintercept = mean, 
                                    color = datalabel), size=1.5)+
  scale_colour_manual("", breaks = c("Optimized solution", "Orders-then-drivers (historical) "),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Orders-then-drivers (historical) "),
                    values = c("#DA853E", "#649BCB")) +
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  labs(x="Delay per completed order\n(% of shortest path)", y="Density")+
  scale_x_continuous(limits=c(0,2), labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position="none") +
  #theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20))

dev.off()

#---------------------------------------------------------------#
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#Read in driver empty distribution
abcgdata <- read_csv('driverempty_abcg.csv')
rrdata <- read_csv('driverempty_rr.csv')

#---------------------------------------------------------------#

vec_mvg <- rep("Optimized solution", 2998)
vec_rr <- rep("Orders-then-drivers (historical) solution", 3316)

mvghist <- abcgdata %>% select(driverid, emptymiles, percentempty)
rrhist <- rrdata %>% select(driverid, emptymiles, percentempty)

mvghist <- mvghist %>%
  mutate(datalabel = vec_mvg)
rrhist <- rrhist %>%
  mutate(datalabel = vec_rr)

alldata = rbind(mvghist, rrhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("Optimized solution", "Orders-then-drivers (historical) solution"), 
                            labels=c("Optimized solution", "Orders-then-drivers (historical) solution"))

alldata <- subset(alldata, !is.na(percentempty))
alldata <- subset(alldata, !is.na(emptymiles))

mean_pct <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(percentempty))

mean_empty <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(emptymiles))

#mean_pct["mean"][mean_pct["datalabel"] == "Orders-then-drivers (historical) solution"] <- 0.24
mean_pct["mean"][mean_pct["datalabel"] == "Optimized solution"] <- 0.13

#---------------------------------------------------------------#

#Driver empty percent
update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="driveremptypct.png",
    width=450, height=600)

alldata %>%
  ggplot(aes(x=percentempty, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3,size=1)+ 
  geom_vline(data = mean_pct, aes(xintercept = mean, 
                                  color = datalabel), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("Optimized solution", "Orders-then-drivers (historical) solution"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Orders-then-drivers (historical) solution"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Empty miles per driver\n(% of miles)", y="Density")+
  scale_x_continuous(limits=c(0,0.5), labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position="none") +
  #theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20))

dev.off()

#---------------------------------------------------------------#

#Read in remaining miles data
abcgdata <- read_csv('orderoutcomes_abcg.csv')
rrdata <- read_csv('orderoutcomes_rr.csv')

#---------------------------------------------------------------#

vec_mvg <- rep("Optimized solution", 831)
vec_rr <- rep("Actual routes solution", 831)

mvghist <- abcgdata %>% select(orderid, remainingmilespct)
rrhist <- rrdata %>% select(orderid, remainingmilespct)

mvghist <- rename(mvghist, remainingmiles = remainingmilespct)
rrhist <- rename(rrhist, remainingmiles = remainingmilespct)

mvghist <- mvghist %>%
  mutate(datalabel = vec_mvg)
rrhist <- rrhist %>%
  mutate(datalabel = vec_rr)

alldata = rbind(mvghist, rrhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("Optimized solution", "Actual routes solution"), 
                            labels=c("Optimized solution", "Actual routes solution"))

alldata <- subset(alldata, !is.na(remainingmiles))

mean_rem <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(remainingmiles))

#---------------------------------------------------------------#

#Order remaining miles
update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="remainingmilespct.png",
    width=450, height=600)


alldata %>%
  ggplot(aes(x=remainingmiles, color=datalabel, fill=datalabel)) +
  #geom_density(alpha=0.3,size=1)+ 
  #geom_histogram(binwidth = 0.05, position = 'identity', alpha=0.3)+
  geom_density(aes(y=..count..), alpha=0.3, size=1)+
  geom_vline(data = mean_rem, aes(xintercept = mean, color = datalabel), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 250)+
  scale_colour_manual("", 
                      breaks = c("Optimized solution", "Actual routes solution"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Actual routes solution"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Remaining miles per non-completed order \n(% of trip remaining)", y="Density")+
  scale_x_continuous(limits=c(0,1), labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position="none") +
  #theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20))

dev.off()

