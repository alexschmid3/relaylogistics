
#install.packages(c("tidyverse", "shiny")) #, "leaflet", "DT", "hash"))
library(shiny)
library(tidyverse)
library(dplyr)
library(formattable)

#---------------------------------------------------------------#

#Read in data
data <- read_csv('orderdata.csv')
data <- data %>% filter(select == 1)

ggplot(data, aes(x=delivery_time))

#---------------------------------------------------------------#

mean1 <- data %>% summarise(mean_miles=mean(miles_travelled))
mean2 <- data %>% summarise(mean_miles=mean(delivery_time))

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 14))

dev.new(width=8, height=4)

png(file="../ordermiles.png",
    width=700, height=650)

#Make the histogram
ggplot(data, aes(x=miles_travelled)) +
  geom_histogram(binwidth=100, color="black", fill="#0072B2") +
  geom_vline(aes(xintercept=mean(miles_travelled)), color="black", linetype="solid", size=1) +
  scale_x_continuous(breaks = c(0,500,1000,1500,2000,2500)) +
  xlab(" Order miles traveled ") +
  ylab(" Order count ") +
  labs(title="") +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) 

dev.off()

#---------------------------------------------------------------#

#create sacrificial data frame
overflowdays = 10
data_sac <- data
data_sac$delivery_time[data_sac$delivery_time > overflowdays ] <- overflowdays

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 14))

dev.new(width=8, height=4)

png(file="../ordertimes.png",
    width=700, height=650)

#Make the histogram
ggplot(data_sac, aes(x=delivery_time)) +
  geom_histogram(binwidth=0.5, color="black", fill="#0072B2") +
  geom_vline(aes(xintercept=mean(delivery_time)), color="black", linetype="solid", size=1) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  xlab(" Order delivery time (days) ") +
  ylab(" Order count ") +
  labs(title="") +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) 

dev.off()
