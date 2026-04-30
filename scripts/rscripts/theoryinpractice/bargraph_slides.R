#install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(dplyr)
library(ggplot2)

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


labels <- c(0.0, 0.3, 0.6, 0.9)
ptp <- c(1270, 1339, 1481, 1615)
relay <- c(1121, 1119, 1118, 1122)
ptp_repositioning <- c(1281, 1292, 1324, 1409)


df <- data.frame(
  label = factor(rep(labels, 2)),
  value = c(ptp, relay),
  group = rep(c("Point-to-point", "Relay"), each = length(labels))
)

df3 <- data.frame(
  label = factor(rep(labels, 3)),
  value = c(ptp, ptp_repositioning, relay),
  group = rep(c("Point-to-point", "Point-to-point with repositioning", "Relay"), each = length(labels))
)


update_geom_defaults("text", list(size = 24))

dev.new(width=8, height=4)

png(file="theoryinpractice_slides.png",
    width=3000, height=2400)

ggplot(df, aes(x = label, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "\n Edge-level imbalance", y = "Kilometers per order \n ", fill = "Legend") +
  scale_fill_manual(values = c("#648FFF", "#FE6100")) +
  theme_linedraw() + 
  theme(legend.position="none", legend.text = element_text(size = 130)) +
  theme(axis.text = element_text(size = 90)) +
  theme(axis.title = element_text(size = 130))

dev.off()



