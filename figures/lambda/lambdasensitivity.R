#install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#

#Read in data
data <- read_csv('lambdasensitivity.csv')

#---------------------------------------------------------------#


