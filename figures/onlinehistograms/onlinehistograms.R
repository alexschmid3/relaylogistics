install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(leaflet)
library(lubridate)
library(dplyr)
library(network)
library(outliers)
library(geosphere)
library(fields)
library(hash)
library(naniar)

#---------------------------------------------------------------#

#Read in data
#abcgdata <- read_csv('convergence_ex2_abcg-ext_finalcomp1_rundate2022-03-21.csv')
#pbcgdata <- read_csv('convergence_ex2_pbcg_finalcomp1_rundate2022-03-21.csv')

abcgdata <- read_csv('convergence/convergence_ex2_abcg_conv_1_rundate2022-03-21.csv')
pbcgdata <- read_csv('convergence/convergence_ex2_pbcg_conv_25_rundate2022-03-22.csv')
onedata <- read_csv('convergence/convergence_ex2_abcg_one_1_rundate2022-03-25.csv')

abcgdata <- read_csv('convergence/convergence_ex4_abcg_BTex4_15_rundate2022-03-30.csv')
pbcgdata <- read_csv('convergence/convergence_ex4_pbcg_BTpbcg_71_rundate2022-03-31.csv')
onedata <- read_csv('convergence/convergence_ex4_abcg_one_55_rundate2022-03-25.csv')

abcgdata <- read_csv('convergence/convergence_ex4_abcg_BTex4_9_rundate2022-03-30.csv')
pbcgdata <- read_csv('convergence/convergence_ex4_pbcg_BTpbcg_65_rundate2022-03-31.csv')
onedata <- read_csv('convergence/convergence_ex4_abcg_one_49_rundate2022-03-25.csv')

abcgdata <- read_csv('convergence/convergence_ex4_abcg_BTex4_10_rundate2022-03-30.csv')
pbcgdata <- read_csv('convergence/convergence_ex4_pbcg_BTpbcg_66_rundate2022-03-31.csv')
onedata <- read_csv('convergence/convergence_ex4_abcg_one_50_rundate2022-03-25.csv')

#---------------------------------------------------------------#

joineddata <- full_join(abcgdata, pbcgdata, by = "pbcgiteration") %>%
  select(pbcgiteration, arc_count.x, arc_count.y, path_count.x, path_count.y,lowerbound.x, lowerbound.y, upperbound.x, upperbound.y)

joineddata <- rename(joineddata, arc_count_mvg = arc_count.x)
joineddata <- rename(joineddata, arc_count_cg = arc_count.y)
joineddata <- rename(joineddata, path_count_mvg = path_count.x)
joineddata <- rename(joineddata, path_count_cg = path_count.y)
joineddata <- rename(joineddata, lowerbound_mvg = lowerbound.x)
joineddata <- rename(joineddata, lowerbound_cg = lowerbound.y)
joineddata <- rename(joineddata, upperbound_mvg = upperbound.x)
joineddata <- rename(joineddata, upperbound_cg = upperbound.y)

joineddata <- full_join(joineddata, onedata, by = "pbcgiteration") %>%
  select(pbcgiteration, arc_count_mvg, arc_count_cg, arc_count, path_count_mvg, path_count_cg, path_count, lowerbound_mvg, lowerbound_cg, lowerbound, upperbound_mvg, upperbound_cg, upperbound)

joineddata <- rename(joineddata, arc_count_svg = arc_count)
joineddata <- rename(joineddata, path_count_svg = path_count)
joineddata <- rename(joineddata, lowerbound_svg = lowerbound)
joineddata <- rename(joineddata, upperbound_svg = upperbound)

#---------------------------------------------------------------#

update_geom_defaults("text", list(size = 10))

dev.new(width=8, height=4)

#Combined path and arc
ggplot(joineddata, aes(x=pbcgiteration)) + 
  geom_ribbon(aes(ymin = 0, ymax = path_count_svg, fill = "SVG paths (induced)"), alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = path_count_cg, fill = "PBCG paths (variables)"), alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = path_count_mvg, fill = "MVG paths (induced)"), alpha = 0.3) +
  geom_line(aes(y = arc_count_svg, colour = "SVG arcs (variables)"), size = 2) +
  geom_line(aes(y = arc_count_cg, colour = "PBCG arcs (induced)"), size = 2) +
  geom_line(aes(y = arc_count_mvg, colour = "MVG arcs (variables)"), size = 2) +
  scale_colour_manual("", 
                      breaks = c("MVG arcs (variables)", "PBCG arcs (induced)", "SVG arcs (variables)"),
                      values = c("#E27532", "#2D5C87", "#DC267F")) +
  scale_fill_manual("",breaks = c("MVG paths (induced)", "PBCG paths (variables)", "SVG paths (induced)"),
                    values = c("#E27532", "#2D5C87", "#DC267F")) +
  xlab(" Iteration ") +
  ylab(" Count ") +
  labs(title="") +
  #theme(legend.position = c(0.4,.9), legend.direction = "horizontal")
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18))+
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))

#Convergence plot
ggplot(joineddata, aes(x=pbcgiteration)) + 
  geom_line(aes(y = lowerbound_svg, colour = "SVG lower bound"), size=1.2) +
  geom_line(aes(y = lowerbound_cg, colour = "PBCG lower bound"), size=1.2) +
  geom_line(aes(y = lowerbound_mvg, colour = "MVG lower bound"), size=1.2) +
  geom_line(aes(y = upperbound_svg, colour = "SVG upper bound"), size=2) +
  geom_line(aes(y = upperbound_cg, colour = "PBCG upper bound"), size=2) +
  geom_line(aes(y = upperbound_mvg, colour = "MVG upper bound"), size=2) +
  scale_colour_manual("", 
                      breaks = c("MVG upper bound", "MVG lower bound", "SVG upper bound", "SVG lower bound", "PBCG upper bound", "PBCG lower bound"),
                      values = c("#E27532", "#E2A47E", "#DC267F", "#D884AD", "#2D5C87", "#649BCB")) +
  xlab(" Iteration ") +
  ylab(" Objective value ") +
  ylim(-1000000,1000000) +
  labs(title="") + 
  #theme(legend.position = c(0.4,1.0), legend.direction = "horizontal")
  theme(legend.position="bottom")+
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18))  
  #guides(color=guide_legend(nrow=2, byrow=TRUE))

#---------------------------------------------------------------#

#Read in tstep data
tstepdata <- read_csv('tstep/processeddata.csv')
tstepdata <- rename(tstepdata, solvetime = time)

#---------------------------------------------------------------#

options(ggplot2.discrete.color = c("orange", "#DA853E"))
update_geom_defaults("text", list(size = 30))

dev.new(width=7, height=4)

#tstep graph - combined

png(file="figures/tstep.png",
width=600, height=450)

ggplot(tstepdata, aes(x=as.character(tstep), y=obj/36000, group = 1)) + 
  geom_bar(aes(x=as.character(tstep), y=solvetime/3600, fill = "CPU time"), stat="identity", size = 1) +
  geom_line(aes(y = obj/36000, colour = "Objective value"), size = 2) +
  geom_point(aes(y = obj/36000, colour = "Objective value"), size = 5) + # Show dots
  geom_point(aes(y = miles/36000, colour = "Total miles"), size = 5) + # Show dots
  geom_label(aes(label=optgap), nudge_y = 0.2, size=5) +
  geom_label(aes(label=milesgap), nudge_y = -0.44, size=5) +
  geom_line(aes(y = miles/36000, colour = "Total miles"), size = 2) +
  scale_colour_manual("", 
                      breaks = c("Objective value", "Total miles"),
                      values = c("#2D5C87", "#649BCB")) +
  scale_fill_manual("", 
                      breaks = c("CPU time"),
                      values = c("#DA853E")) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Time (hours)",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*36000, name="Cost")
  ) +
  xlab(" Time discretization (hours) ") +
  labs(title="") +
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18))

#---------------------------------------------------------------#

#tstep graph - objective
ggplot(tstepdata, aes(x=as.character(tstep), y=obj, group = 1)) + 
  geom_line(aes(y = obj, colour = "Objective value"), size = 2) +
  geom_point(aes(y = obj, colour = "Objective value"), size = 5) + # Show dots
  geom_point(aes(y = miles, colour = "Total miles"), size = 5) + # Show dots
  geom_label(aes(label=optgap), nudge_y = 10000) +
  geom_line(aes(y = miles, colour = "Total miles"), size = 2) +
  scale_colour_manual("", 
                      breaks = c("Objective value", "Total miles"),
                      values = c("#2D5C87", "#649BCB")) +
  xlab(" Time discretization (hours) ") +
  ylim(0,100000) +
  ylab(" Vehicle miles traveled ") +
  labs(title="") +
  theme(legend.position="bottom", legend.text = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 20))  

#tstep graph - time
ggplot(tstepdata, aes(x=as.character(tstep), y=time, group = 1)) + 
  geom_line(aes(y = time/3600, colour = "CPU time"), size = 2) +
  geom_point(aes(y = time, colour = "CPU time")) + # Show dots
  scale_colour_manual("", 
                      breaks = c("CPU time"),
                      values = c( "#DA853E")) +
  xlab(" Time discretization (hours) ") +
  ylab(" CPU time (hours)") +
  labs(title="") +
  theme(legend.position="bottom", legend.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 15))  

#---------------------------------------------------------------#

#Read in detour and delay distribution data
abcgdata <- read_csv('orderhist_abcg.csv')
otddata <- read_csv('orderhist_otd.csv')
rrdata <- read_csv('orderhist_rr.csv')

#---------------------------------------------------------------#

joineddata <- full_join(abcgdata, otddata, by = "orderid") %>%
  select(orderid, orderdetour.x, orderdetour.y)

joineddata <- rename(joineddata, detour_abcg = orderdetour.x)
joineddata <- rename(joineddata, detour_otd = orderdetour.y)

joineddata <- full_join(joineddata, rrdata, by = "orderid") %>%
  select(orderid, detour_abcg, detour_otd, orderdetour)
         
joineddata <- rename(joineddata, detour_rr = orderdetour)

vec_mvg <- rep("MVG", 831)
vec_otd <- rep("Orders-then-drivers", 831)

vec_rr <- rep("Actual routes", 831)

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

alldata = rbind(mvghist, otdhist) #, rrhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("MVG", "Orders-then-drivers", "Actual routes"), 
                            labels=c("MVG", "Orders-then-drivers", "Actual routes"))

alldata <- subset(alldata, !is.na(orderdetour))

mean_detour <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(orderdetour))

#---------------------------------------------------------------#

#Detour
png(file="plots/detourdensity.png",
    width=600, height=450)

alldata %>%
  ggplot(aes(x=orderdetour, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3,size=1)+ 
  geom_vline(data = mean_detour, aes(xintercept = mean, 
                                         color = datalabel), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("MVG", "Orders-then-drivers"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("MVG", "Orders-then-drivers"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Detour (as % of shortest path)", y="Density")+
  scale_x_continuous(limits=c(0,0.15), labels = scales::percent_format(accuracy = 1)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 26))

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
vec_rr <- rep("Actual routes", 831)

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

alldata = rbind(mvghist, otdhist) #, rrhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("Optimized solution", "Orders-then-drivers", "Actual routes"), 
                            labels=c("Optimized solution", "Orders-then-drivers", "Actual routes"))

alldata <- subset(alldata, !is.na(orderdelay))

mean_delay <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(orderdelay))

#---------------------------------------------------------------#

#Delay
png(file="plots/delaydensity.png",
    width=450, height=600)

alldata %>%
  ggplot(aes(x=orderdelay, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3, size=1)+ 
  geom_vline(data = mean_delay, aes(xintercept = mean, 
                                     color = datalabel), size=1.5)+
  scale_colour_manual("", breaks = c("Optimized solution", "Orders-then-drivers"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Orders-then-drivers"),
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

#Read in driver empty distribution
abcgdata <- read_csv('online/driverempty_abcg.csv')
otddata <- read_csv('online/driverempty_otd.csv')

#---------------------------------------------------------------#

vec_mvg <- rep("Optimized solution", 2998)
vec_otd <- rep("Orders-then-drivers solution", 2337)

mvghist <- abcgdata %>% select(driverid, emptymiles, percentempty)
otdhist <- otddata %>% select(driverid, emptymiles, percentempty)

mvghist <- mvghist %>%
  mutate(datalabel = vec_mvg)
otdhist <- otdhist %>%
  mutate(datalabel = vec_otd)

alldata = rbind(mvghist, otdhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("Optimized solution", "Orders-then-drivers solution"), 
                            labels=c("Optimized solution", "Orders-then-drivers solution"))

alldata <- subset(alldata, !is.na(percentempty))
alldata <- subset(alldata, !is.na(emptymiles))

mean_pct <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(percentempty))

mean_empty <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(emptymiles))

mean_pct["mean"][mean_pct["datalabel"] == "Orders-then-drivers solution"] <- 0.24
mean_pct["mean"][mean_pct["datalabel"] == "Optimized solution"] <- 0.13

#---------------------------------------------------------------#

#Driver empty percent
png(file="plots/driveremptypct.png",
    width=450, height=600)

alldata %>%
  ggplot(aes(x=percentempty, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3,size=1)+ 
  geom_vline(data = mean_pct, aes(xintercept = mean, 
                                    color = datalabel), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("Optimized solution", "Orders-then-drivers solution"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Orders-then-drivers solution"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Empty miles per driver\n(% of miles)", y="Density")+
  scale_x_continuous(limits=c(0,0.5), labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position="none") +
  #theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20))

dev.off()


#Driver empty miles
png(file="plots/driverempty.png",
    width=450, height=600)

alldata %>%
  ggplot(aes(x=emptymiles, color=datalabel, fill=datalabel)) +
  geom_density(alpha=0.3,size=1)+ 
  geom_vline(data = mean_empty, aes(xintercept = mean, 
                                    color = datalabel), size=1.5)+
  #xlim(0, 0.00000001)+
  #ylim(0, 100)+
  scale_colour_manual("", 
                      breaks = c("Optimized solution", "Orders-then-drivers solution"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Orders-then-drivers solution"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Empty miles per driver", y="Density")+
  #scale_x_continuous(limits=c(0,1), labels = scales::percent_format(accuracy = 1)) +
  #theme(legend.position="none") +
  theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 26))

dev.off()

#---------------------------------------------------------------#

#Read in remaining miles data
abcgdata <- read_csv('online/orderoutcomes_abcg.csv')
otddata <- read_csv('online/orderoutcomes_otd.csv')

#---------------------------------------------------------------#

vec_mvg <- rep("Optimized solution", 831)
vec_otd <- rep("Orders-then-drivers solution", 831)

mvghist <- abcgdata %>% select(orderid, remainingmilespct)
otdhist <- otddata %>% select(orderid, remainingmilespct)

mvghist <- rename(mvghist, remainingmiles = remainingmilespct)
otdhist <- rename(otdhist, remainingmiles = remainingmilespct)

mvghist <- mvghist %>%
  mutate(datalabel = vec_mvg)
otdhist <- otdhist %>%
  mutate(datalabel = vec_otd)

alldata = rbind(mvghist, otdhist)

alldata$datalabel <- factor(alldata$datalabel, 
                            levels=c("Optimized solution", "Orders-then-drivers solution"), 
                            labels=c("Optimized solution", "Orders-then-drivers solution"))

alldata <- subset(alldata, !is.na(remainingmiles))

mean_rem <- alldata %>%
  group_by(datalabel) %>%
  summarize(mean=mean(remainingmiles))

#---------------------------------------------------------------#

#Order remaining miles
png(file="plots/remainingmilespct.png",
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
                      breaks = c("Optimized solution", "Orders-then-drivers solution"),
                      values = c("#DA853E", "#649BCB"))+
  scale_fill_manual("", breaks = c("Optimized solution", "Orders-then-drivers solution"),
                    values = c("#DA853E", "#649BCB")) +
  labs(x="Remaining miles per non-completed order \n(% of trip remaining)", y="Density")+
  scale_x_continuous(limits=c(0,1), labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position="none") +
  #theme(legend.position="bottom", legend.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20))

dev.off()

