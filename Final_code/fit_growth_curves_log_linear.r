library(tidyverse)
rm(list = ls())
setwd("/Users/apple/Desktop/")
source("/Users/apple/Desktop/Colony counting/log_linear_functions.R")
OD_path = "co1_bigger_JC.csv"

# read in the tecan data
OD = read.table(OD_path, sep = ",", as.is = TRUE, row.names = 1)
# transpose it so wells are columns, then convert to dataframe

OD = t(OD)
OD = data.frame(OD)
names(OD)[1:3] = c("cycle","seconds","temp")

# convert to long format by "gathering"
OD = OD %>%
  gather(well, OD, -cycle, -seconds, -temp) %>%
  mutate(hour = seconds / 60 / 60)

#Plot raw growth curves
ggplot(OD,aes(x=hour,y=OD,col=well))+geom_point()

# fit piece-wise log-linear
growth_rate = OD %>%
  group_by(well) %>%
  summarize(model_fit = fit_loglinear(hour, log(OD), 200),
            fit_variable = c("growth_rate", "y0", "start", "end")) %>%
  ungroup() %>%
  pivot_wider(names_from = fit_variable, values_from = model_fit)

all_data = left_join(OD, growth_rate) %>%
  group_by(well) %>%
  mutate(OD_pred = log_linear(hour, growth_rate[1], y0[1],start[1], end[1])) %>%
  ungroup

all_data %>%
  #filter(hour < 10) %>%
  ggplot(aes(x = hour, y = exp(OD_pred), color = well))+
  #geom_line(aes(y = OD), linetype = "dashed", size = 1)+
  geom_line()+
  scale_y_log10()

write.table(growth_rate,"XX03062024_R4_final_CFP_RIF_gr.txt",sep="\t")

