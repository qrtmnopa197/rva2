#SOURCE THE POWER ANALYSIS FUNCTION BEFORE RUNNING THIS

library(tidyverse)
library(nlme)
library(MuMIn)
library(emmeans)
library(sigmoid)
library(future)
library(future.apply)

# Set up parallel processing
plan(multisession)

disorg_power <- future_lapply(1:10, power_rva_2, sample_sizes=c(55,70,85),sims=1000)
dp_df <- do.call(rbind,disorg_power)
power_df <- dp_df %>%
  group_by(sample_size) %>%
  summarize(across(everything(), mean))
view(power_df)





