# making data for aquaSync project

# this script organizes the data into size and site metadata

# justin pomeranz, October 2023

library(tidyverse)

# size data
dat <- read_csv("data/estimated_dw.csv")
names(dat)

# rename columns and add missing columns for data format
dat_format <- dat %>%
  mutate(sampling_method = "surber",
         sampling_area = 0.06,
         body_weight_units = "g",
         count = 1,
         multiplier = 1 / 0.06,
         organism_group = "Invertebrate"
  ) %>%
  separate(surber, into = c("surber","sample"),1) %>%
  rename(
    taxon = taxa,
    body_mass = dw,
    body_length = linear_meas,
    body_length_units = Unit
  ) %>%
  select(
    site,
    sampling_method,
    sample,
    sampling_area,
    organism_group,
    taxon,
    body_mass,
    body_length,
    body_weight_units,
    body_length_units,
    count,
    multiplier)
dat_format

dat_format <- dat_format %>%
  mutate(year = 2016,
         month = 01) %>%
  unite(site_date, site, year, month,
        remove = FALSE) %>%
  select(site, 
         year, 
         month, 
         site_date, 
         sampling_method,
         sample, 
         sampling_area,
         organism_group,
         taxon,
         body_mass,
         body_length,
         body_weight_units,
         body_length_units,
         count,
         multiplier)
#View(dat_format)

write_csv(dat_format, "data/aqua-sync-amd-size-data.csv")

