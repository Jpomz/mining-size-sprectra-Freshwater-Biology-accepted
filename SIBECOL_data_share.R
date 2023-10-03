# Data Sharing for SIBECOL

library(tidyverse)

dat <- read_csv("data/estimated_dw.csv")

no_spp <- dat %>%
  select(site, taxa) %>%
  distinct() %>%
  group_by(site) %>%
  summarize(count = n())

write_csv(no_spp, "data/sibecol_no_spp.csv")

names(dat)

share_name <- data.frame(sibe = c("Pomeranz_1",
                                  "Pomeranz_2",
                                  "Pomeranz_3",
                                  "Pomeranz_4",
                                  "Pomeranz_5",
                                  "Pomeranz_6",
                                  "Pomeranz_7",
                                  "Pomeranz_8",
                                  "Pomeranz_9",
                                  "Pomeranz_10",
                                  "Pomeranz_11",
                                  "Pomeranz_12",
                                  "Pomeranz_13",
                                  "Pomeranz_14",
                                  "Pomeranz_15",
                                  "Pomeranz_16",
                                  "Pomeranz_17",
                                  "Pomeranz_18",
                                  "Pomeranz_19",
                                  "Pomeranz_20",
                                  "Pomeranz_21",
                                  "Pomeranz_22",
                                  "Pomeranz_23",
                                  "Pomeranz_24",
                                  "Pomeranz_25"),
                         site = c("Burke",
                                  "Burnet Tributary",
                                  "Burnett",
                                  "Cannel",
                                  "Charming",
                                  "Chasm",
                                  "Coal",
                                  "Coalbrookedale",
                                  "Coorang",
                                  "Granity",
                                  "Hot",
                                  "Italian",
                                  "Kiwi",
                                  "Lankey",
                                  "Miller",
                                  "Mine",
                                  "Murray",
                                  "One Horse Creek",
                                  "Packtrack",
                                  "Portal",
                                  "Stream X",
                                  "Stream Y",
                                  "Sullivan East",
                                  "Sullivan West",
                                  "Warne"))

full_join(dat, share_name) %>% write_csv("data/sibecol_body.csv")


# location ----------------------------------------------------------------

cords <- read_csv("data/AMD_northing_easting.csv")

cords %>% 
  select(site, northing, easting) %>%
  distinct() %>% 
  write_csv("data/cords.csv")
