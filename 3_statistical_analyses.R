# 3_statistical_analyses

library(MuMIn)
library(plyr)
library(dplyr)
library(quantreg)
library(ggplot2)

# read in binned size spectra data
binned <- read.csv("results/binned_size_spectra.csv",
                   stringsAsFactors = FALSE)

# global models ####
# full quadratic
global.quadratic <- lm(log_count_corrected~ 
                         log_mids_center *pca1 +
                         I(log_mids_center^2) + 
                         I(log_mids_center^2):pca1,
                       data=binned,
                       na.action = "na.fail")
# full linear model
global.linear <- lm(log_count_corrected ~ 
                      log_mids_center*pca1,
                    data=binned,
                    na.action = "na.fail")
# compare quadratic and linear models
AIC(global.quadratic, global.linear)
# quadratic much better than linear
# quad AIC = 407.179
# linear AIC = 528.801
# move forward with quad model

# systematically test simplified quadratic models using MuMIn::dredge 
dredge.models <- dredge(global.quadratic,
                    beta = "none",
                    extra = "R^2")

# table 2 for MS #### 
# table with results for all simplified models
sink("results/table_2.txt")
dredge.models
sink()

# pick best model based on AICc using MuMIn::get.models
top.models <- get.models(dredge.models,
                         subset = delta < 2)
# single top model selected

# save top.model for making figure in script 4
saveRDS(top.models[[1]],
        "results/top_quadratic_model.RDS")


# Mrange ####
# Mrange ####
# calculate range of M values
mrange <- binned %>% 
  group_by(site, pca1) %>%
  summarise(mrange = max(log_mids_center) - min(log_mids_center))
# save Mrange as csv
write.csv(mrange, "results/mrange_data.CSV",
          row.names = FALSE)

# mrange ~ gradient linear model
m.mod <- lm(mrange ~ pca1, data = mrange)
summary(m.mod)

# Quantile regression
M.quant <- rq(log_mids~pca1, data = binned, tau = c(0.05, 0.95))
summary(M.quant)


