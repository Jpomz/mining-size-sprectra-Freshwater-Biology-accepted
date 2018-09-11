# 5_Supplementary_material

library(plyr)
library(dplyr)
library(ggplot2)
library(MuMIn)

# number of bins
# function to bin and center data
bin_and_center <- function(data, var, breaks, ...){
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks, # need to predefine breaks
                     # e.g. Log2 breaks = 2^seq(min, max) 
                     # Log10 breaks = 10^seq(min, max)
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  break_width = breaks_offset - breaks_orig
  count = binned_hist$counts 
  dataout = data.frame(
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / break_width),
    # original midpoint of bin log10 transformed
    log_mids = log10(binned_hist$mids))
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # recenter data at x=0
  mid_row = ceiling(nrow(dataout)/2)
  # subtract value of mid row from all mids
  dataout$log_mids_center = 
    dataout[,"log_mids"] - dataout[mid_row,"log_mids"]  
  dataout
}

dataset <- read.csv("data/estimated_dw.csv",
                    stringsAsFactors = FALSE)
pca.axis <- read.csv("results/pca1_axis.csv",
                     stringsAsFactors = FALSE)

# log2 breaks ####
breaks2 <- 2^seq(-22,-2)
binned2 <- ddply(dataset, .(site),
                 bin_and_center, var = "dw", breaks = breaks2)
# adding pca axes
binned2 <- merge(pca.axis, binned2, by = "site")

# number of bins
binned2 %>% group_by(site) %>%
  summarize(no.bins = n(),
            min.bin = min(log_mids),
            max.bin = max(log_mids)) %>%
  as.data.frame
# 18 total bins; 6 - 16 unique bins per site
# which bins are not empty?

# fig S1 ####
png(filename = "results/fig_S1.png",
    width = 420, height = 200, units = "mm", res =300)
binned2 %>%
  filter(log_count_corrected != 0 ) %>%
  ggplot((aes(x = reorder(site, pca1),
              y = log_mids))) +
  geom_point(size = 4) +
  theme_classic(base_size = 17) +
  labs(y = "Bin",
       x = "Site") +
  theme(axis.text.x=element_text(size = 12,
                                 hjust = 1,
                                 angle = 90))
dev.off()

# log10 breaks ####
breaks10 <- 10^seq(-8, 1)
binned10 <- ddply(dataset, .(site),
                  bin_and_center, var = "dw", breaks = breaks10)
binned10 <- merge(pca.axis, binned10, by = "site")
# number of bins
binned10 %>% group_by(site) %>%
  summarize(no.bins = n(),
            min.bin = min(log_mids),
            max.bin = max(log_mids)) %>%
  as.data.frame
# 6 total bins; 2-6 unique bins per site
# which bins are not empty?
# figure S2 ####
png(filename = "results/fig_S2.png",
    width = 420, height = 200, units = "mm", res =300)
binned10 %>%
  filter(log_count_corrected != 0 ) %>%
  ggplot((aes(x = reorder(site, pca1),
              y = log_mids))) +
  geom_point(size = 4) +
  theme_classic(base_size = 17) +
  labs(y = "Bin",
       x = "Site") +
  theme(axis.text.x=element_text(size = 12,
                                 hjust = 1,
                                 angle = 90))
dev.off()
# compare models using log2 / log10 width bins
# global models
quad2 <- lm(log_count_corrected~ 
              log_mids_center *pca1 +
              I(log_mids_center^2) + 
              I(log_mids_center^2):pca1,
            data=binned2,
            na.action = "na.fail")
quad10 <- lm(log_count_corrected~ 
               log_mids_center *pca1 +
               I(log_mids_center^2) + 
               I(log_mids_center^2):pca1,
             data=binned10,
             na.action = "na.fail")
# information in Table S1
# dredge models
quad2.dredge <-  dredge(quad2,
                    beta = "none",
                    extra = "R^2")
# quad2.dredge has coefficient estimates, adj R2, delta AIC and model weight
top.quad2 <- get.models(quad2.dredge,
                        subset = delta < 2)
# look at summary of top model for signifigance of coefficients
summary(top.quad2[[1]])
# log 10
quad10.dredge <- dredge(quad10,
                    beta = "none",
                    extra = "R^2")
# quad10.dredge has coefficient estimates, adj R2, delta AIC and model weight
top.quad10 <- get.models(quad10.dredge,
                         subset = delta < 2)
# look at summary of top models for signifigance of coefficients
summary(top.quad10[[1]])
summary(top.quad10[[2]])

# fish presence / absence ####
unim <- binned2 %>% filter(pca1 < -1.5) %>% as.tbl()
fish.sites <- data.frame(site = c("Italian", "Coorang",
                                  "Burke", "Lankey",
                                  "Murray", "Kiwi"),
                         fish = "Y", 
                         stringsAsFactors = FALSE)
unim <- merge(unim, fish.sites, all.x = TRUE)
unim$fish[is.na(unim$fish)] <- "N"

# does fish significantly impact size spectra?
summary(lm(formula = log_count_corrected ~ 
             log_mids_center*fish + I(log_mids_center^2),
           data = unim, na.action = "na.fail"))
# fishY factor not significant (p > 0.05)
# log_mids_center:fishY also non-significant (p > 0.05)
# therefore, the presence of fish does not change the intercept, or the slope of the size spectra, respectively

# figure S3 ####
png(filename = "results/fig_S3.png",
    width = 420, height = 200, units = "mm", res =300)
ggplot(unim, aes(y = log_count_corrected,
                 x = log_mids_center,
                 shape = fish,
                 color = fish))+
  geom_jitter(size = 3,
              width = 0.05) +  
  scale_color_manual(values = c("black", "red")) +
  stat_smooth(alpha = 0,
              method = "lm",
              formula = y~x + I(x^2)) +
  labs(x = "Log10 M", y = "Log10 N") +
  theme_classic(base_size = 20)
dev.off()

# linear vs quadratic (log2 bins)
# global linear model
lin2 <- lm(log_count_corrected ~ 
             log_mids_center*pca1,
           data=binned2,
           na.action = "na.fail")
# quad 2 from above
AIC(lin2, quad2)
# linear model AIC = 615
# quadratic model AIC = 473
lin2.dredge <- dredge(lin2, 
                  beta = "none",
                  extra = "R^2")
top.lin2 <- get.models(lin2.dredge,
                       subset = delta < 2)
# information in table S2
# summary of top linear model
summary(top.lin2[[1]])
# summary of top quad model
summary(top.quad2[[1]])

# compare intercept and density ####
dens.summ <- dataset %>% group_by(site, surber) %>%
  count() 
dens.summ <- left_join(dens.summ, pca.axis, by = "site")

# relationship between density and mining gradient?
summary(lm(log10(n) ~ pca1, data = dens.summ))

# table S3 ####
# calc max, mean, min weight of dominant (> 0.1 relative abundance) taxa
dom.taxa <- dataset %>%
  filter(taxa != "Acari", taxa!= "Plecoptera",
         taxa != "Ephemeroptera", taxa != "Diptera") %>%
  group_by(site, taxa) %>%
  summarize(n = n(),
            
            ffg = first(FFG), 
            max.dw.mg = max(dw*1000),
            mean.dw.mg = mean(dw*1000),
            min.dw.mg = min(dw*1000)) %>%
  mutate(tot.ab = sum(n),
         rel.ab = n / tot.ab) %>%
  filter(rel.ab > 0.1)

dom.taxa <- left_join(dom.taxa, pca.axis, by = "site") %>% 
   arrange(pca1)

dom.taxa %>%
  select(site, pca1, taxa, rel.ab, ffg,
         max.dw.mg, mean.dw.mg, min.dw.mg) %>%
  write.csv("results/table_S3.csv",
            row.names = FALSE)
