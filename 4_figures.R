# 4_figures

library(dplyr)
library(ggplot2)
library(gridExtra)

# read in binned size spectra data
binned <- read.csv("results/binned_size_spectra.csv",
                   stringsAsFactors = FALSE)
# read in top quadratic model
mod <- readRDS("results/top_quadratic_model.RDS")
# read in mrange data
mrange <- read.csv("results/mrange_data.CSV")

# fig_2 for MS ####
# make new data frame for prediction
# e.g. fit line
newdata <- binned %>% 
  select(site, pca1, log_mids_center,
         log_count_corrected, log_mids)

# predict values with "predict" interval and cbind to newdata 
newdata <- cbind(newdata,
                 predict(mod,
                         newdata,
                         interval = "predict"))

# figure 2 in MS ####
tiff("results/fig_2.tiff",
     units="mm", width=180, height=150, res=300)
ggplot(binned, aes(y = log_count_corrected, x = log_mids)) + 
  labs(x = expression(Log[10]~M), y = expression(Log[10]~N)) +
  facet_wrap(~reorder(site, pca1)) +
  geom_point(size = 0.7) +
  geom_line(data = newdata, aes(y=fit, x = log_mids, group = site), 
            size = 0.5) +
  geom_ribbon(data = newdata, aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 8, margin = margin(.09, 0, .09, 0, "cm"))) 
dev.off()


# figure 3 ####
# slope coefficient
# make new dataframe for calculation
slopedata <- data.frame(
  pca1 = unique(binned$pca1),
  log_mids_center = 1) #mid point of PCA gradient

# coefficients of model
c.mod <- coef(mod)

#predict values (with CI) for new data
etta <- data.frame(
  etta = predict(mod, slopedata, interval = "confidence"),
  pca1 = slopedata$pca1,
  slope.x = slopedata$log_mids_center)

#mutate predicted values (+/- CI)
#formula from talk with D Stouffer 17 jan 17
#(B + EP) = (etta - (g + rP) - sC^2 ) / C
slope.CI <- etta %>%
  mutate(estimate = (etta.fit - (c.mod["(Intercept)"] + 
                                   c.mod["pca1"] * pca1) - 
                       c.mod["I(log_mids_center^2)"] * slope.x^2) / slope.x,
         lo = (etta.lwr - (c.mod["(Intercept)"] + c.mod["pca1"]*pca1) - 
                 c.mod["I(log_mids_center^2)"] * slope.x^2) / slope.x,
         hi = (etta.upr - (c.mod["(Intercept)"] + c.mod["pca1"]*pca1) - 
                 c.mod["I(log_mids_center^2)"] * slope.x^2) / slope.x) %>%
  select(pca1, estimate, lo, hi)

tiff(filename = "results/fig_3.tiff", width = 100, height = 80,
    units = "mm", res =300)
ggplot(slope.CI, aes(y = estimate, x = pca1)) +
  geom_line(color = "black", size = 1) + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_point(aes(x = pca1, y = estimate), size = 2) +
  theme_classic()+ 
  labs(y = "Estimated slope \ncoefficient", x = "Mining gradient")
dev.off()

# table 4 ####
binned %>%
  select(site, pca1) %>%
  unique() %>%
  mutate(Intercept = c.mod[1] + c.mod[4]*pca1,
         MPC1 = c.mod[5]*pca1 + c.mod[2],
         M2 = c.mod[3]) %>%
  arrange (pca1) %>%
  write.csv("results/table_4.png")

# figure 4 in MS ####
fig4a <- ggplot(mrange, aes(y = mrange, x = pca1)) +
  geom_point() +
  stat_smooth(method = "lm", colour = "black") +
  theme_classic() +
  labs(y = expression(Range~of~Log[10]~M),
       x = "Mining gradient") +
  annotate(geom = "text", x = 5, y = 4.75, label = "A", size = 6)

fig4b <- binned %>% 
  ggplot(aes(y = log_mids, x = pca1)) +
  geom_point() +
  theme_classic() +
  geom_quantile(quantiles = c(0.05, 0.95),
                color = "black", size = 1)+
  annotate(geom = "text", x = 5, y = -1, label = "B", size = 6) +
  labs(y = expression(Log[10]~M),
       x = "Mining gradient")

tiff(filename = "results/fig_4.tiff", width = 225,
    height = 100, units = "mm", res =300)
grid.arrange(fig4a, fig4b, ncol = 2)
dev.off()


# figure 5 ####
dataset <- read.csv("data/estimated_dw.csv",
                    stringsAsFactors = FALSE)
pca.axis <- read.csv("results/pca1_axis.csv",
                     stringsAsFactors = FALSE)
# join pca.axis and data set
dataset <- left_join(dataset, pca.axis, by = "site") %>%
  arrange(pca1)
# abbreviate site names for visual display
dataset$mod.site <- recode(dataset$site,
                      "Sullivan East" = "Sull.E",
                      "Burnet Tributary" = "Burn.Trib",
                      "One Horse Creek" = "One Horse",
                      "Coalbrookedale" = "Cbdale",
                      "Sullivan West" = "Sull.W")

ffg.data <- dataset %>%
  filter(taxa != "Acari" & !is.na(FFG)) %>%
  group_by(mod.site, FFG) %>%
  summarize(n = n(), pc1 = first(pca1)) %>%
  mutate(tot.ab = sum(n),
         rel.ab = n / tot.ab)

tiff(filename = "results/fig_5.tiff",
    width = 420, height = 200, units = "mm", res =300)
ggplot(ffg.data, aes(x = reorder(mod.site, pc1), 
                     fill = FFG)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9",
                               "#009E73", "#0072B2", "#D55E00",
                               "#CC79A7")) +
  theme_classic(base_size = 20)+
  labs(x = "Site", y = "Proportion")+
  theme(axis.text.x=element_text(size = 17,
                                 hjust = 1,
                                 vjust = 0.5,
                                 angle = 90,
                                 color = "black"))
dev.off()

