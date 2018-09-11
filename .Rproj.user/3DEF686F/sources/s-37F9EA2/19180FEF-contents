# 1_calculate_gradient

library(dplyr)
library(ggplot2)
library(ggrepel)

# read in chemisty data to calculcate AMD gradient
# Cond = conductivity
# element concentrations in ug/l
raw.dat <- read.csv("data/gradient_chemistry.csv",
                stringsAsFactors = FALSE,
                row.names = 1)

dat <- raw.dat
# log10 (x +1) transform Cond and element variables
dat[,2:13] <- log10(dat[,2:13] + 1)

# PCA to determine gradient
gradient <- prcomp(dat, center = T, scale. = T)

# figure 1 in MS ####
# function to make biplot and have labels not overlap
PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x),
    PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) +
    geom_point(shape = 2, size = 2)
  datapc <- data.frame(varnames=rownames(PC$rotation),
                       PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/
       (max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/
       (max(datapc[,x])-min(datapc[,x])))
    )
  datapc <- transform(datapc,
                      v1 = 1 * mult * (get(x)),
                      v2 = 0.5 * mult * (get(y))
  )
  x_limits = c(1,NA)
  plot <- 
    plot +
    coord_equal() +
    geom_text_repel(data=datapc,
              aes(x=v1,
                  y=v2,
                  label=varnames),
              size = 5, color="black",
              direction = "y",
              nudge_x = 0.8,
              segment.alpha = 0.1,
              box.padding = 0.2,
              point.padding = NA,
              hjust = 1
              ) +
    geom_segment(data=datapc,
                 aes(x=0, y=0, xend=v1, yend=v2),
                 arrow=arrow(length=unit(0.3,"cm")),
                 alpha=1, color="black") +
    theme_classic()
  plot
}

tiff("results/fig_1.tiff",
     units="mm", width=200, height=100, res=300)
PCbiplot(gradient)
dev.off()

# correlation of chemical variables with PC1
data.frame(variable = rownames(gradient$rotation),
           pca1 =cor(dat,
                     gradient$x[,1])) %>%
  arrange(pca1)

# PC1 axis for future analyses
pca.axis <- data.frame(site = row.names(dat),
                       pca1= gradient$x[,1],
                       stringsAsFactors = FALSE) %>%
  arrange(pca1)

# save pca.axis
write.csv(pca.axis, "results/pca1_axis.csv",
          row.names = FALSE)

# table 1 in MS
raw.dat$site <- row.names(raw.dat)
tab1 <- left_join(pca.axis, raw.dat, by = "site")
names(tab1)[2] <- "PC1"
# convert elemental concentrations from ug/l to mg/l
tab1[,5:15] <- tab1[,5:15] / 1000 
write.csv(tab1, "results/table_1.csv",
          row.names = FALSE)

