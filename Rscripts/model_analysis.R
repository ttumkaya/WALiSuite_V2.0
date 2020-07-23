d1 = read.delim("AllData_cleaned/Or67d_wTSALE_reduced.txt", header=T, as.is=T, check.names = FALSE)
d2 = read.delim("AllData_cleaned/Or85d_wTSALE_reduced.txt", header=T, as.is=T, check.names = FALSE)
d3 = read.delim("AllData_cleaned/Or67d-Or85d_wTSALE_reduced.txt", header=T, as.is=T, check.names = FALSE)

####################################################################
#### Data structure -- Fly ID, two parents + one offspring per file
####################################################################
d1$Neuron = "Or67d"
d2$Neuron = "Or85d"
d3$Neuron = "Or67d-Or85d"

### replace "wind" to "air" wherever it appears
d1$`Wind status` = gsub("Wind", "Air", d1$`Wind status`)
d2$`Wind status` = gsub("Wind", "Air", d2$`Wind status`)
d3$`Wind status` = gsub("Wind", "Air", d3$`Wind status`)

### Read data and compute effect size per condition
get.effect.size = function(X, intensity, wind, satiety, sex) {
  ### filter data
  id = X$`Light Intensity(uW/mm2)` == intensity & X$`Wind status` == wind & X$Satiety == satiety & X$Sex == sex
  X = X[id, ]
  ord = order(X$`Fly ID`)
  X = X[ord, ]
  keep.id = !is.na(X$weighted_TSALE_P10)
  X = X[keep.id, ]
  ct = table(X$`Fly ID`)
  #keep.fly = names(ct)[ct == 3]
  keep.fly = names(ct)[(ct %% 3) == 0]
  keep.fly = as.numeric(keep.fly)
  X = X[X$`Fly ID` %in% keep.fly, ]
  
  ### make it into effect size
  ufly = unique(X$`Fly ID`)
  nfly = length(ufly)
  Z = X[X$Status == "Offspring", ]
  for(i in 1:nfly) {
    pid = which(X$`Fly ID` == ufly[i] & X$Status == "Parent")
    oid = which(X$`Fly ID` == ufly[i] & X$Status == "Offspring")
    Z$weighted_TSALE_P10[i] = mean(X$weighted_TSALE_P10[oid]) - mean(X$weighted_TSALE_P10[pid])
  }
  Z
}

#### By intensity
d1_14 = get.effect.size(d1, intensity="14uW", wind="NoAir", satiety="fed", sex="male")
d1_42 = get.effect.size(d1, intensity="42uW", wind="NoAir", satiety="fed", sex="male")
d1_70 = get.effect.size(d1, intensity="70uW", wind="NoAir", satiety="fed", sex="male")

d2_14 = get.effect.size(d2, intensity="14uW", wind="NoAir", satiety="fed", sex="male")
d2_42 = get.effect.size(d2, intensity="42uW", wind="NoAir", satiety="fed", sex="male")
d2_70 = get.effect.size(d2, intensity="70uW", wind="NoAir", satiety="fed", sex="male")

d3_14 = get.effect.size(d3, intensity="14uW", wind="NoAir", satiety="fed", sex="male")
d3_42 = get.effect.size(d3, intensity="42uW", wind="NoAir", satiety="fed", sex="male")
d3_70 = get.effect.size(d3, intensity="70uW", wind="NoAir", satiety="fed", sex="male")

d = rbind(d1_14, d1_42, d1_70, d2_14, d2_42, d2_70, d3_14, d3_42, d3_70)  ### append all together

####################################################################
#### Plotting the effect size per neuron and light intensity
####################################################################

pdf(sprintf("%s_boxplot.pdf", d3$Neuron))
par(mar=c(12,4,2,2)) 
bp <- boxplot(d$weighted_TSALE_P10 ~ d$Neuron * d$`Light Intensity(uW/mm2)`, las=2)
abline(h=0, lty=2, col=2)
dev.off()

####################################################################
#### Regression analysis
####################################################################
#### Get B bootstrap samples of K flies of each neuron (B = 10000) at each light intensity.
#### Sort each group of flies by their activation state (wTSALE) from the lowest to the highest
#### and then form trios of (neuron1, neuron2, neuron1+neuron2)
#### Next, run regression of Or47a-Or67d against Or47a and Or67d
####     Model: Or47a-Or67d = beta0 + beta1 * Or47a + beta2 * Or67d + beta3 * (Or47a+Or67d)

d14 = d[d$`Light Intensity(uW/mm2)` == "14uW", ]
d42 = d[d$`Light Intensity(uW/mm2)` == "42uW", ]
d70 = d[d$`Light Intensity(uW/mm2)` == "70uW", ]

neuron1 = d1$Neuron[1]
neuron2 = d2$Neuron[1]
neuron12 = d3$Neuron[1]

#### tmpd is the data
bootstrap_reg = function(X, B=10000, neuron1, neuron2, neuron12) {
  ###### Remove NaN's
  X = X[!is.na(X$weighted_TSALE_P10),]
  K = min(table(X$Neuron))
  set.seed(2222)
  beta1 = beta2 = beta3 = rep(NA,  B)
  ind = 1:nrow(X)
  
  nid1 = which(X$Neuron == neuron1)
  nid2 = which(X$Neuron == neuron2)
  nid12 = which(X$Neuron == neuron12)

  for(i in 1:B) {
    id1 = sample(nid1, K, replace=TRUE)
    id2 = sample(nid2, K, replace=TRUE)
    id12 = sample(nid12, K, replace=TRUE)
    y1 = X$weighted_TSALE_P10[id1]
    y2 = X$weighted_TSALE_P10[id2]
    y12 = X$weighted_TSALE_P10[id12]
    y1 = y1[order(y1)]
    y2 = y2[order(y2)]
    y12 = y12[order(y12)]
    tmp.fit = lm(y12 ~ y1 + y2)
    tmp.coef = coef(tmp.fit)
    beta1[i] = tmp.coef[2]
    beta2[i] = tmp.coef[3]
    beta3[i] = tmp.coef[4]
  }
  list(beta1=beta1, beta2=beta2, beta3=beta3)
}

res14 = bootstrap_reg(d14, B=10000, neuron1=d1$Neuron[1], neuron2=d2$Neuron[1], neuron12=d3$Neuron[1])
res42 = bootstrap_reg(d42, B=10000, neuron1=d1$Neuron[1], neuron2=d2$Neuron[1], neuron12=d3$Neuron[1])
res70 = bootstrap_reg(d70, B=10000, neuron1=d1$Neuron[1], neuron2=d2$Neuron[1], neuron12=d3$Neuron[1])

clean.samples = function(x, cut) x[abs(x) <= cut]

pdf(sprintf("%s_distribution.pdf", d3$Neuron))
par(mfrow=c(3,2))
hist(clean.samples(res14$beta1,5), breaks=1000, xlim=c(-2,2), main="14uW", xlab=d1$Neuron[1])
hist(clean.samples(res14$beta2,5), breaks=1000, xlim=c(-2,2), main="14uW", xlab=d2$Neuron[1])
#hist(clean.samples(res14$beta3,5), breaks=1000, xlim=c(-2,2), main="14uW", xlab="Interaction")

hist(clean.samples(res42$beta1,5), breaks=1000, xlim=c(-2,2), main="42uW", xlab=d1$Neuron[1])
hist(clean.samples(res42$beta2,5), breaks=1000, xlim=c(-2,2), main="42uW", xlab=d2$Neuron[1])
#hist(clean.samples(res42$beta3,5), breaks=1000, xlim=c(-2,2), main="42uW", xlab="Interaction")

hist(clean.samples(res70$beta1,5), breaks=1000, xlim=c(-2,2), main="70uW", xlab=d1$Neuron[1])
hist(clean.samples(res70$beta2,5), breaks=1000, xlim=c(-2,2), main="70uW", xlab=d2$Neuron[1])
#hist(clean.samples(res70$beta3,5), breaks=1000, xlim=c(-2,2), main="70uW", xlab="Interaction")
dev.off()

############### Write medians to csv #############################################

beta1_14_med <- quantile(clean.samples(res14$beta1,5))[3]
beta2_14_med <- quantile(clean.samples(res14$beta2,5))[3]

beta1_42_med <- quantile(clean.samples(res42$beta1,5))[3]
beta2_42_med <- quantile(clean.samples(res42$beta2,5))[3]

beta1_70_med <- quantile(clean.samples(res70$beta1,5))[3]
beta2_70_med <- quantile(clean.samples(res70$beta2,5))[3]

df <- matrix(c(beta1_14_med, beta2_14_med, beta1_42_med, beta2_42_med, beta1_70_med, beta2_70_med), ncol = 6, byrow = TRUE)
colnames(df) <- c('14uW/mm2_B1', '14uW/mm2_B2','42uW/mm2_B1', '42uW/mm2_B2','70uW/mm2_B1', '70uW/mm2_B2')
rownames(df) <- c(neuron12)

fname <- (sprintf("%s_regression_coef.csv", neuron12))
write.csv(df, fname)

################# Plot the combo coefficients ######################################

#### Read the B1 and B2 from csv files
library(dplyr)

setwd("C:/Users/tumkayat/Desktop/Valence_rules/Results")
temp <- list.files(path = "C:/Users/tumkayat/Desktop/Valence_rules/Results", pattern = ".csv")
dff <- lapply(temp, read.csv) %>% bind_rows()

ORN_labels <- matrix(c(rep(dff$X, 3)), ncol=1)
Intensities <- matrix(c(rep("14uW/mm2", length(dff$X)), rep("42uW/mm2", length(dff$X)), rep("70uW/mm2", length(dff$X))), ncol=1)
Beta1 <- matrix(c(dff$X14uW.mm2_B1, dff$X42uW.mm2_B1, dff$X70uW.mm2_B1), ncol=1)
Beta2 <- matrix(c(dff$X14uW.mm2_B2, dff$X42uW.mm2_B2, dff$X70uW.mm2_B2), ncol=1)

dfff <- cbind.data.frame(ORN_labels, Intensities, Beta1, Beta2)

## Calculate coefficients distance from the diagonal
# https://www.intmath.com/plane-analytic-geometry/perpendicular-distance-point-line.php
# https://stackoverflow.com/questions/54509365/how-to-calculate-a-distance-matrix-from-a-diagonal-line

eucladian_dist = abs(dfff$Beta1 - dfff$Beta2) / sqrt(2)
signed_dist = (dfff$Beta2 - dfff$Beta1) / sqrt(2)
dfff$Eucladian_Dist = eucladian_dist
dfff$Signed_dist = signed_dist
write.csv(dfff, "Median_Beta_weights.csv")

##### Plot the eucladian distances
# https://rkabacoff.github.io/datavis/Bivariate.html
# https://psyteachr.github.io/msc-data-skills/ggplot.html

library(ggbeeswarm)
library(scales)
library(Hmisc)
library(wesanderson)

ggplot(dfff,
       aes(x = factor(Intensities), 
           y = Eucladian_Dist, 
           label= ORN_labels,
           colour=Intensities)) +
  
  geom_point(show.legend=FALSE, size=1.5) +
  
  # geom_text(data=subset(dfff, Intensities == "14uW/mm2"),show.legend=FALSE, nudge_x = -0.30, size = 3) +
  
  # geom_line(aes(group=ORN_labels),show.legend=FALSE, size = .5) +

  stat_summary(position = position_jitterdodge(jitter.width = 1),
    fun.data = "mean_cl_boot",show.legend=FALSE,
    shape = "+",
    size = 1) +
  
  # geom_hline(yintercept=0, size=0.25) +
  
  ylab("Abs(Diagonal Distance)")+
  
  xlab("") +

  theme_classic() +
  
  theme(plot.margin = margin(2, 2, 2, 2, "cm"), aspect.ratio=1)+
  
  
  scale_colour_manual(values = wes_palette("Rushmore1", 3, type = "continuous"))

## Calculate summary stats for the eucladian distances

require(dplyr)

summary_dfff <- dfff %>% group_by(Intensities) %>% summarise(ave=mean(Eucladian_Dist),
                                             SD= sd(Eucladian_Dist), N = n(), 
                                             LB = ave - 1.96*(SD/sqrt(N)),
                                             UB = ave + 1.96*(SD/sqrt(N)))
write.csv(summary_dfff, "Summary_Abs_Eucladian_Distance.csv")

################ Plot the scatter plot for Beta weights
library(ggplot2)
library(ggthemes)
require("ggrepel")
par(mar=c(12,14,12,12)) 

ggplot(data = dfff, aes(x=Beta1, y=Beta2, col=Intensities, label = ORN_labels)) +
  
  geom_point(shape=19,
             alpha=1,
             size=3,
             stroke=0) +
  
  geom_text_repel(size = 3) +
  
  facet_wrap(facets = vars(Intensities)) +
  
  xlim(0,1.5) + ylim(0,1.5) +
  
  geom_abline(slope=1, intercept = 0) +
  
  theme_minimal() +
  
  theme(plot.margin = margin(2, 2, 2, 2, "cm"), aspect.ratio=1)

################# Exploration space#################################################
neuron1="Or42b"
neuron2="Or67d"
neuron12="Or67d-Or42b"

id1 = which(d14$Neuron == neuron1)
id2 = which(d14$Neuron == neuron2)
id12 = which(d14$Neuron == neuron12)

y1 = d14$weighted_TSALE_P10[id1]
y2 = d14$weighted_TSALE_P10[id2]
y12 = d14$weighted_TSALE_P10[id12]

y1 = sample(y1[order(y1)], 26, replace=TRUE)
y2 = sample(y2[order(y2)], 26, replace=TRUE)
y12 = sample(y12[order(y12)], 26, replace=TRUE)

y1_sq = y1**2
y2_sq = y2**2
y12_sq = y12**2

y1_sq = y1_sq[order(y1_sq)]
y2_sq = y2_sq[order(y2_sq)]
y12_sq = y12_sq[order(y12_sq)]
  
length(y12)
x<- y1_sq
y<- y2_sq
z <- outer(x, y, f)
z[is.na(z)] <- 1
op <- par(bg = "white")

f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }

persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")


wireframe(y12~y1*y2,colorkey=TRUE,drape=TRUE, zlim=c(0,2),scales=list(arrows=FALSE))

cloud(y12~y1*y2,colorkey=TRUE,drape=TRUE, zlim=c(-1,1),scales=list(arrows=FALSE))

densityplot(~weighted_TSALE_P10, groups=Neuron, data=d)

plot(lm(y12 ~ y1 + y2))
cone <- function(x, y){
  sqrt(x^2+y^2)
}
x <- y <- seq(-1, 1, length= 20)
z <- outer(x, y, cone)
persp(x, y, z)
persp(x, y, z,
      main="Perspective Plot of a Cone",
      zlab = "Height",
      theta = 30, phi = 15,
      col = "springgreen", shade = 0.5)


eth.lo <- loess(y12~y1*y2, span = 1/3, parametric = "C",
                drop.square = "C", family="symmetric")

eth.marginal <- list(y1 = seq(min(y1), max(y1), length.out = 25),
                     y2 = seq(min(y2), max(y2), length.out = 25))
eth.grid <- expand.grid(eth.marginal)
eth.fit <- predict(eth.lo, eth.grid)

wireframe(eth.fit ~ eth.grid$y1 * eth.grid$y2,
          shade=TRUE,scales=list(arrows=FALSE),
          screen = list(z = 45, x = -60, y=0),
          distance = .1,
          xlab = "ORN1", ylab = "ORN2", zlab = "combo")

##### Time series analyses #######################################################
data("AirPassengers")
AP <- AirPassengers

plot(AP, ylab="Passengers (1000s)", type="o", pch =20)

AP.decompM <- decompose(AP, type = "multiplicative")
plot(AP.decompM)









