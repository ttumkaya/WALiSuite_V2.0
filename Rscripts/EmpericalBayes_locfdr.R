## locdfr emprical bayes method

library(locfdr)

z_scores <- read.csv('C:/Users/tumkayat/Desktop/deltadelta_ZScores.csv')
Zi <- z_scores$Zi

#data(lfdrsim)

#zex <- lfdrsim[, 2]

w <- locfdr(Zi)
str(w)
w$fdr
which(w$fdr<.20)

## Hyungwon's EBprot library
### Tutorial ###
data(UPS1)

### removing contaminants and reversed peptide results ###
rid = UPS1$Reverse == "+" | UPS1$Contaminant == "+"
UPS1 = UPS1[!rid, ]

#### creating the ratios for analysis ####
UPS1$R300vs100fmol = log(1/UPS1$Ratio.H.L.normalized.300fmol)-log(1/UPS1$Ratio.H.L.normalized.100fmol)

#### extracting relevant information from dataset ####
UPS1 = UPS1[,c(1,29,30,98)]

## renaming the dataset in the required manner ##
names(UPS1) = c("Peptide.seq","Protein.id", "Gene.name", "R300vs100fmol")

### remove proteins with missing ratios ###
id = which(is.na(UPS1[,4]))
UPS1 =UPS1[-id,]

str(UPS1)

### independent design ###
par(mfrow=c(1,2))
output.pepLvl <- EBprot(UPS1, c("R300vs100fmol"), nullK=3, verbose=TRUE, design=1, min.rep = 1,
                        min.k=5, BFDR.threshold=0.05, fold.require=c(1/2,2), outlier.step=TRUE, report.sig=FALSE,
                        refit = FALSE, recalibrate = FALSE, plotBox = FALSE,plot.range = c(-2,2))

### My DATA ###
library(EBprot)
allData <- read.csv('C:/Users/tumkayat/Desktop/deltadelta_ZScores_Air_Effect_Only_Used.csv')
#allData <- allData[allData$ORNs != 'Gr66a',]
#allData_70uw <- allData[allData$LightInt == "70uW",]

str(allData)
allData_filtered <- allData[,c(2,3,4)]
#allData_70uw_filtered <- allData[allData$LightInt == "70uW",][,c(2,3,4)]

#names(allData_70uw_filtered) = c("Peptide.seq","Protein.id","Zi")
names(allData_filtered) = c("Peptide.seq","Protein.id","Zi")

########## Separtely
data1 <- read.csv('C:/Users/tumkayat/Desktop/male_fed_NoAir.csv')
data2 <- read.csv('C:/Users/tumkayat/Desktop/male_fed_Air.csv')
data3 <- read.csv('C:/Users/tumkayat/Desktop/male_starved_NoAir.csv')
data4 <- read.csv('C:/Users/tumkayat/Desktop/male_starved_Air.csv')

str(data2)
data_filtered1 <- data1[,c(3,2,4)]
data_filtered2 <- data2[,c(3,2,4)]
data_filtered3 <- data3[,c(3,2,4)]
data_filtered4 <- data4[,c(3,2,4)]

data_filtered1 <- data_filtered1[data_filtered1$ORNs != 'Gr66a',]
data_filtered2 <- data_filtered2[data_filtered2$ORNs != 'Gr66a',]
data_filtered3 <- data_filtered3[data_filtered3$ORNs != 'Gr66a',]
data_filtered4 <- data_filtered4[data_filtered4$ORNs != 'Gr66a',]

names(data_filtered1) = c("Peptide.seq","Protein.id","ES")
which(output.pepLvl$ES.BFDR<.3)

total1<-merge(data_filtered1, data_filtered2,by=c("Peptide.seq","Protein.id"))
total2<-merge(data_filtered3, data_filtered4,by=c("Peptide.seq","Protein.id"))
grand_total<-merge(total1, total2, by=c("Peptide.seq","Protein.id"))

str(grand_total)
names(grand_total) = c("Peptide.seq","Protein.id","ES1","ES2","ES3","ES4")
#######

#c("ES.x.x","ES.y.x","ES.x.y","ES.y.y")
par(mfrow=c(1,2))
output.pepLvl <- EBprot(allData_filtered, c("Zi"), nullK=3, verbose=TRUE, design=1, min.rep = 1,
                        min.k=5, BFDR.threshold=0.05, fold.require=c(1/2,2), outlier.step=TRUE, report.sig=FALSE,
                        refit = FALSE, recalibrate = FALSE, plotBox = FALSE,plot.range = c(-2,2),null.mean.zero=TRUE)

which(output.pepLvl$ES.Probability > 0.5) 
output.pepLvl[output.pepLvl$ZI.BFDR < 0.3,]
dim(output.pepLvl[output.pepLvl$ES.X.X.Probability > 0.5 | output.pepLvl$ES.X.X.Probability < -0.5,])
dim(allData_filtered) 

write.csv(output.pepLvl, file = "C:/Users/tumkayat/Desktop/deltadelta_All_EBprot.csv")

### FOR XYZ ##############
data <- read.csv('C:/Users/tumkayat/Desktop/TSAR_ORNdata.csv')
str(data_filtered)
data_filtered <- data[,c(2,1,3)] 
names(data_filtered) = c("Peptide.seq","Protein.id","ES")
dim(data_filtered)

output.pepLvl <- EBprot(data_filtered, c("ES"), nullK=3, verbose=TRUE, design=1, min.rep = 1,
                        min.k=5, BFDR.threshold=0.05, fold.require=c(1/2,2), outlier.step=TRUE, report.sig=FALSE,
                        refit = FALSE, recalibrate = FALSE, plotBox = FALSE,plot.range = c(-2,2),null.mean.zero=FALSE)

dim(output.pepLvl[output.pepLvl$ES.BFDR < 0.3,])
output.pepLvl[output.pepLvl$ES.Probability > 0.5 | output.pepLvl$ES.Probability < -0.5,]

write.csv(output.pepLvl, file = "C:/Users/tumkayat/Desktop/TSAR_RESULTS.csv")

