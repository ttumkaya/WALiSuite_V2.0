d = read.csv("AllData/Or67d_weighted_TSALE_values.csv", header=T, as.is=T, check.names=FALSE)
cid = c("Fly ID","Light Intensity(uW/mm2)","Wind status", "Satiety", "Sex", "Genotype","Status","weighted_TSALE_P10")
d = d[,cid]
write.table(d, "AllData_cleaned/Or67d_wTSALE_reduced.txt", sep="\t", row.names=F)


d = read.csv("AllData/Or85d_weighted_TSALE_values.csv", header=T, as.is=T, check.names=FALSE)
cid = c("Fly ID","Light Intensity(uW/mm2)","Wind status", "Satiety", "Sex", "Genotype","Status","weighted_TSALE_P10")
d = d[,cid]
write.table(d, "AllData_cleaned/Or85d_wTSALE_reduced.txt", sep="\t", row.names=F)


d = read.csv("AllData/Or67d_Or85d_weighted_TSALE_values.csv", header=T, as.is=T, check.names=FALSE)
cid = c("Fly ID","Light Intensity(uW/mm2)","Wind status", "Satiety", "Sex", "Genotype","Status","weighted_TSALE_P10")
d = d[,cid]
write.table(d, "AllData_cleaned/Or67d-Or85d_wTSALE_reduced.txt", sep="\t", row.names=F)


