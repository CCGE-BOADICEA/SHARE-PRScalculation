library(data.table)

###############
##USER INPUTS##

#Ancestry proportions file contains the raw PGS and the ancestry proportions
propsFile <- "Examples/example_valuesA.csv"

#PGS info file contains SNPs frequencies and logOR, and PGS model parameters
parsFile  <- "./BCAC_307_PRS-UKB_all.prs"



#################
##INPUT READING##

#sample-specific parameters; raw PRS and population proportion, in one line
tmpTab  <- read.table(propsFile, header = TRUE, sep = ",", fill = TRUE)
rawPRS  <- tmpTab[1,1]
sampleProp <- as.numeric(tmpTab[1,2:5])
if(sum(sampleProp) != 1) stop("Error in the input file!")

#SNPs frequencies and log(OR)
SNPtab  <- read.table(parsFile, header = TRUE, skip= 8, sep = ",", fill = TRUE)
SNPtabR <- SNPtab[,c("eaf_eur","eaf_afr","eaf_eas","eaf_sas")]

#ethic parameters; europeans,african,eastAsian,southAsian
tmpTab     <- read.table(parsFile, header = FALSE, skip= 4, nrows=4, sep = ",", comment.char = "#", fill = TRUE)
monoThresh <- as.numeric(tmpTab[1,2:5])
POPmeans   <- as.numeric(tmpTab[2,2:5])
POPsigmas  <- as.numeric(tmpTab[3,2:5])
alphas     <- as.numeric(tmpTab[4,2:5])



################
##CALCULATIONS##

monoFlag <- (sampleProp >= monoThresh)
if (any(monoFlag)){
  #note: if the sample were truly mono-Ancestry, then the results from this if clause and those from the else clause would be identical
  monoIndex<- which(monoFlag)
  Z_score  <- (rawPRS-POPmeans[monoIndex])/POPsigmas[monoIndex]
  alphaPar <- alphas[monoIndex]
  
} else {
  
  ##calculate the possible overall sigma
  SNPtmp1 <- colSums((t(SNPtabR**2)*sampleProp))
  SNPtmp2 <- colSums((t(SNPtabR   )*sampleProp))
  
  stDev <- sum((POPsigmas**2)*sampleProp)
  stDev <- stDev + 4*sum((SNPtmp1 - (SNPtmp2**2)) * (SNPtab$Log_Odds_Ratio**2))
  stDev <- sqrt(stDev)
  
  ##calculate the weighted alphas and PRS
  tmp1 <- sum( alphas*sampleProp*(rawPRS-POPmeans)/POPsigmas )
  tmp2 <- sum( alphas*sampleProp/POPsigmas )
  tmp3 <- sum( alphas*sampleProp*POPsigmas )
  
  Z_score  <- (tmp1/tmp2)/stDev
  alphaPar <- tmp3/stDev
  
}

print(sprintf("alpha= %f; PRS_Zscore= %f", alphaPar, Z_score))
