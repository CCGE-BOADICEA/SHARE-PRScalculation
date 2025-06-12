library(data.table)

###############
##USER INPUTS##
###############
#multiple options; comment out the ones you don't intend to use

##option 1: hardcoded paths
propsFile  <- "Examples/MixedProportions/exA1.csv"           #Raw PGS and the ancestry proportions
parsFile   <- "PRSmodels_CanRisk/BCAC_309_PRS-UKB_all.prs"   #SNPs frequencies and logOR, and PGS model parameters
outputFile <- "Res-example_valuesA.csv"                      #Results of the calculation

##option 2: interactive use (e.g. Rstudio)
propsFile  <- readline(prompt="File for raw PGS and the ancestry proportions: ")
parsFile   <- readline(prompt="File for SNPs freq. and logOR, and PGS model pars: ")
outputFile <- readline(prompt="File for storing the results of the calculation: ")

##option 3: script called externally
args <- commandArgs(trailingOnly = TRUE)
propsFile  <- args[1]
parsFile   <- args[2]
outputFile <- args[3]


#################
##INPUT READING##
#################

#sample-specific parameters; raw PRS and population proportion, in one line
tmpTab  <- read.table(propsFile, header = TRUE, sep = ",", fill = TRUE)
rawPRS  <- tmpTab[1,1]
sampleProp <- as.numeric(tmpTab[1,2:5])

total_sum <- sum(sampleProp)
if (total_sum < 0.99 | total_sum > 1.01) {stop("Error in the input file!")}
sampleProp <- sampleProp / total_sum  # Scale values to sum exactly to 1

#SNPs frequencies and log(OR)
tmpTab  <- readLines(parsFile); tmpInd  <- 1
while(substr(tmpTab[tmpInd],1,1)=='#'){
  tmpInd <- tmpInd+1
}
tmpTab  <- tmpTab[-(1:(tmpInd-1))]
SNPtab  <- read.table(text=tmpTab, header = TRUE, skip= 4, sep = ",", fill = TRUE)
SNPtabR <- SNPtab[,c("eaf_eur","eaf_afr","eaf_eas","eaf_sas")]

#ancestry parameters (europeans,african,eastAsian,southAsian)
tmpTab     <- read.table(parsFile, header = FALSE, skip= 4, nrows=4, sep = ",", comment.char = "#", fill = TRUE)
monoThresh <- as.numeric(tmpTab[1,2:5])
POPmeans   <- as.numeric(tmpTab[2,2:5])
POPsigmas  <- as.numeric(tmpTab[3,2:5])
alphas     <- as.numeric(tmpTab[4,2:5])

monoFlag <- (sampleProp >= monoThresh)


################
##CALCULATIONS##
################

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

if ((outputFile=="") | (outputFile==" ") | (outputFile==".")){
  #Printing on screen
  print(sprintf("alpha= %f; PRS_Zscore= %f", alphaPar, Z_score))
} else {
  #Writing to CSV
  output_data <- data.frame(alpha = alphaPar, PRS_Zscore = Z_score)
  write.csv(output_data, outputFile, row.names = FALSE)
}