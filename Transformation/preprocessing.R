#' ---
#' title: "Prepocessing lung cancer data included pooled samples"
#' author: "Kajsa Mollersen"
#' date: "02/03/2017"
#' output: pdf_document
#' ---

#' We're doing here several preprocessings. The preprocessing has two different things going on: (1) Transform the expression
#' from the raw value to something that is i.i.d. (2) Reduce the number of probes. Note that (1) and (2) are
#' independent of each other, so the execution order doesn't matter. 
#' 
#' (1) I will do this in two manners:   
#' (Normexp) According to the original code from Hege B\o velstad: Hege_preprocessing.Rmd, title: "lung-cancer-data-preprocessing", author: 
#' "Hege MB", date: "4 Feb 2016".  
#' (Kajsa) which I haven't decided on yet, but this is in case I desagree with Hege. 
#' There is an interesting paper, 'Optimizing the noise versus bias trade-off...', Shi et al (2010), that analyses pre-processing methods. 
#' One method they use is the same as Hege is using. Their conclusion is that there is no best pre-processing method, and they are all the 
#' same if you add the right offset. But there is no suggestion on how to find the right offset.
#' 
#' (2) We will create datasets of four different sizes. (i) Full sized: no probes removed. (ii) Illumina size: only Bad probes and No match
#' (according to IlluminaHumanv4) are removed. (iii) DetectionP: Those probes with a "detection p-value" < 0.01, and presens < 0.01 are 
#' removed (in addition to Bad and No match). (iv) Variance: Those probes with variance < noise variance are removed (in addition to
#' Bad and No match). I haven't figured out just yet how to do (iv).
#' 
#' The data set includes the 12 pooled samples. We keep them for quality control.

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=TRUE)
options(width=80)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60), comment="")
Sys.setlocale('LC_ALL','UTF-8') 
library(printr)

library(genefilter)
library(limma)
library(lumi)
library(illuminaHumanv4.db)

#' This data is created by the 'create_lumi_object.R' file, and contains the 12 pooled samples. I have run nowaclean on them, but no 
#' obvious outliers were detected, hence none removed.
#+ readdata
datapath <- "/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/ReadData/"
load(paste(datapath,"LungCancerPooled.RData", sep = ""))
ls()
ctrlData <- controlData(lobj)
negCtrl <- ctrlData[ctrlData[,1]=="NEGATIVE",]
negCtrl <- negCtrl[,-(1:2)]
dim(negCtrl)

write_datapath <- "/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/Transformation/"
setwd(write_datapath)

#' This transformation includes background correction using negative controls, normexp transformation, and quantile 
#' normalisation (not using positive controls, why?). 
#+ Normexp_i
totalData <- rbind(exprs(lobj),negCtrl) # use the neg controls for background correction
status <- c(rep("regular",nrow(exprs(lobj))), rep("negative", nrow(negCtrl))) # assign status
totalData.nec <- nec(totalData, status)
exprs.nec <- totalData.nec[which(status=="regular"),] # exclude the negative controls
nec.obj <- lobj
# make a lumi object where only the expression values have changed
exprs(nec.obj) <- exprs.nec 
# quantile normalisation of the backgr.corr. + normexp exprs values
norm.obj <- lumiN(nec.obj,method= "quantile") 
dim(norm.obj)

#' We plot the three versions of the expression values. I am not sure what I am looking for
#+ plotNormexp
par(mfrow = c(1,3))
Qlobj <- lumiQ(lobj)
plot(Qlobj ,what = "density", legend = FALSE, col = 1, lwd = .05, xlim = c(4,12), 
     main = "Raw data")
Qnobj <- lumiQ(nec.obj)
plot(Qnobj ,what = "density", legend = FALSE, col = 1, lwd = .05, xlim = c(4,12), 
     main = "Backgr. corr. + normexp")
normQ <- lumiQ(norm.obj)
plot(normQ ,what = "density", legend = FALSE, col = 1, lwd = .05, xlim = c(4,12), 
     main = "Normalised")

#' Saving Normexp(i)  
#+ saveNormexp_Fullsize
save(norm.obj, file="Normexp_Full.RData")

#' This is Normexp(ii). We are now removing 'Bad probes' and 'No match' as defined by illuminaHumanv4. 
#+ Normexp_ii
probes <- nuID2IlluminaID(as.character(featureNames(lobj)), lib.mapping=NULL, 
                          species = "Human",idType = 'Probe')
probe.quality <- unlist(mget(as.character(probes), illuminaHumanv4PROBEQUALITY, 
                             ifnotfound = NA))
which(is.na(probe.quality)) # just checking for NA's
length(which(probe.quality=="Bad")) # Lots of Bad quality probes
length(which(probe.quality=="No match")) # Just a few without a match
ok.quality <- !((probe.quality == "Bad")|(probe.quality == "No match"))
# Extracting the ok quality probes from the pre-processed lumi object
ok.obj <- norm.obj[ok.quality,] 
dim(ok.obj) # Checking the dimension

#' Saving Normexp(ii):  
#+ saveNormexp_Illimunasize
save(ok.obj, file="Normexp_Illumina.RData")

#' Remove those probes with detection p-value < 0.01 for less than 1% of the individuals. "The intensity associated with a certain detection p value 
#' corresponds to the respective upper quantile (percentile) of the (probability) distribution of the intensities 
#' of the negative controls." This quote is from the internet. I have written an e-mail to Illumina tech support to get 
#' the exact definition.
#+ Normexp_iii
detection.pval <- 0.01  # Set the limits
present.limit <- 0.01
presentLim <- present.limit*ncol(ok.obj)
# At least this many ladies must have detection p-value above the limit for a probe to be 
# "present".
presentLim   
present_count <- detectionCall(ok.obj, Th=detection.pval, type ='probe') 
present <- (present_count > presentLim)
present.obj <- norm.obj[present,]
dim(present.obj)

#' Saving Normexp(iii):  
#+ saveNormexp_Presentsize
save(present.obj, file="Normexp_Present.RData")

#' Have added a check after Therese discovered that I had managed to save the same lumi object to two different files.
rm(list = ls())
load("Normexp_Full.RData")
load("Normexp_Illumina.RData")
load("Normexp_Present.RData")
ls()
dim(norm.obj)
dim(ok.obj)
dim(present.obj)
