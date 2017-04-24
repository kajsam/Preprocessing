#' ---
#' title: "Prepocessing lung cancer data, Aberdeen style"
#' author: "Kajsa Mollersen"
#' date: "22/04/2017"
#' output: pdf_document
#' ---

#' The preprocessing done in 'preprocessing.R' led to some peculiar results. We made three different data sets, one keeping all 47'
#' probes, one excluding the Illumina labelled "Bad Probes" and "No match", and one excluding based on the wrongfully named "P 
#' value detection". Then Therese ranked the probes according to p-value (case-ctrl) with som limma function, presumably a 
#' straightforward t-test. The two first data sets gave fairly similar results, looking at the top 10 probes. Some of the top ranked
#' probes were excluded by the "P value detection", which is alarming. The Pvaldet doesn't make any sense statistically, so for the 
#' future I will just ignore it. 
#' 
#' I did have a look at the negative control probes, that supposedly measure the independent noise, as they are designed not to match
#' any human mRNA. In this context, "independent" means independent of signal, and of course independent of each other. Suspicious as 
#' I am, I normally don't trust any claims of independence, so I checked it. I used the 'background_noise.R' to extract the neg ctrl 
#' probe values, and then I switched to MatLab, where I wrote 'neg_noise_est.m' and 'cross_hybr.m'. The neg ctrl probe values clearly 
#' measure some other things than independent noise, namely cross-hybridization. Cross-hybridization is when mRNA clings to a probe it 
#' was not meant to cling to, and then it sends out a signal. It seems to be weak, but definetely not independent. 
#' 
#' I've decided on using only the non-hybridized (is that a word?) neg ctrl probes for independent noise estimation. First I excluded,
#' at random, 1/3 of the 770 neg ctrl probes, because I want to use the as quality control later on. Then I identified around 100 cross-
#' hybridized neg ctrl probes (details in 'neg_noise_est.m'), which will be excluded from the noise estimation. 
#' 
#' Hege B\o velstad used normexp transformation in her analysis, see Hege_preprocessing.Rmd, title: "lung-cancer-data-preprocessing", 
#' author: "Hege MB", date: "4 Feb 2016". I discussed preprocessing shortly with Claus-Dieter Mayer here in Aberdeen, and I agree with
#' when he says that it really shouldn't matter that much. So it should be safe to go with some standard procedure. 

#' There is an interesting paper, 'Optimizing the noise versus bias trade-off...', Shi et al (2010), that analyses pre-processing 
#' methods. One method they use is the same as Hege is using. Their conclusion is that there is no best pre-processing method, and they 
#' are all the same if you add the right offset. But there is no suggestion on how to find the right offset.
#' 
#' I will create datasets of three different sizes. (i) Full sized: no probes removed. (ii) Illumina size: only Bad probes and No match
#' (according to IlluminaHumanv4) are removed. (iii) DetectionP: Those probes with a "detection p-value" < threshold, and presens < 0.01 
#' are removed (in addition to Bad and No match). Hege set the threshold at 0.01. I will increase it to match her number of removed 
#' probes. 
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
idx_negCtrl <- 1:dim(negCtrl)[1]
datapath <- "/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/"

#' In 'neg_noise_est.m' and 'cross_hybr.m' I randomly chose a test set of neg ctrl probes for quality control, then I identified the
#' cross-hybridized probes from the rest. The remaining neg ctrl probes will be used for Pvaldet-based removal of regular probes and,
#' more importantly, normexp tranformation of probe values. 
#' 
#+ indepedent_negCtrlprobes
idx_testnegCtrl <- read.csv(paste(datapath,"probenr_test.csv", sep = ""), header = FALSE)
idx_testnegCtrl <- as.integer(as.matrix(idx_testnegCtrl))
idx_hybrnegCtrl <- read.csv(paste(datapath,"probenr_crosshybr.csv", sep = ""), header = FALSE)
idx_hybrnegCtrl <- as.integer(as.matrix(idx_hybrnegCtrl))
idx_indepnegCtrl <- setdiff(idx_negCtrl,union(idx_testnegCtrl,idx_hybrnegCtrl))

length(idx_indepnegCtrl) + length(idx_testnegCtrl) + length(idx_hybrnegCtrl) # Just checkin'

write_datapath <- "/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/Aberdeen/"
setwd(write_datapath)

#' This transformation includes background correction using negative controls, normexp transformation, and quantile 
#' normalisation (not using positive controls, why?). 
#+ Normexp_i
indnegCtrl <- negCtrl[idx_indepnegCtrl,]
totalData <- rbind(exprs(lobj),indnegCtrl) # use the independent neg controls for background correction
status <- c(rep("regular",nrow(exprs(lobj))), rep("negative", nrow(indnegCtrl))) # assign status
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

#' Remove those probes with detection p-value < 0.01 for less than 1% of the individuals. "The intensity associated with a certain 
#' detection p value corresponds to the respective upper quantile (percentile) of the (probability) distribution of the intensities 
#' of the negative controls." This quote is from the internet. I have written an e-mail to Illumina tech support to get the exact 
#' definition. Seem to be correct. 
#+ Normexp_iii
detection.pval <- 0.01  # Set the limits
present.limit <- 0.01
presentLim <- present.limit*ncol(ok.obj)
# At least this many ladies must have detection p-value above the limit for a probe to be "present".
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
