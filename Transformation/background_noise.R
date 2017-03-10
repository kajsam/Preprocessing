#' ---
#' title: "Estimating the noise variance"
#' author: "Kajsa Mollersen"
#' date: "07/03/2017"
#' output: pdf_document
#' ---

#' I will play around with the negative controls and try to estimate the variance

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=TRUE)
options(width=80)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60), comment="")

library(printr)
library(genefilter)
library(limma)
library(lumi)
library(illuminaHumanv4.db)

#' This data is created by the 'create_lumi_object.R' file, and contains the 12 pooled samples. 
#' I have run nowaclean on them, but no obvious outliers were detected, hence none removed.
#+ readdata
datapath <- "/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/ReadData/"
load(paste(datapath,"LungCancerPooled.RData", sep = ""))
ls()
ctrlData <- controlData(lobj)
negCtrl <- ctrlData[ctrlData[,1]=="NEGATIVE",]
negCtrl <- negCtrl[,-(1:2)]
n_sam <- dim(negCtrl)[2]
n_obs <- dim(negCtrl)[1]

houseCtrl <- ctrlData[ctrlData[,1]=="HOUSEKEEPING",]
houseCtrl <- houseCtrl[,-(1:2)]
write.csv(houseCtrl, file="house_ctrls.csv") # so that I can do some of it in MatLab

strCtrl <- ctrlData[ctrlData[,1]=="LOW_STRINGENCY_HYB",]
strCtrl <- strCtrl[,-(1:2)]
write.csv(strCtrl, file="str_ctrls.csv") # so that I can do some of it in MatLab

probCtrl <- exprs(lobj)
write.csv(probCtrl, file="prob_ctrls.csv") # so that I can do some of it in MatLab


# write.csv(negCtrl, file="neg_ctrls.csv") # so that I can do some of it in MatLab

#' #' I suspect that the negative controls are not Gaussian distributed, as 'everyone' assumes without stating
#' #' it explicitly
#' #+ visual_inspection
#' n_ran <- 5
#' ran <- sample(1:n_sam, n_ran, replace = FALSE, prob = NULL)
#' ran
#' par(mfrow = c(1,n_ran))
#' for (j in 1:n_ran){
#'   hist(negCtrl[,ran[j]])
#'   }


