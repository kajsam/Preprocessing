#' ---
#' title: "Lung cancer lumi object included pooled samples"
#' author: "Kajsa Mollersen"
#' date: "02/03/2017"
#' output: pdf_document
#' ---

#' I found, in the R package 'nowac' a file named 'generate_all_datasets.R'. There I discovered that there
#' are 12 pooled samples that have been excluded when creating the lumi object. I will include them again.
#' This file is more or less copied code from 'generate_all_datasets.R'. I comment all changes, starting with 
#' kam025. I don't have the full overview of what's going on here. Maybe some day.

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=TRUE)
options(width=80)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60), comment="")
Sys.setlocale('LC_ALL','UTF-8') 
library(printr)

#' kam025: libraries
#+ libraries, include = FALSE
library(Biobase)

#' Creates lumi object from expression and
#' control data file. Removes user-defined bad saples, remove suffixes
#' and returns the object.
#' @export
#' @import lumi
#' @import Biobase
#' @import plyr
make_lumi <-
  function(input_path,
           exprs_filename,
           ctrl_data_filename,
           bad_samples = NULL,
           remove_suffix = NULL,
           output_path="~/") {
    lobj <-
      lumi::lumiR(
        paste0(input_path, exprs_filename),
        convertNuID = TRUE,
        lib.mapping = "lumiHumanIDMapping",
        inputAnnotation = FALSE,
        QC = FALSE
      )
    
    # control data name without file extension (csv or txt)
    ctrl_probe_name = unlist(strsplit(ctrl_data_filename, "\\."))[1]
    path = unlist(strsplit(input_path, "/"))
    dataset_name = path[length(path)]
    ctrl_probe_name = paste0(dataset_name, "_", ctrl_probe_name)
    
    # Fixes control data file. Some of there have an unnamed
    # column with row numbers crashing the addControlData
    # function below.
    control_filename = fix_troubled_file(paste0(input_path, ctrl_data_filename),
                                         path = output_path,
                                         sep = "\t",
                                         dataset_name = ctrl_probe_name)
    
    lobj <- lumi::addControlData2lumi(control_filename, lobj)
    
    if (!is.null(bad_samples)) {
      to_remove = Biobase::sampleNames(lobj) %in% bad_samples
      lobj = lobj[,-which(to_remove)]
    }
    
    if (!is.null(remove_suffix)) {
      bad_suffixes = paste(remove_suffix, collapse = "|")
      Biobase::sampleNames(lobj) = gsub(bad_suffixes, "", Biobase::sampleNames(lobj))
    }
    
    return(lobj)
  }

#' Creates a single RData file from expression file, control file, overview,
#' samplesheet etc.
generate_dataset_file <- function(input_path,
                                  exprs_filename,
                                  overview_filename,
                                  samplesheet_filename,
                                  ctrl_data_filename,
                                  bad_samples = NULL,
                                  remove_suffix = NULL,
                                  output_path,
                                  output_filename) {
  lobj <- make_lumi(input_path,
                    exprs_filename,
                    ctrl_data_filename,
                    bad_samples,
                    remove_suffix,
                    output_path)
  
  samplesheet.data <-
    read.table(
      paste(input_path, samplesheet_filename,
            sep = ""),
      sep = ",",
      header = TRUE,
      check.names = FALSE,
      skip = 7
    )
  
  samplesheet.data <- samplesheet.data[,-c(4, 5, 6)]
  overview.data <-
    read.table(
      paste(input_path, overview_filename,
            sep = ""),
      sep = ";",
      header = TRUE,
      check.names = FALSE
    )
  
  # only samples from lumi obj
  overview.data <- overview.data[overview.data$Sample_ID %in%
                                   sampleNames(lobj),]
  
  # select only negative controls
  CData.all <- getControlData(lobj)
  negativeCtrls <- CData.all[CData.all$controlType == "NEGATIVE",]
  
  output = paste0(output_path, output_filename)
  save(
    list = c(
      "lobj",
      "overview.data",
      "samplesheet.data",
      "negativeCtrls"
    ),
    file = output
  )
  
  cat(paste("Created dataset", output))
  
}

#' Check to see if the file contains numbers as row names. These files will
#' break the addControlData2lumi function call and must be fixed. We remove
#' all row names and return the filename of the new fixed file. If the file
#' is fine we don't do anything and return the original filename.
fix_troubled_file <- function(filename,
                              path = "~/",
                              sep = ";",
                              dataset_name="") {
  # Read file line by line, get number of columns and number of
  # data entries in the second line.
  file = readLines(filename)
  rows = strsplit(file, split = "\n")
  num_columns = length(unlist(strsplit(rows[[1]], split = sep)))
  num_data_entries = length(unlist(strsplit(rows[[2]], split = sep)))
  
  if (num_columns == num_data_entries) {
    cat("Control data file was ok.")
    return(filename)
  }
  
  # Pretty regex! What it does is that it substitutes the first
  # occurance of a number (any length) followed by a ';' on
  # every line in the file.  This way we're removing the row
  # names from the file.
  data_rows = rows[2:length(unlist(rows))]
  without_row_names = sub(paste0("[0-9]*", sep), "", data_rows)
  
  # Write contents to file and return filename
  new_filename = paste0(path, dataset_name, ".csv")
  
  cat(rows[[1]], sep = "\n", file = new_filename)
  cat(without_row_names,
      sep = "\n",
      file = new_filename,
      append = TRUE)
  
  cat(paste0("Wrote control data file to ", new_filename, "\n"))
  
  return(new_filename)
}

#' kam025: Instead of creating lumi objects from all datasets, I create only from the lung cancer data 
#' set. I have also changed 'input_path' and 'output_path'.   
#'
#' kam025: These are the "bad samples" that were removed in the earlier version: 
#' c("POOL1",  "POOL10", "POOL11", "POOL12", "POOL2",  "POOL3",
#'   "POOL4",  "POOL5"  ,"POOL6"  ,"POOL7",  "POOL8" , "POOL9")  
#'   
# Lung cancer data set 

#+ finnish, message = FALSE
generate_dataset_file(input_path="/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/ReadData/",
                      exprs_filename="SCannGCF373_Sample_Probe_Profile.txt",
                      overview_filename="OVERSIKT_RNA_cRNA_GCF373_2015_SendesTromso.csv", 
                      samplesheet_filename="GCF_2015_373_SampleSheet_Tromso.csv",
                      ctrl_data_filename="TableControl.txt",
                      remove_suffix=NULL,
                      bad_samples=NULL, 
                      output_path="/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/ReadData/",
                      output_filename="LungCancerPooled.RData")
