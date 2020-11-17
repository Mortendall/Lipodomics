library(tidyverse)
library(lipidr)
library(ggplot2)
library(limma)
library(vroom)
library(openxlsx)
library(here)
library(fs)



#' Lipodomics data loader
#'
#' @param file_type the name of the count file
#'
#' @return a trimmed count matrix

count_file_loader <- function(file_type){
  count_file <- fs::dir_ls(here("data_raw/"),
                           regexp = file_type,
                           recurse = TRUE)
  count_data <- read.xlsx(count_file)
  rownames(count_data) <- count_data$`mol%`
  count_data <- count_data %>% dplyr::select(-`mol%`)
  return(count_data)
}

#' Load metadata
#'
#' @param file_name the name of the metadata file
#'
#' @return metadata file

load_metadata <- function(file_name) {
  data_file <- fs::dir_ls(here::here("data_raw/"),
                          regexp = file_name,
                          recurse = T)
  metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
  return(metadata)
}
