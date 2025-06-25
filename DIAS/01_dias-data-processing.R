## DIAS Processing
# Processing of Distortion-Interaction Activation Strain Model data.

# STEP 1: 
# Script to process output from autoDIAS python script

# author: Natalia Labadie
#--------------------------------------------------------------------

# packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)


# Global variables

hartreetokcal <- 627.5

variable_names <- c("RxnStep", "Distance1", "Distance2", "Distance3", "Angle", 
                    "Dihedral", "TotalE", "EInt", "EDistT", "EDist1", "EDist2", 
                    "scfC", "scf1", "scf2")

# Data processing function

procdata <- function(path, filename_in, filename_out, 
                     scfbasal1, scfbasal2) {
  
  # imports multiple autoDIAS ouput (.txt) and saves processed data in one .csv file
  dataproc <- read_delim(file.path(path, filename_in), 
                         delim = "    ",
                         col_names = variable_names,
                         col_types = "ccccccccccccccccccccc") %>%
    slice(-1:-4) %>%
    select_if(function(x){!all(is.na(x))}) %>%
    mutate_if(is.character, as.numeric) %>%
    select(matches("Step"), 
           starts_with("Distance"), 
           starts_with("scf")) %>%
    mutate(Etotal = scfC - (scfbasal1 + scfbasal2), 
           EDist1 = scf1 - scfbasal1, 
           EDist2 = scf2 - scfbasal2, 
           EDistT = EDist1 + EDist2, 
           EInt = Etotal - EDistT)  %>%
    mutate_at(vars(starts_with("E")), 
              function(x){x = x*hartreetokcal}) %>%
    write_csv(paste0(path, '/', filename_out, ".csv"), col_names = TRUE)
  
  return(dataproc)
  
}


###### Processing of multiple reaction outputs (for phd Thesis)


# Strings for building paths to different reaction results

df_names <- c("abpin", "aco2me", "vbpin", "meacr") # names of compounds

# folder names with files
df_folders <- c("allenylbpin-DIAS", "allenylco2me-DIAS",
                "vinylbpin-DIAS", "meacrylate-DIAS") 

names(df_folders) <- df_names

# levels of theory used
lev_theory <- c("m062x-631+Gd-gp", "m062x-631+Gd-tol")

# folder name for output
res_folder <- "results"


# SCF values for optimized basal geometries of compounds 
# each vector has values for both levels of theory, ordered accordingly

scf_list <- list(scf_cp = c(-194.011887749, -194.019261168),
                 scf_abpin = c(-527.133738146, -527.143118375),
                 scf_aco2me = c(-344.395985534, -344.403495522),
                 scf_vbpin = c(-489.074323486, -489.081974101),
                 scf_meacr = c(-306.342572142, -306.348247844))

# set names to vector values according to level of theory
scf_list <- lapply(scf_list, function(x)
  setNames(x, lev_theory))


## Process DIAS data

# create log file to track processing
log_file <- file("dias_R_proc.log", open='a')

# in this loop, process DIAS results for multiple reactions (multiple txt files)
# clean data 
# and save processed data in .csv files
# the loop saves 1 csv for each combination of reaction and level of theory
for (df in df_names) {
  for (lev in lev_theory) {
    
    # build path to results for corresp df and level
    path <- file.path('.', df_folders[df], lev, res_folder)
    # get .txt files from results directory
    files <-  list.files(path, pattern = "\\DIAS.txt$")
    
    for (file in files) {
      # get reaction name
      reaction <- str_extract(file, pattern = "ts[^_]*(?=_)")
      # create output filename
      filename_out <- paste(df, lev, reaction, 'DIAS_proc', 
                            sep = '_')
      
      # get basal scf energy values for the corresp iteration
      scfbas1 <- unname(scf_list$scf_cp[lev])
      
      scf2_nam <- paste0('scf_', df)
      scfbas2 <- unname(scf_list[[scf2_nam]][lev])
      
      # process with proc_data
      cat(paste("Processing", file, "in", path), file = log_file,
          append = TRUE, sep = '\n')
      
      procdata(path, file, filename_out, 
               scfbasal1 = scfbas1, scfbasal2 = scfbas2)
      
      cat(paste("Created output", paste0(filename_out, '.csv'), 
                "in", path),
          file = log_file, append = TRUE, sep = '\n')
      
    }
  }
}


