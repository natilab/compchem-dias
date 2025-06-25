## DIAS Processing
# Processing of Distortion-Interaction Activation Strain Model data.

# STEP 2: 
# Script to clean processed output

# author: Natalia Labadie
#--------------------------------------------------------------------


# packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

### Import processed results

# import all .csv files in loops and append to generate only one data.frame (tibble)
# add a column named "reaction" and assign a value that indicates which reaction the results correspond to
# Additional columns are also added to indicate the solvent and the level of theory used in the calculations.

alldata <- NULL

for (df in df_names) {
  for (lev in lev_theory) {
    
    # build path to results for corresp df and level
    path <- file.path('.', df_folders[df], lev, res_folder)
    # get .txt files from results directory
    files <-  list.files(path, pattern = "\\DIAS_proc.csv$")
    
    for (file in files) {
      
      # split filename to get data
      file_split <- str_split(file, "_")
      
      # get substrate name
      substrate <- file_split[[1]][1]
      # get reaction name
      reaction <- file_split[[1]][3]
      
      # split level of theory to get data
      lev_split <- str_split(lev, "-")
      
      # get functional
      func <- lev_split[[1]][1]
      # get basis set
      basis_set <- lev_split[[1]][2]
      # get solvent
      solvent <- lev_split[[1]][3]
      
      # import data
      data <- read_csv(file.path(path, file)) %>% 
        mutate(substrate = substrate,
               reaction = reaction,
               functional = func,
               basis = basis_set,
               solvent = solvent)
      
      alldata <- rbind(alldata, data)
    }
  }
}


### Clean data

# Get levels and set labels for factor variables

rxn_lev <- unique(alldata$reaction)
rxn_lab <- c("Product 5N", "Product 5X", "Product 6E", "Product 6Z",
             "Product 5N", "Product 5X", "Product 6E", "Product 6Z",
             "Product 3N", "Product 3X")

subs_lev <- unique(alldata$substrate)
subs_lab <- c("AllenylBPin", "AllenylCO2Me", "VinylBPin", "Methylacrylate")

func_lev <- unique(alldata$functional)
func_lab <- c("M062X")

bas_lev <- unique(alldata$basis)
bas_lab <- c("6-31+G*")

sv_lev <- unique(alldata$solvent)
sv_lab <- c("in vacuo", "toluene")


# transform dataset in order to:
# calculate mean distance
# turn qualitative vars into factors
# unite some columns in order to filter by condition

alldata <- alldata %>% 
  mutate(DistanceFinal = (Distance3 + Distance2)/2,
         reaction = factor(reaction,
                           levels = rxn_lev,
                           labels = rxn_lab),
         substrate = factor(substrate,
                            levels = subs_lev,
                            labels = subs_lab),
         functional = factor(functional,
                             levels = func_lev,
                             labels = func_lab),
         basis = factor(basis,
                        levels = bas_lev,
                        labels = bas_lab),
         solvent = factor(solvent,
                          levels = sv_lev,
                          labels = sv_lab)) %>% 
  unite(full_rxn, substrate, reaction, sep = "-", remove = FALSE) %>% 
  unite(level_theory, functional, basis, sep = "/", remove = FALSE) %>%   
  unite(full_condition, level_theory, solvent,
        sep = " ", remove = FALSE)   %>% 
  mutate(full_condition = ifelse(solvent == "toluene",
                                 str_replace(full_condition, 
                                             "toluene", 
                                             "in toluene"),
                                 full_condition))


### Get max energy 

# generate a df with the rows of the max energies, grouped by rxn and solvent
# with mutate add a column with the energies that is named Ets
maxs <- alldata %>% 
  group_by(full_rxn, full_condition) %>% 
  filter(Etotal == max(Etotal)) %>%
  mutate(Ets = Etotal)

# add the column Ets from mxs to alldata by joining both dfs
alldata <- left_join(alldata, maxs)

# add columns for the E dist and int in the TS
alldata <- mutate(alldata,
                  EtsD = ifelse(!is.na(Ets), EDistT, NA),
                  EtsI = ifelse(!is.na(Ets), EInt, NA),
                  EtsD1 = ifelse(!is.na(Ets), EDist1, NA),
                  EtsD2 = ifelse(!is.na(Ets), EDist2, NA))

###### Save clean data for future use

write_csv(alldata, "dias-2022_alldata_clean.csv") #generate csv with all processed data
saveRDS(alldata, file = "dias-2022_alldata_clean.rds") #save alldata object

###### Get energies at certain reaction points and save output

# Get energies at TS
data_ts_energies <- alldata %>% 
  filter(!is.na(Ets))

write_csv(data_ts_energies, "all_energies_ts.csv")

# Get energies at a mean distance of 2.3 Angstrom
data_23A <- alldata %>% 
  filter(round(DistanceFinal, 1) == 2.3)

write_csv(data_23A, "all_energies_23A.csv")
