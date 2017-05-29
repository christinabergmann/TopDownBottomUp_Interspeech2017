#### Read in data and clean up ####

# This script was written to pre-process MetaLab compatible spreadsheets on behavioral word segmentation studies 
# Author: Christina Bergmann
# chbergma'at'gmail.com

#### libraries ####

library(tidyverse)

#### Data read in ####

#For now focusing on the eye tracking data
db = read.csv("data/InWordDB.csv")

#### Some cleanup ####

db = db %>%
  #filter(method == "HPP") %>%
  filter(test_lang == "native") %>%
  filter(datasubset == "Content Words") %>%
  filter(infant_type == "typical")