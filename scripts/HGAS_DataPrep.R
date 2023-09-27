## HgAS_DataPrep
## Author(s): Brian Kielstra
## Originated: 2022-06-04
##
##
## Prepares Hg and As data for analysis
##
##
## Script information:
## Major comments with two hashtags (procedural)
## Minor comments with one hashtag (functional)
## Major code chunks delineated by:
##    ##***********
##    ## NAME"----" (delete quotes to be navigable in R Studio Document outline)
##    ## **********
## Code lines should not exceed 80 characters

## **********
## SETUP ----
## ********** 

## Load libraries 
library(readxl)
library(lme4)
## install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
library(effects)
library(sjPlot)
library(rstanarm)

## *******
## Hg ----
## *******

## Read the saved file
Hg <- readxl::read_excel("./data/Ontario Inland 3 Species Hg 2008-2020-12-16.xlsx")
nrow(Hg) #39319

sum( (Hg$WEIGHT_GRAM > 900 & Hg$WEIGHT_GRAM < 1100) , na.rm = TRUE) # 3843 fish

## Log transform
Hg$LENGTH_CM_LOG <- log(Hg$LENGTH_CM) 
Hg$WEIGHT_GRAM_LOG <- log(Hg$WEIGHT_GRAM)
Hg$VALUE_LOG <- log(Hg$VALUE)

Hg$WATERBODY_CODE_SAMPLE_YEAR <- paste(Hg$WATERBODY_CODE, Hg$SAMPLE_YEAR, sep = "_")

## Subset by "SKINLESS, BONELESS FILLET" OR "SKINLESS, BONELESS FILLET PLUG"
Hg_sub <- subset(Hg, PORTION_TYPE_DESC %in% 
                   c("SKINLESS, BONELESS FILLET (STANDARD MOE DORSAL FILLET)", 
                     "SKINLESS, BONELESS FILLET PLUG (SKIN-OFF)"))

## Subset and remove those lengths and weights that are NA
Hg_sub <- subset(Hg_sub, !is.na(LENGTH_CM_LOG))
Hg_sub <- subset(Hg_sub, !is.na(WEIGHT_GRAM_LOG))

length(unique(Hg_sub$LOCATION_NAME)) # 1197 lakes
length(unique(Hg_sub$WATERBODY_CODE_SAMPLE_YEAR)) # 1915 sampling events

## 1.x) Data exploration for VALUE (Hg), WEIGHT_GRAM, and LENGTH_CM ----
hist(Hg_sub$VALUE); summary(Hg_sub$VALUE)
dotchart(Hg_sub$VALUE) # majority of data < 1 but some outliers 
dotchart(log(Hg_sub$VALUE)) # few outliers and majority of data b/w -1 and 0 as expected 
hist(log(Hg_sub$VALUE))

hist(Hg_sub$WEIGHT_GRAM); summary(Hg_sub$WEIGHT_GRAM)
dotchart(Hg_sub$WEIGHT_GRAM) # majority of data < 2000 but some outliers 
dotchart(log(Hg_sub$WEIGHT_GRAM)) # few outliers and majority of data b/w -1 and 0 as expected 
hist(log(Hg_sub$WEIGHT_GRAM))

hist(Hg_sub$LENGTH_CM); summary(Hg_sub$LENGTH_CM)
dotchart(Hg_sub$LENGTH_CM) # majority of data < 2000 but some outliers 
dotchart(log(Hg_sub$LENGTH_CM)) # few outliers and majority of data b/w -1 and 0 as expected 
hist(log(Hg_sub$LENGTH_CM))

## 1.x) All species w/ dataset-based median length, mass ----

median_LENGTH_CM <- sapply(unique(Hg_sub$SPECIES_CODE), function(x) {
  Hg_sub_sub <- subset(Hg_sub, SPECIES_CODE == x)
  res <- median(Hg_sub_sub$LENGTH_CM, na.rm = T)
  return(res)
}, USE.NAMES = TRUE)

hist(Hg_sub$LENGTH_CM, breaks = 100); abline(v = median_LENGTH_CM, lty = 2, col = "red")

median_WEIGHT_GRAM <- sapply(unique(Hg_sub$SPECIES_CODE), function(x) {
  Hg_sub_sub <- subset(Hg_sub, SPECIES_CODE == x)
  res <- median(Hg_sub_sub$WEIGHT_GRAM, na.rm = TRUE)
  return(res)
}, USE.NAMES = TRUE)

## Get median size for all species and construct prediction dataframe
Hg_median_frame <- data.frame(
  SPECIES_CODE = unique(Hg_sub$SPECIES_CODE),
  LENGTH_CM = median_LENGTH_CM,
  WEIGHT_GRAM = median_WEIGHT_GRAM,
  LENGTH_CM_LOG = log(median_LENGTH_CM), 
  WEIGHT_GRAM_LOG = log(median_WEIGHT_GRAM), 
  stringsAsFactors = FALSE
)

Hg_LT <- subset(Hg_sub, SPECIES_NAME == "Lake Trout")
Hg_NP <- subset(Hg_sub, SPECIES_NAME == "Northern Pike")
Hg_WE <- subset(Hg_sub, SPECIES_NAME == "Walleye")

## Some sampling statistics - how many lakes have more than one sampling event? 
nrow(Hg_LT); length(with(Hg_LT, unique(WATERBODY_CODE))); 
median(table(Hg_LT$WATERBODY_CODE_SAMPLE_YEAR))
min(table(Hg_LT$WATERBODY_CODE_SAMPLE_YEAR))
max(table(Hg_LT$WATERBODY_CODE_SAMPLE_YEAR))

summary(Hg_LT$WEIGHT_GRAM); quantile(Hg_LT$WEIGHT_GRAM, probs = c(0.025, 0.975))

nrow(Hg_NP); length(with(Hg_NP, unique(WATERBODY_CODE))); 
median(table(Hg_NP$WATERBODY_CODE_SAMPLE_YEAR))
min(table(Hg_NP$WATERBODY_CODE_SAMPLE_YEAR))
max(table(Hg_NP$WATERBODY_CODE_SAMPLE_YEAR))

summary(Hg_NP$WEIGHT_GRAM); quantile(Hg_NP$WEIGHT_GRAM, probs = c(0.025, 0.975))

nrow(Hg_WE); length(with(Hg_WE, unique(WATERBODY_CODE))); 
median(table(Hg_WE$WATERBODY_CODE_SAMPLE_YEAR))
min(table(Hg_WE$WATERBODY_CODE_SAMPLE_YEAR))
max(table(Hg_WE$WATERBODY_CODE_SAMPLE_YEAR))

summary(Hg_WE$WEIGHT_GRAM); quantile(Hg_WE$WEIGHT_GRAM, probs = c(0.025, 0.975))

## *******
## AS ----
## *******

## Read the saved file
As <- read.csv("./data/Fish_As_2021.12.01.csv")
As_sub <- subset(As, System_Type == "Lake")
As_sub <- subset(As_sub, Data_source == "MECP")

## Subset by "SKINLESS, BONELESS FILLET" OR "SKINLESS, BONELESS FILLET PLUG"
As_sub <- subset(As_sub, PORTION_TYPE_DESC %in% 
                   c("SKINLESS, BONELESS FILLET (STANDARD MOE DORSAL FILLET)", 
                     "SKINLESS, BONELESS FILLET PLUG (SKIN-OFF)"))

nrow(unique(As_sub[,c("Waterbody", "System_Type")])) #125

## Some data cleanup and variable creation to align with Hg analysis
As_sub$VALUE <- As_sub$As_ug_ww
As_sub$LENGTH_CM <- As_sub$TLEN/10
As_sub$WEIGHT_GRAM <- As_sub$RWT
As_sub$SAMPLE_YEAR <- as.factor(as.numeric(sapply(strsplit(As_sub$Sampling_Date, split = "-"), function(x){x[[1]]})))
As_sub$WEIGHT_GRAM_LOG <- log(As_sub$WEIGHT_GRAM)
As_sub$LENGTH_CM_LOG <- log(As_sub$LENGTH_CM)
As_sub$VALUE_LOG <- log(As_sub$VALUE)
As_sub$WATERBODY_CODE <- As_sub$Waterbody
As_sub$SPECIES_CODE <- As_sub$Taxon 
As_sub$WATERBODY_CODE_SAMPLE_YEAR <- paste(As_sub$WATERBODY_CODE, As_sub$SAMPLE_YEAR, sep = "_")

median_LENGTH_CM <- sapply(unique(As_sub$SPECIES_CODE), function(x) {
  As_sub_sub <- subset(As_sub, SPECIES_CODE == x)
  res <- median(As_sub_sub$LENGTH_CM, na.rm = T)
  return(res)
}, USE.NAMES = TRUE)

hist(As_sub$LENGTH_CM, breaks = 100); abline(v = median_LENGTH_CM, lty = 2, col = "red")

median_WEIGHT_GRAM <- sapply(unique(As_sub$SPECIES_CODE), function(x) {
  As_sub_sub <- subset(As_sub, SPECIES_CODE == x)
  res <- median(As_sub_sub$WEIGHT_GRAM, na.rm = TRUE)
  return(res)
}, USE.NAMES = TRUE)

hist(As_sub$WEIGHT_GRAM, breaks = 100); abline(v = median_WEIGHT_GRAM, lty = 2, col = "red")

## Get median size for all species and construct prediction dataframe
As_median_frame <- data.frame(
  SPECIES_CODE = unique(As_sub$SPECIES_CODE),
  LENGTH_CM = median_LENGTH_CM,
  WEIGHT_GRAM = median_WEIGHT_GRAM,
  LENGTH_CM_LOG = log(median_LENGTH_CM), 
  WEIGHT_GRAM_LOG = log(median_WEIGHT_GRAM), 
  stringsAsFactors = FALSE
)

## 2.1) Data exploration for VALUE (Hg), WEIGHT_GRAM, and LENGTH_CM ----
hist(As_sub$VALUE); summary(As_sub$VALUE)
dotchart(As_sub$VALUE)
dotchart(log(As_sub$VALUE)) 
hist(log(As_sub$VALUE))

hist(As_sub$WEIGHT_GRAM); summary(As_sub$WEIGHT_GRAM)
dotchart(As_sub$WEIGHT_GRAM)  
dotchart(log(As_sub$WEIGHT_GRAM)) 
hist(log(As_sub$WEIGHT_GRAM))

hist(As_sub$LENGTH_CM); summary(As_sub$LENGTH_CM)
dotchart(As_sub$LENGTH_CM) 
dotchart(log(As_sub$LENGTH_CM)) 
hist(log(As_sub$LENGTH_CM))

## 2.2) Subsetting by species ----
As_LT <- subset(As_sub, Taxon == "LT")
As_NP <- subset(As_sub, Taxon == "NP")
As_WE <- subset(As_sub, Taxon == "WALL") 

## Some sampling statistics - how many lakes have more than one sampling event? 
nrow(As_LT); length(with(As_LT, unique(WATERBODY_CODE))); 
median(table(As_LT$WATERBODY_CODE_SAMPLE_YEAR))
min(table(As_LT$WATERBODY_CODE_SAMPLE_YEAR))
max(table(As_LT$WATERBODY_CODE_SAMPLE_YEAR))

summary(As_LT$WEIGHT_GRAM); quantile(As_LT$WEIGHT_GRAM, probs = c(0.025, 0.975))

nrow(As_NP); length(with(As_NP, unique(WATERBODY_CODE))); 
median(table(As_NP$WATERBODY_CODE_SAMPLE_YEAR))
min(table(As_NP$WATERBODY_CODE_SAMPLE_YEAR))
max(table(As_NP$WATERBODY_CODE_SAMPLE_YEAR))

summary(As_LT$WEIGHT_GRAM); quantile(As_LT$WEIGHT_GRAM, probs = c(0.025, 0.975))

nrow(As_WE); length(with(As_WE, unique(WATERBODY_CODE))); 
median(table(As_WE$WATERBODY_CODE_SAMPLE_YEAR))
min(table(As_WE$WATERBODY_CODE_SAMPLE_YEAR))
max(table(As_WE$WATERBODY_CODE_SAMPLE_YEAR))

summary(As_LT$WEIGHT_GRAM); quantile(As_LT$WEIGHT_GRAM, probs = c(0.025, 0.975))

## ***********
## EXPORT ----
## *********** 

saveRDS(list(Hg_LT, Hg_NP, Hg_WE), "./out_workspaces/HGAS_Hg_LTNPWEData.rds")
saveRDS(list(As_LT, As_NP, As_WE), "./out_workspaces/HGAS_As_LTNPWEData.rds")
