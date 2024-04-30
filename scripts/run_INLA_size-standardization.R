
#### Running INLA prediction model with the MECP mercury database ####
#' This script is provided so researchers wishing to apply INLA on their dataset can produce predicted concentrations 
#' for a set fish size using INLA models and the MECP mercury database

### Step 1: Import packages ----

## Install the INLA package, if it has not already been done
## install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

## Load libraries 
library(readxl)
library(lme4)
library(INLA)
library(effects)
library(tidyverse)
library(ggpubr)

### Step 2: Import the mercury database ----
#' you can replace this database with any alternative you wish, or combine it with other datasets as long as the other data includes these fields:
#' - WATERBODY_CODE: A unique identifier for each lake/waterbody in the dataset
#' - SPECIES_NAME: the common name of the fish OR a species code (ensure you are using the same ones in your data)
#' - WEIGHT_GRAM_LOG: the log-transformed weight of each fish (You could also use length, but you will need to validate the predictions using a test-train approach)
#' - VALUE_LOG: the log-transformed mercury concentration of the fish
#' - SAMPLE_YEAR: the year of the sample collection


hg <- readxl::read_excel("./data/Ontario Inland 3 Species Hg 2008-2020-12-16.xlsx") %>% 
  ## filter data to only fish from lakes, produced by MECP, Walleye, Laketrout and Northern Pike and only ones that are from dorsal filet samples
  filter(SPECIES_NAME %in% c("Lake Trout", "Walleye", "Northern Pike"), 
         PORTION_TYPE_DESC %in% c("SKINLESS, BONELESS FILLET (STANDARD MOE DORSAL FILLET)", 
                                  "SKINLESS, BONELESS FILLET PLUG (SKIN-OFF)")) %>% 
  mutate(WEIGHT_GRAM_LOG = log(WEIGHT_GRAM), LENGTH_CM_LOG = log(LENGTH_CM), 
         VALUE_LOG = log(VALUE)) %>% 
  dplyr::select(WATERBODY_CODE, SAMPLE_YEAR, SPECIES_NAME, CONTAMINANT, VALUE, VALUE_LOG, LENGTH_CM, LENGTH_CM_LOG, WEIGHT_GRAM, WEIGHT_GRAM_LOG) %>% 
  filter(!is.na(LENGTH_CM_LOG), 
         !is.na(WEIGHT_GRAM_LOG))

### Step 3: Import your data ----
#' Make sure your data includes these fields:
#' - WATERBODY_CODE: A unique identifier for each lake/waterbody in the dataset
#' - SPECIES_NAME: the common name of the fish OR a species code (ensure you are using the same ones in the database)
#' - WEIGHT_GRAM_LOG: the log-transformed weight of each fish (You could also use length, but you will need to validate the predictions using a test-train approach)
#' - VALUE_LOG: the log-transformed mercury concentration of the fish
#' - SAMPLE_YEAR: the year of the sample collection

data <- readxl::()

data <- hg[sample(500, 1:nrow(hg)), ]

### Step 4: Select the sizes you want to predict concentrations for (in g) ----

sizes = c(500, 1000)

### Step 5: Create prediction data ----

# Find the unique event combinations for the analysis

unique.sets <- data %>% 
  dplyr::select(WATERBODY_CODE, SAMPLE_YEAR, SPECIES_NAME) %>% 
  distinct()

pred.data <- lapply(sizes, function(x){temp <- unique.sets; temp$WEIGHT_GRAM_LOG <- log(x); temp$VALUE_LOG <- NA; return(temp)}) %>% bind_rows()

pred.data <- pred.data %>% rbind(data) %>% rbind(hg %>% dplyr::select(WATERBODY_CODE, SAMPLE_YEAR, SPECIES_NAME, WEIGHT_GRAM_LOG, VALUE_LOG)) %>%
  mutate(EVENT = paste0(WATERBODY_CODE, SAMPLE_YEAR)) ## Note: you should account for the year effect where possible, but even if you only have one year for each of your lakes, the EVENT will account for the year-specific random effects in the model training.


### Step 6: Run the INLA model ----

INLA_MASS <-inla(VALUE_LOG ~WEIGHT_GRAM_LOG+ 
       #' The next two lines code for the random slopes effect due to waterbody, - basically codes for the random effect of the 
       #' two variables, expecting covariance of these variables  
       #' See the documentation for the description of this implementation 'inla.doc("iid2d")'
       f(WATERBODY_CODE1, n = 2*unique(x$n_waterbody), model = "iid2d") + 
       f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
       f(EVENT, model = "iid"),
     data = pred.data, 
     control.predictor = list(
       compute = TRUE, 
       quantiles = c(0.025, 0.5, 0.975) ### You can include more output quantiles if you wish to expand your predictions, or plot broader distributions of probable concentrations
     ),
     control.compute = list(
       cpo = TRUE
     )
)

### Step 7: extract predictions ----

## We extract the median, and the percentiles of the distribution of probable predictions.
#' Note: in our paper, we used the "median" to represent a single predicted value - something that is often required for downstream modelling, and what we validated this tool for
#' but you could also export more quantiles and use simulated distributions in downstream modelling. We have given you a three quantile example so you can choose which direction you want to take your analysis.

# Extract the median predicted value
pred.data$INLA_posterior_q50 = INLA_MASS$summary.fitted.values[, "0.5quant"]
# Extract the 2.5th percentile 
pred.data$INLA_posterior_q2p5 = INLA_MASS$summary.fitted.values[, "0.025quant"]
# Extract the 97.5th percentile
pred.data$INLA_posterior_q97p5 = INLA_MASS$summary.fitted.values[, "0.975quant"]
