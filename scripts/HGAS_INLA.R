## HGAS_INLA
## Author(s): Brian Kielstra
## Originated: 2022-06-04
##
##
## Runs approximate Bayesian inference regression models
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
library(dplyr)
library(INLA)

## *************
## ANALYSIS ----
## ************* 

## FUNCTIONS ---- 

## Structure derived from: 
## https://rstudio-pubs-static.s3.amazonaws.com/93118_c8d057fc0f754e7dafbd61c17ac5708c.html

INLA_mod <- function(dat, filter_statement){
  
  sdat <- dat %>% 
    dplyr::filter(eval(rlang::parse_expr(filter_statement)))
  
  sdat$WATERBODY_CODE <- factor(sdat$WATERBODY_CODE) # for inla model
  sdat$WATERBODY_CODE1 <- as.integer(sdat$WATERBODY_CODE) # for inla model
  sdat$WATERBODY_CODE2 <- sdat$WATERBODY_CODE1 + max(sdat$WATERBODY_CODE1) # for inla model
  n_waterbody <- dplyr::n_distinct(sdat$WATERBODY_CODE) # for inla model
  
  ## Run model
  mod <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                f(WATERBODY_CODE1, n = 2 * n_waterbody, model = "iid2d") +   
                f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
              data = sdat, 
              control.predictor = list(
                compute = TRUE, 
                quantiles = c(0.025, 0.5, 0.975)
              ),
              control.compute = list(
                cpo = TRUE
              )
  )
  summary(mod)
  
  ## Now predictions 
  
  ## Predictions in INLA are made by appending your prediction information to the dataframe used to run the model 
  ## Except, the VALUE_LOG values are set to NA
  ## Then you extract the posterior information 
  
  mod_predict <<- unique(sdat[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")])
  mod_predict$WEIGHT_GRAM_LOG <- log(1000)
  mod_predict$VALUE_LOG <- NA
  
  sdat_prediction <- bind_rows(sdat, mod_predict)
  nrow(sdat)
  nrow(sdat_prediction)
  
  mod_prediction <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                           f(WATERBODY_CODE1, n = 2 * sdat$n_waterbody, model = "iid2d") +   
                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                         data = sdat_prediction, 
                         control.predictor = list(
                           compute = TRUE, 
                           quantiles = c(0.025, 0.5, 0.975)
                         ),
                         control.compute = list(
                           cpo = TRUE
                         )
  )
  
  mod_prediction_posteriors <- mod_prediction$summary.fitted.values[
    (nrow(sdat) + 1):nrow(sdat_prediction),
    c("0.025quant", "0.5quant", "0.975quant")
  ]
  
  nrow(mod_prediction_posteriors)
  
  mod_predict <- cbind(mod_predict, mod_prediction_posteriors)
  mod_predict <- mod_predict[!duplicated(mod_predict),]
  
  message("Generate results")
  ret_list <- list(model = mod_prediction, pred_frame = mod_predict)
  
}

## *******
## HG ----
## ******* 

Hg_dat <- readRDS("./out_workspaces/HGAS_Hg_LTNPWEData.rds")
Hg_dat_full <- dplyr::bind_rows(Hg_dat)

## Model and diagnostics
Hg_INLA_LT_MASS <- INLA_mod(dat = Hg_dat_full, 
                            filter_statement = "SPECIES_NAME == 'Lake Trout'")

Hg_INLA_NP_MASS <- INLA_mod(dat = Hg_dat_full, 
                            filter_statement = "SPECIES_NAME == 'Northern Pike'")

Hg_INLA_WE_MASS <- INLA_mod(dat = Hg_dat_full, 
                            filter_statement = "SPECIES_NAME == 'Walleye'")


Hg_INLA_results <- list(Hg_INLA_LT_MASS, Hg_INLA_NP_MASS, Hg_INLA_WE_MASS)

saveRDS(Hg_INLA_results, "./out_workspaces/HGAS_Hg_INLA_results.rds")

## ******* 
## AS ----
## ******* 

As_dat <- readRDS("./out_workspaces/HGAS_As_LTNPWEData.rds")
As_dat_full <- dplyr::bind_rows(As_dat)

## Model and diagnostics
As_INLA_LT_MASS <- INLA_mod(dat = As_dat_full, 
                            filter_statement = "Taxon == 'LT'")

## Model and diagnostics
As_INLA_NP_MASS <- INLA_mod(dat = As_dat_full, 
                            filter_statement = "Taxon == 'NP'")

## Model and diagnostics
As_INLA_WE_MASS <- INLA_mod(dat = As_dat_full, 
                            filter_statement = "Taxon == 'WALL'")

As_INLA_results <- list(As_INLA_LT_MASS, As_INLA_NP_MASS, As_INLA_WE_MASS)

saveRDS(As_INLA_results, "./out_workspaces/HGAS_As_INLA_results.rds")