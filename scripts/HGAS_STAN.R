## HGAS_STAN
## Author(s): Brian Kielstra
## Originated: 2022-06-04
##
##
## Runs Bayesian inference regression models
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

library(dplyr)
library(rstanarm)

## *************
## ANALYSIS ----
## *************

## FUNCTIONS ----

STAN_mod <- function(dat, filter_statement, mod_save){
  
  sdat <- dat %>% 
    dplyr::filter(eval(rlang::parse_expr(filter_statement)))
  
  system.time(mod <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                              (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                                              (1|WATERBODY_CODE_SAMPLE_YEAR), data = sdat, 
                                            cores = 4, chains = 4, iter = 2000))
  saveRDS(mod, mod_save)
  
  message("Generate predictions")
  ## Generate predictions for each waterbody based on the median size of LT in whole dataset 
  mod_predict <<- unique(sdat[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
  mod_predict$WEIGHT_GRAM_LOG <- log(1000)
  
  pred <- rstanarm::posterior_predict(mod,
                                      newdata = mod_predict)
  
  pred <- apply(pred,
                MARGIN = 2,
                function(x) {
                  quantile(x, probs = c(0.5, 0.025, 0.975))
                }
  )
  
  pred <- (t(pred))
  colnames(pred) <- c("STAN_fit", "STAN_lwr", "STAN_upr")
  pred_frame <- cbind(mod_predict, pred)

  message("Generate results")
  ret_list <- list(model = mod, pred_frame = pred_frame)
}
  
## *******
## HG ----
## ******* 

Hg_dat <- readRDS("./out_workspaces/HGAS_Hg_LTNPWEData.rds")
Hg_dat_full <- dplyr::bind_rows(Hg_dat)

## Model and diagnostics
Hg_STAN_LT_MASS <- STAN_mod(dat = Hg_dat_full, 
                            filter_statement = "SPECIES_NAME == 'Lake Trout'", 
                            mod_save = "./out_workspaces/Hg_STAN_LT_MASS.rds")

Hg_STAN_NP_MASS <- STAN_mod(dat = Hg_dat_full, 
                            filter_statement = "SPECIES_NAME == 'Northern Pike'",
                            mod_save = "./out_workspaces/Hg_STAN_NP_MASS.rds")

Hg_STAN_WE_MASS <- STAN_mod(dat = Hg_dat_full, 
                            filter_statement = "SPECIES_NAME == 'Walleye'",
                            mod_save = "./out_workspaces/Hg_STAN_WE_MASS.rds")

Hg_STAN_results <- list(Hg_STAN_LT_MASS, Hg_STAN_NP_MASS, Hg_STAN_WE_MASS)

saveRDS(Hg_STAN_results, "./out_workspaces/HGAS_Hg_STAN_results.rds")

## ******* 
## AS ----
## ******* 

As_dat <- readRDS("./out_workspaces/HGAS_As_LTNPWEData.rds")
As_dat_full <- dplyr::bind_rows(As_dat)

## Model and diagnostics
As_STAN_LT_MASS <- STAN_mod(dat = As_dat_full, 
                            filter_statement = "Taxon == 'LT'",
                            mod_save = "./out_workspaces/As_STAN_LT_MASS.rds")

As_STAN_LT_MASS$model

## Model and diagnostics
As_STAN_NP_MASS <- STAN_mod(dat = As_dat_full, 
                            filter_statement = "Taxon == 'NP'",
                            mod_save = "./out_workspaces/As_STAN_NP_MASS.rds")

## Model and diagnostics
As_STAN_WE_MASS <- STAN_mod(dat = As_dat_full, 
                            filter_statement = "Taxon == 'WALL'",
                            mod_save = "./out_workspaces/As_STAN_WE_MASS.rds")

As_STAN_results <- list(As_STAN_LT_MASS, As_STAN_NP_MASS, As_STAN_WE_MASS)

saveRDS(As_STAN_results, "./out_workspaces/HGAS_As_STAN_results.rds") 