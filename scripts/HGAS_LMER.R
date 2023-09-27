## HGAS_LMER
## Author(s): Brian Kielstra
## Originated: 2022-06-04
##
##
## Runs maximum likelihood inference regression models
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
library(lme4)
library(effects)

## *************
## ANALYSIS ----
## *************

## FUNCTIONS ----

## Some functions for LMER bootstrap processing

## Bootstrap
LMER_boot_est <- function(.) {
  c(beta=fixef(.), 
    as.data.frame(VarCorr(.))$sdcor[c(1:3,5)])
}

## Bootstrap initiate parallel processing
LMER_boot_initiate <- function(varlist, nclust){
  closeAllConnections()
  clust <- parallel::makeCluster(nclust)
  parallel::clusterEvalQ(clust, library("lme4"))
  parallel::clusterExport(cl = clust, varlist = varlist)
  showConnections()
  clust
}

## Summarize model coefficients
LMER_boot_est_summary <- function(merBoot){
  
  t <- as.data.frame(merBoot$t)
  
  names(t) <- c("Intercept", "Slope", "WATERBODY_CODE_SAMPLE_YEAR", 
                "WATERBODY_CODE", "WATERBODY_CODE_SLOPE", "RESIDUAL")
  
  est = apply(t, 2, function(x) as.numeric(quantile(x, probs = 0.5, na.rm = T)))
  upr = apply(t, 2, function(x) as.numeric(quantile(x, probs = 0.975, na.rm = T)))
  lwr = apply(t, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = T)))
  
  ret_tab <- as.data.frame(dplyr::bind_rows(est, lwr, upr))
  rownames(ret_tab) <- c("Estimate", "Lower", "Upper")  
  
  return(ret_tab)  
}

## Bootstrap predictions
LMER_boot_pred_summary <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

## Full model runs w/ bootMER, use.u = TRUE
LMER_mod_run <- function(dat, filter_statement, cores){
  
  message("Generate sub data")
  ## Generate subdata
  sdat <- dat %>% 
    dplyr::filter(eval(rlang::parse_expr(filter_statement)))
  
  message("Run model")
  ## Run model 
  mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + (1|WATERBODY_CODE_SAMPLE_YEAR), data = sdat)
  
  message("Evaluate model")
  ## Evaluate residuals
  hist(resid(mod), breaks = 100) # seems good 
  plot(resid(mod) ~ fitted(mod)) # seems reasonable; a few very high residuals 
  
  message("Generate predictions")
  ## Generate predictions for each waterbody based on the median size of LT in whole dataset 
  mod_predict <<- unique(sdat[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
  mod_predict$WEIGHT_GRAM_LOG <- log(1000)
  
  ## Bootstrap confidence intervals from bootMer
  clust <- LMER_boot_initiate(varlist = "mod_predict", nclust = cores)
  system.time(mod_boot_est <- lme4::bootMer(mod, LMER_boot_est, 
                                            nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                            parallel = "snow", 
                                            cl = clust, 
                                            ncpus = cores)) #17s for 100
  LMER_boot_est_summary(mod_boot_est) ## 889s
  
  ## Bootstrap predictions based on 2000 simulations
  mod_boot_pred <- function(., newdata) {
    predict(., newdata=mod_predict)
  }
  
  clust <- LMER_boot_initiate(varlist = "mod_predict", nclust = cores)
  system.time(mod_boot_pred_res <- lme4::bootMer(mod, mod_boot_pred, 
                                                 nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                 parallel = "snow", 
                                                 cl = clust, 
                                                 ncpus = cores)) # 150 
  
  (mod_boot_pred_res_sum <- LMER_boot_pred_summary(mod_boot_pred_res))
  
  ## Results 
  LMER_boot_est_summary(mod_boot_est)
  head(mod_boot_pred_res_sum)
  names(mod_boot_pred_res_sum) <- c("LMER_fit", "LMER_lwr", "LMER_upr")
  
  pred_frame <- cbind(mod_predict, mod_boot_pred_res_sum)
  
  message("Generate results")
  ret_list <- list(model = mod, pred_frame = pred_frame)
}

## Full model runs w/ bootMER, use.u = optional
LMER_mod_run <- function(dat, filter_statement, cores, use.u){
  
  message("Generate sub data")
  ## Generate subdata
  sdat <- dat %>% 
    dplyr::filter(eval(rlang::parse_expr(filter_statement)))
  
  message("Run model")
  ## Run model 
  mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + (1|WATERBODY_CODE_SAMPLE_YEAR), data = sdat)
  
  message("Evaluate model")
  ## Evaluate residuals
  hist(resid(mod), breaks = 100) # seems good 
  plot(resid(mod) ~ fitted(mod)) # seems reasonable; a few very high residuals 
  
  message("Generate predictions")
  ## Generate predictions for each waterbody based on the median size of LT in whole dataset 
  mod_predict <<- unique(sdat[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
  mod_predict$WEIGHT_GRAM_LOG <- log(1000)
  
  ## Generate confidence intervals 
  merTools::FEsim(mod)
  
  
  ## Bootstrap confidence intervals from bootMer
  clust <- LMER_boot_initiate(varlist = "mod_predict", nclust = cores)
  system.time(mod_boot_est <- lme4::bootMer(mod, LMER_boot_est, 
                                            nsim=2000, use.u = use.u, .progress = "txt",  ## 99 before
                                            parallel = "snow", 
                                            cl = clust, 
                                            ncpus = cores)) #17s for 100
  LMER_boot_est_summary(mod_boot_est) ## 889s
  
  ## Bootstrap predictions based on 2000 simulations
  mod_boot_pred <- function(., newdata) {
    predict(., newdata=mod_predict)
  }
  
  clust <- LMER_boot_initiate(varlist = "mod_predict", nclust = cores)
  system.time(mod_boot_pred_res <- lme4::bootMer(mod, mod_boot_pred, 
                                                 nsim=2000, use.u = use.u, .progress = "txt",  ## 99 before
                                                 parallel = "snow", 
                                                 cl = clust, 
                                                 ncpus = cores)) # 150 
  
  (mod_boot_pred_res_sum <- LMER_boot_pred_summary(mod_boot_pred_res))
  
  ## Results 
  LMER_boot_est_summary(mod_boot_est)
  head(mod_boot_pred_res_sum)
  names(mod_boot_pred_res_sum) <- c("LMER_fit", "LMER_lwr", "LMER_upr")
  
  pred_frame <- cbind(mod_predict, mod_boot_pred_res_sum)
  
  message("Generate results")
  ret_list <- list(model = mod, pred_frame = pred_frame)
}

## Full model runs, merTools predictInterval()
LMER_mod_merTools <- function(dat, filter_statement, cores = 1){
  
  message("Generate sub data")
  ## Generate subdata
  sdat <- dat %>% 
    dplyr::filter(eval(rlang::parse_expr(filter_statement)))
  
  message("Run model")
  ## Run model 
  mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + (1|WATERBODY_CODE_SAMPLE_YEAR), data = sdat)
  
  message("Evaluate model")
  ## Evaluate residuals
  hist(resid(mod), breaks = 100) # seems good 
  plot(resid(mod) ~ fitted(mod)) # seems reasonable; a few very high residuals 
  
  ## Coefficient confidence intervals 
  feEx <- merTools::FEsim(mod, 2000)
  
  message("Generate predictions")
  ## Generate predictions for each waterbody based on the median size of LT in whole dataset 
  mod_predict <<- unique(sdat[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
  mod_predict$WEIGHT_GRAM_LOG <- log(1000)
  
  pred <- merTools::predictInterval(mod, newdata = mod_predict, 
                          level = 0.95, 
                          n.sims = 2000, stat = "median", include.resid.var = TRUE, 
                          )
  names(pred) <- c("LMER_fit", "LMER_lwr", "LMER_upr")
  
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
Hg_LMER_LT_MASS <- LMER_mod_merTools(dat = Hg_dat_full, 
                                     filter_statement = "SPECIES_NAME == 'Lake Trout'")

Hg_LMER_NP_MASS <- LMER_mod_merTools(dat = Hg_dat_full, 
                                     filter_statement = "SPECIES_NAME == 'Northern Pike'")

Hg_LMER_WE_MASS <- LMER_mod_merTools(dat = Hg_dat_full, 
                                     filter_statement = "SPECIES_NAME == 'Walleye'")


Hg_LMER_results <- list(Hg_LMER_LT_MASS, Hg_LMER_NP_MASS, Hg_LMER_WE_MASS)

saveRDS(Hg_LMER_results, "./out_workspaces/HGAS_Hg_LMER_results.rds")

## ******* 
## AS ----
## ******* 

As_dat <- readRDS("./out_workspaces/HGAS_As_LTNPWEData.rds")
As_dat_full <- dplyr::bind_rows(As_dat)

## Model and diagnostics
As_LMER_LT_MASS <- LMER_mod_merTools(dat = As_dat_full, 
                                filter_statement = "Taxon == 'LT'")

## Model and diagnostics
As_LMER_NP_MASS <- LMER_mod_merTools(dat = As_dat_full, 
                                filter_statement = "Taxon == 'NP'")

## Model and diagnostics
As_LMER_WE_MASS <- LMER_mod_merTools(dat = As_dat_full, 
                                filter_statement = "Taxon == 'WALL'")

As_LMER_results <- list(As_LMER_LT_MASS, As_LMER_NP_MASS, As_LMER_WE_MASS)

saveRDS(As_LMER_results, "./out_workspaces/HGAS_As_LMER_results.rds")