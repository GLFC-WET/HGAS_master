## HGAS_ModelDiagnosticsSummaries
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

library(dplyr)
library(ggplot2)
library(patchwork)
library(influence.ME)
library(rstanarm)
library(ggplot2)

## ********
## SER ----
## ******** 

Hg_SER <- readRDS("./out_workspaces/HGAS_Hg_dat_full_SER_results.rds")
As_SER <- readRDS("./out_workspaces/HGAS_As_dat_full_SER_results.rds")

## RMSE and R2, event
Hg_SER_summary <- Hg_SER %>% 
  group_by(SPECIES_CODE) %>% 
  mutate(RMSE = as.numeric(RMSE), 
         r2 = as.numeric(r2), 
         slp = as.numeric(slp), 
         int = as.numeric(int)) %>%
  summarise(RMSE_mean = mean(RMSE),
            RMSE_2p5 = quantile(RMSE, probs = 0.025),
            RMSE_97p5 = quantile(RMSE, probs = 0.975),
            R2_mean = mean(r2), 
            R2_2p5 = quantile(r2, probs = 0.025),
            R2_97p5 = quantile(r2, probs = 0.975), 
            sd_slope = sd(slp, na.rm = T), 
            sd_int = sd(int, na.rm = T))

As_SER_summary <- As_SER %>% 
  group_by(SPECIES_CODE) %>% 
  mutate(RMSE = as.numeric(RMSE), 
         r2 = as.numeric(r2), 
         slp = as.numeric(slp), 
         int = as.numeric(int)) %>%
summarise(RMSE_mean = mean(RMSE),
          RMSE_2p5 = quantile(RMSE, probs = 0.025),
          RMSE_97p5 = quantile(RMSE, probs = 0.975),
          R2_mean = mean(r2), 
          R2_2p5 = quantile(r2, probs = 0.025),
          R2_97p5 = quantile(r2, probs = 0.975), 
          sd_slope = sd(slp, na.rm = T), 
          sd_int = sd(int, na.rm = T))

## Some investigative plots
Hg_wb <- as.numeric(Hg_SER$int)
Hg_slope <- as.numeric(Hg_SER$slp) 

p1 <- ggplot() + geom_histogram(aes(Hg_wb), bins = 100) + 
  theme_classic() + xlab("WATERBODY")
p2 <- ggplot() + geom_histogram(aes(Hg_slope), bins = 100) + 
  xlab("SLOPE") + theme_classic()

(p1 + p2)

As_wb <- as.numeric(As_SER$int)
As_slope <- as.numeric(As_SER$slp) 

p3 <- ggplot() + geom_histogram(aes(As_wb), bins = 100) + 
  theme_classic() + xlab("WATERBODY")
p4 <- ggplot() + geom_histogram(aes(As_slope), bins = 100) + 
  xlab("SLOPE") + theme_classic()

(p3 + p4)


## **************
## ML - LMER ----
## ************** 

Hg_ML <- readRDS("./out_workspaces/HGAS_Hg_LMER_results.rds")
As_ML <- readRDS("./out_workspaces/HGAS_As_LMER_results.rds")

names(Hg_ML) <- c("LT", "NP", "WE")
names(As_ML) <- c("LT", "NP", "WE")

## Some investigative plots
lapply(Hg_ML, function(xx){
  
  mod_resid <- (resid(xx$model, type = "pearson"))
  mod_fitted <- fitted(xx$model)
  ranef_wb <- (ranef(xx$model)$WATERBODY_CODE$`(Intercept)`)
  ranef_wbsy <- (ranef(xx$model)$WATERBODY_CODE_SAMPLE_YEAR$`(Intercept)`)
  ranef_slope <- (ranef(xx$model)$WATERBODY_CODE$WEIGHT_GRAM_LOG) 
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})

lapply(As_ML, function(xx){
  
  mod_resid <- (resid(xx$model, type = "pearson"))
  mod_fitted <- fitted(xx$model)
  ranef_wb <- (ranef(xx$model)$WATERBODY_CODE$`(Intercept)`)
  ranef_wbsy <- (ranef(xx$model)$WATERBODY_CODE_SAMPLE_YEAR$`(Intercept)`)
  ranef_slope <- (ranef(xx$model)$WATERBODY_CODE$WEIGHT_GRAM_LOG) 
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})

## Models look reasonable to proceed

## **************
## AB - INLA ----
## **************

Hg_B <- readRDS("./out_workspaces/HGAS_Hg_INLA_results.rds")
As_B <- readRDS("./out_workspaces/HGAS_As_INLA_results.rds")

names(Hg_B) <- c("LT", "NP", "WE")
names(As_B) <- c("LT", "NP", "WE")

lapply(Hg_B, function(xx){
  
  summary(xx$model)
  
  mod_frame <- xx$model$.args$data
  
  mod_fitted <- xx$model$summary.fitted.values$`0.5quant`
  mod_resid <- mod_frame$VALUE_LOG - mod_fitted
  
  ranef_wb <- xx$model$summary.random$WATERBODY_CODE1[min(mod_frame$WATERBODY_CODE1):max(mod_frame$WATERBODY_CODE1),]$`0.5quant`
  ranef_slope <- xx$model$summary.random$WATERBODY_CODE2[min(mod_frame$WATERBODY_CODE2):max(mod_frame$WATERBODY_CODE2),]$`0.5quant`
  ranef_wbsy <- xx$model$summary.random$WATERBODY_CODE_SAMPLE_YEAR$`0.5quant`
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})

lapply(As_B, function(xx){
  
  summary(xx$model)
  
  mod_frame <- xx$model$.args$data
  
  mod_fitted <- xx$model$summary.fitted.values$`0.5quant`
  mod_resid <- mod_frame$VALUE_LOG - mod_fitted
  
  ranef_wb <- xx$model$summary.random$WATERBODY_CODE1[min(mod_frame$WATERBODY_CODE1):max(mod_frame$WATERBODY_CODE1),]$`0.5quant`
  ranef_slope <- xx$model$summary.random$WATERBODY_CODE2[min(mod_frame$WATERBODY_CODE2):max(mod_frame$WATERBODY_CODE2),]$`0.5quant`
  ranef_wbsy <- xx$model$summary.random$WATERBODY_CODE_SAMPLE_YEAR$`0.5quant`
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})

## Models look reasonable to proceed

## *****************
## B - RSTANARM ----
## *****************

## https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html

Hg_B <- readRDS("./out_workspaces/HGAS_Hg_STAN_results.rds")
As_B <- readRDS("./out_workspaces/HGAS_As_STAN_results.rds")

names(Hg_B) <- c("LT", "NP", "WE")
names(As_B) <- c("LT", "NP", "WE")

lapply(Hg_B, function(xx){
  
  summary(xx$model)
  
  mod_frame <- xx$model$data
  mod_fitted <- xx$model$fitted.values
  mod_resid <- xx$model$residuals
  
  sims <- as.matrix(xx$model)
  
  a_quant <- apply(
    X = sims,
    MARGIN = 2,
    FUN = quantile,
    probs = c(0.025, 0.50, 0.975)
  )
  a_quant <- data.frame(t(a_quant))
  names(a_quant) <- c("Q2_5", "Q50", "Q97_5")
  
  a_quant$param <- rownames(a_quant)
  
  ranef_wb <- a_quant[grep("b\\[\\(Intercept\\) WATERBODY_CODE\\:", a_quant$param),]$Q50
  ranef_slope <- a_quant[grep("b\\[WEIGHT_GRAM_LOG WATERBODY_CODE:", a_quant$param),]$Q50
  ranef_wbsy <- a_quant[grep("b\\[\\(Intercept\\) WATERBODY_CODE_SAMPLE_YEAR\\:", a_quant$param),]$Q50
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})

lapply(As_B, function(xx){
  
  summary(xx$model)
  
  mod_frame <- xx$model$data
  mod_fitted <- xx$model$fitted.values
  mod_resid <- xx$model$residuals
  
  sims <- as.matrix(xx$model)
  
  a_quant <- apply(
    X = sims,
    MARGIN = 2,
    FUN = quantile,
    probs = c(0.025, 0.50, 0.975)
  )
  a_quant <- data.frame(t(a_quant))
  names(a_quant) <- c("Q2_5", "Q50", "Q97_5")
  
  a_quant$param <- rownames(a_quant)
  
  ranef_wb <- a_quant[grep("b\\[\\(Intercept\\) WATERBODY_CODE\\:", a_quant$param),]$Q50
  ranef_slope <- a_quant[grep("b\\[WEIGHT_GRAM_LOG WATERBODY_CODE:", a_quant$param),]$Q50
  ranef_wbsy <- a_quant[grep("b\\[\\(Intercept\\) WATERBODY_CODE_SAMPLE_YEAR\\:", a_quant$param),]$Q50
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})

## Models look reasonable to proceed

## *******************************************************************
## Calculate RMSE (event), RMSE (global), R2 (event), R2 (global) ----
## ******************************************************************* 

## Event is at sampling event level 
## Global is at whole analysis level

### SER ----- 
Hg_SER_summary
As_SER_summary

## Generate fitted values, calculate R2 at lake level 

### ML - LMER ---- 

Hg_ML_mod_summ <- lapply(Hg_ML, function(xx){
  
  mod_frame <- xx$model@frame
  mod_frame$fitted <- fitted(xx$model)
  mod_frame$resid <- resid(xx$model)
  
  res <- lapply(unique(mod_frame$WATERBODY_CODE_SAMPLE_YEAR), function(yy){
    
    mod_frame_sub <- subset(mod_frame, WATERBODY_CODE_SAMPLE_YEAR == yy) 
    
    R2 <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame_sub))$r.squared
    (RMSE <- sqrt(mean(mod_frame_sub$resid^2)))
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, 
               R2, RMSE)
    
  })
  
  res_bind <- bind_rows(res)
  res_bind
  
  ## 
  (R2_global <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame))$r.squared)
  (RMSE_global <-sqrt(mean(mod_frame$resid^2)))
  
  list(res_bind, data.frame(R2_global, RMSE_global))
  
 })

As_ML_mod_summ <- lapply(As_ML, function(xx){
  
  mod_frame <- xx$model@frame
  mod_frame$fitted <- fitted(xx$model)
  mod_frame$resid <- resid(xx$model)
  
  res <- lapply(unique(mod_frame$WATERBODY_CODE_SAMPLE_YEAR), function(yy){
    
    mod_frame_sub <- subset(mod_frame, WATERBODY_CODE_SAMPLE_YEAR == yy) 
    
    R2 <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame_sub))$r.squared
    (RMSE <- sqrt(mean(mod_frame_sub$resid^2)))
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, 
               R2, RMSE)
    
  })
  
  res_bind <- bind_rows(res)
  res_bind
  
  ## 
  (R2_global <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame))$r.squared)
  (RMSE_global <-sqrt(mean(mod_frame$resid^2)))
  
  list(res_bind, data.frame(R2_global, RMSE_global))
  
})

### AB - INLA ----

Hg_B_mod_summ <- lapply(Hg_B, function(xx){
  
  mod_frame <- xx$model$.args$data
  mod_frame$fitted <- xx$model$summary.fitted.values$`0.5quant`
  mod_frame$resid <- mod_frame$VALUE_LOG - mod_frame$fitted
  
  res <- lapply(unique(mod_frame$WATERBODY_CODE_SAMPLE_YEAR), function(yy){
    
    mod_frame_sub <- subset(mod_frame, WATERBODY_CODE_SAMPLE_YEAR == yy) 
    
    (R2 <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame_sub))$r.squared)
    (RMSE <- sqrt(mean(mod_frame_sub$resid^2, na.rm = T)))
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, 
               R2, RMSE)
    
  })
  
  res_bind <- bind_rows(res)
  res_bind
  
  ## 
  (R2_global <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame))$r.squared)
  (RMSE_global <-sqrt(mean(mod_frame$resid^2, na.rm = T)))
  
  list(res_bind, data.frame(R2_global, RMSE_global))
  
})

As_B_mod_summ <- lapply(As_B, function(xx){
  
  mod_frame <- xx$model$.args$data
  mod_frame$fitted <- xx$model$summary.fitted.values$`0.5quant`
  mod_frame$resid <- mod_frame$VALUE_LOG - mod_frame$fitted
  
  res <- lapply(unique(mod_frame$WATERBODY_CODE_SAMPLE_YEAR), function(yy){
    
    mod_frame_sub <- subset(mod_frame, WATERBODY_CODE_SAMPLE_YEAR == yy) 
    
    (R2 <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame_sub))$r.squared)
    (RMSE <- sqrt(mean(mod_frame_sub$resid^2, na.rm = T)))
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, 
               R2, RMSE)
    
  })
  
  res_bind <- bind_rows(res)
  res_bind
  
  ## 
  (R2_global <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame))$r.squared)
  (RMSE_global <-sqrt(mean(mod_frame$resid^2, na.rm = T)))
  
  list(res_bind, data.frame(R2_global, RMSE_global))
  
})

### B - RSTANARM ---- 

Hg_B_mod_summ <- lapply(Hg_B, function(xx){
  
  mod_frame <- xx$model$data
  mod_frame$fitted <- xx$model$fitted.values
  mod_frame$resid <- xx$model$residuals
  
  res <- lapply(unique(mod_frame$WATERBODY_CODE_SAMPLE_YEAR), function(yy){
    
    mod_frame_sub <- subset(mod_frame, WATERBODY_CODE_SAMPLE_YEAR == yy) 
    
    (R2 <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame_sub))$r.squared)
    (RMSE <- sqrt(mean(mod_frame_sub$resid^2, na.rm = T)))
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, 
               R2, RMSE)
    
  })
  
  res_bind <- bind_rows(res)
  res_bind
  
  ## 
  (R2_global <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame))$r.squared)
  (RMSE_global <-sqrt(mean(mod_frame$resid^2, na.rm = T)))
  
  list(res_bind, data.frame(R2_global, RMSE_global))
  
})

As_B_mod_summ <- lapply(As_B, function(xx){
  
  mod_frame <- xx$model$data
  mod_frame$fitted <- xx$model$fitted.values
  mod_frame$resid <- xx$model$residuals
  
  res <- lapply(unique(mod_frame$WATERBODY_CODE_SAMPLE_YEAR), function(yy){
    
    mod_frame_sub <- subset(mod_frame, WATERBODY_CODE_SAMPLE_YEAR == yy) 
    
    (R2 <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame_sub))$r.squared)
    (RMSE <- sqrt(mean(mod_frame_sub$resid^2, na.rm = T)))
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, 
               R2, RMSE)
    
  })
  
  res_bind <- bind_rows(res)
  res_bind
  
  ## 
  (R2_global <- summary(lm(VALUE_LOG ~ fitted, data = mod_frame))$r.squared)
  (RMSE_global <-sqrt(mean(mod_frame$resid^2, na.rm = T)))
  
  list(res_bind, data.frame(R2_global, RMSE_global))
  
})

### Generate summary table ----

## STOPPED HERE

## *******************************************************************
## Calculate RMSE and R2 for predictions for 1000 g fish +/- 50 g ---- 
## ******************************************************************* 

## Generate prediction test by species 

Hg_SER$WATERBODY_CODE_SAMPLE_YEAR <- paste(Hg_SER$WATERBODY_CODE, 
                                           Hg_SER$SAMPLE_YEAR, 
                                           sep = "_")

## Species of interest
spc <- c("LT", "NP", "WE")

## Run series of analysis to collect observations and predictions for fish 
## >950 g and <1050 g per species 
xx = "LT"
obspred_list <- lapply(spc, function(xx){
  
  print(xx)
  ## Get species dataframe 
  mod_frame <- Hg_AB[[xx]]$model$.args$data
  
  ## Find which observations meet size criteria
  ind <- which(mod_frame$WEIGHT_GRAM > 950 & mod_frame$WEIGHT_GRAM < 1050)
  mod_sub <- mod_frame[ind,]
  
  ## Generate species code
  if(xx == "LT"){
    spc_code <- "081"
  } else if (xx == "NP"){
    spc_code <- "131"
  } else {
    spc_code <- "334"
  }
  
  ## generate SER predictions meeting criteria
  SER_sub <- subset(Hg_SER, (WATERBODY_CODE_SAMPLE_YEAR %in% mod_sub$WATERBODY_CODE_SAMPLE_YEAR) &
                      SPECIES_CODE == spc_code )
  
  mod_sub_SER <- merge(mod_sub, SER_sub[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])
  mod_sub_SER$fit_LOG <- log(as.numeric(mod_sub_SER$fit))
  
  (SER_frame <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = mod_sub_SER$WATERBODY_CODE_SAMPLE_YEAR, 
                          OBS = mod_sub_SER$VALUE_LOG,
                          PRED = mod_sub_SER$fit_LOG))
  
  ## generate ML predictions meeting criteria
  mod_sub_ML <- merge(mod_sub, Hg_ML[[xx]]$pred_frame[,c("WATERBODY_CODE_SAMPLE_YEAR", "LMER_fit")])
  
  ML_frame <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = mod_sub_ML$WATERBODY_CODE_SAMPLE_YEAR, 
                         OBS = mod_sub_ML$VALUE_LOG,
                         PRED = mod_sub_ML$LMER_fit)
  
  ## generate AB predictions meeting criteria
  mod_sub_AB <- merge(mod_sub, Hg_AB[[xx]]$pred_frame[,c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")])
  
  AB_frame <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = mod_sub_AB$WATERBODY_CODE_SAMPLE_YEAR, 
                         OBS = mod_sub_AB$VALUE_LOG,
                         PRED = mod_sub_AB$`0.5quant`)
  
  ## generate B predictions meeting criteria
  mod_sub_B <- merge(mod_sub, Hg_B[[xx]]$pred_frame[,c("WATERBODY_CODE_SAMPLE_YEAR", "STAN_fit")])
  
  B_frame <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = mod_sub_B$WATERBODY_CODE_SAMPLE_YEAR, 
                        OBS = mod_sub_B$VALUE_LOG,
                        PRED = mod_sub_B$STAN_fit)
  
  
  ## Bring together for whole species
  spc_obspred_list <- list(SER_frame, ML_frame, AB_frame, B_frame)
  names(spc_obspred_list) <- c("SER", "ML", "AB", "B")
  
  ## Because of low total numbers of individuals within size class 1000 g +/-50 g 
  ## in each lake, decided to just use Global RMSE/R2 for prediction test
  table(table(spc_obspred_list[[4]]$WATERBODY_CODE_SAMPLE_YEAR))
  
  ## Global RMSE and R2 
  spc_obspred_list_global_stats <- lapply(spc_obspred_list, function(yy){
    
    (RMSE_global <- sqrt( sum( (yy$OBS - yy$PRED)^2 ) / nrow(yy) ))
    (R2_global <- summary(lm(yy$OBS ~ yy$PRED, data = mod_frame))$r.squared)
    
    data.frame(R2_global, RMSE_global)
    
  })
  
  return(spc_obspred_list_global_stats)
  
})
names(obspred_list) <- spc

## Build summary table

Hg_SER_summary

Hg_ML_mod_summ$LT[[1]] %>%
  mutate(SPECIES_CODE)



Hg_SER_summary <- Hg_SER %>% 
  group_by(SPECIES_CODE) %>% 
  mutate(RMSE = as.numeric(RMSE), 
         r2 = as.numeric(r2), 
         slp = as.numeric(slp), 
         int = as.numeric(int)) %>%
  summarise(RMSE_mean = mean(RMSE),
            RMSE_2p5 = quantile(RMSE, probs = 0.025),
            RMSE_97p5 = quantile(RMSE, probs = 0.975),
            R2_mean = mean(r2), 
            R2_2p5 = quantile(r2, probs = 0.025),
            R2_97p5 = quantile(r2, probs = 0.975), 
            sd_slope = sd(slp, na.rm = T), 
            sd_int = sd(int, na.rm = T))



Model | SER




par(mfrow = c(2,2))
lapply(obspred_list, function(yy){
  
  plot(yy$OBS ~ yy$PRED)
  
})





## Prediction test ---- 

## Hg ----

hg_sub_preds_mass$WATERBODY_CODE_SAMPLE_YEAR <- paste(hg_sub_preds_mass$WATERBODY_CODE, 
                                                      hg_sub_preds_mass$SAMPLE_YEAR, 
                                                      sep = "_")
Hg_LT_pred_list <- list()
Hg_NP_pred_list <- list()
Hg_WE_pred_list <- list()

## LT ----

Hg_LT_ind <- which(Hg_LT$WEIGHT_GRAM > 950 & Hg_LT$WEIGHT_GRAM < 1050)
Hg_LT_ind <- Hg_LT[Hg_LT_ind,]
nrow(Hg_LT_ind) # 773
nrow(Hg_LT) # 13870



## Compare with LMER

Hg_LT_LMER_test <- merge(Hg_LT_ind, cbind(Hg_LMER_LT_MASS_predict[, "WATERBODY_CODE_SAMPLE_YEAR"], Hg_LMER_LT_MASS_boot_pred_res_sum))
nrow(Hg_LT_LMER_test)

plot(Hg_LT_LMER_test$fit, Hg_LT_LMER_test$VALUE_LOG)
summary(lm(Hg_LT_LMER_test$VALUE_LOG ~ Hg_LT_LMER_test$fit)) # 0.75

Hg_LT_LMER_test_lm <- lm(Hg_LT_LMER_test$VALUE_LOG ~ Hg_LT_LMER_test$Hg_LMER_pred_LOG)
summary(Hg_LT_LMER_test_lm)
sqrt(mean(Hg_LT_LMER_test_lm$residuals^2))

Hg_LT_pred_list[["ML"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_LT_LMER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = Hg_LT_LMER_test$VALUE_LOG,
                                      PRED = Hg_LT_LMER_test$fit)

## Compare with INLA 
Hg_INLA_LT_MASS_predict

Hg_LT_INLA_test <- merge(Hg_LT_ind, Hg_INLA_LT_MASS_predict[, c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")] )
Hg_LT_INLA_test$Hg_INLA_pred_LOG <- log(Hg_LT_INLA_test$`0.5quant`)

Hg_LT_INLA_test_lm <- lm(Hg_LT_INLA_test$VALUE_LOG ~ Hg_LT_INLA_test$Hg_INLA_pred_LOG)
summary(Hg_LT_INLA_test_lm)
sqrt(mean(Hg_LT_INLA_test_lm$residuals^2))
nrow(Hg_LT_INLA_test)

Hg_LT_pred_list[["AB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_LT_INLA_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = Hg_LT_INLA_test$VALUE_LOG,
                                      PRED = Hg_LT_INLA_test$Hg_INLA_pred_LOG)

## Compare with STAN 

Hg_STAN_LT_MASS_predictions

Hg_LT_STAN_test <- merge(Hg_LT_ind, Hg_STAN_LT_MASS_predictions[, c("WATERBODY_CODE_SAMPLE_YEAR", "Hg_STAN_0.5quant")] )
Hg_LT_STAN_test$Hg_STAN_pred_LOG <- log(Hg_LT_STAN_test$Hg_STAN_0.5quant)

Hg_LT_STAN_test_lm <- lm(Hg_LT_STAN_test$VALUE_LOG ~ Hg_LT_STAN_test$Hg_STAN_pred_LOG)
summary(Hg_LT_STAN_test_lm)
sqrt(mean(Hg_LT_STAN_test_lm$residuals^2))
nrow(Hg_LT_STAN_test)

Hg_LT_pred_list[["MB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_LT_STAN_test$WATERBODY_CODE_SAMPLE_YEAR,
                                      OBS = Hg_LT_STAN_test$VALUE_LOG,
                                      PRED = Hg_LT_STAN_test$Hg_STAN_pred_LOG)

par(mfrow = c(2,4))
with(Hg_LT_pred_list[["SER"]], plot(OBS ~ PRED))
with(Hg_LT_pred_list[["ML"]], plot(OBS ~ PRED))
with(Hg_LT_pred_list[["AB"]], plot(OBS ~ PRED))
with(Hg_LT_pred_list[["AB"]], plot(OBS ~ PRED))

Hg_LT_pred_list

## Global
lapply(Hg_LT_pred_list, function(xx){
  
  sqrt( sum( (xx$OBS - xx$PRED)^2 ) / nrow(xx) )
  
})

## Local 
xx <- Hg_LT_pred_list[[1]]
yy <- un_event[1]
lapply(Hg_LT_pred_list, function(xx){
  
  un_event <- unique(xx$WATERBODY_CODE_SAMPLE_YEAR)
  
  res <- lapply(un_event, function(yy){
    
    print(yy)
    sub <- subset(xx, WATERBODY_CODE_SAMPLE_YEAR %in% yy)
    sub
    
    rmse <- sqrt( sum( (sub$OBS - sub$PRED)^2 ) / nrow(sub) )
    r2 <- summary(lm(OBS ~ PRED, data = sub))$r.squared
    
    data.frame(WATERBODY_CODE_SAMPLE_YEAR = yy, rmse, r2)
    
  })
  res
  rm(rmse, r2)
  
})

## NP ----

Hg_NP_ind <- which(Hg_NP$WEIGHT_GRAM > 950 & Hg_NP$WEIGHT_GRAM < 1050)
Hg_NP_ind <- Hg_NP[Hg_NP_ind,]
nrow(Hg_NP_ind) # 773
nrow(Hg_NP) # 13870

Hg_NP_pt_SER <- subset(hg_sub_preds_mass, (WATERBODY_CODE_SAMPLE_YEAR %in% Hg_NP_ind$WATERBODY_CODE_SAMPLE_YEAR) &
                         SPECIES_CODE == "131")

Hg_NP_SER_test <- merge(Hg_NP_ind, Hg_NP_pt_SER[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])
nrow(Hg_NP_SER_test)

Hg_NP_SER_test$fit_LOG <- log(as.numeric(Hg_NP_SER_test$fit))
plot(Hg_NP_SER_test$fit_LOG, Hg_NP_SER_test$VALUE_LOG)

Hg_NP_SER_test_lm <- lm(Hg_NP_SER_test$VALUE_LOG ~ Hg_NP_SER_test$fit_LOG)
summary(Hg_NP_SER_test_lm)
sqrt(mean(Hg_NP_SER_test_lm$residuals^2))

Hg_NP_pred_list[["SER"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_NP_SER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                       OBS = Hg_NP_SER_test$VALUE_LOG,
                                       PRED = Hg_NP_SER_test$fit_LOG)

## Compare with LMER

Hg_NP_LMER_test <- merge(Hg_NP_ind, cbind(Hg_LMER_NP_MASS_predict[, "WATERBODY_CODE_SAMPLE_YEAR"], Hg_LMER_NP_MASS_boot_pred_res_sum))
nrow(Hg_NP_LMER_test)

plot(Hg_NP_LMER_test$fit, Hg_NP_LMER_test$VALUE_LOG)
summary(lm(Hg_NP_LMER_test$VALUE_LOG ~ Hg_NP_LMER_test$fit)) # 0.75

Hg_NP_LMER_test_lm <- lm(Hg_NP_LMER_test$VALUE_LOG ~ Hg_NP_LMER_test$Hg_LMER_pred_LOG)
summary(Hg_NP_LMER_test_lm)
sqrt(mean(Hg_NP_LMER_test_lm$residuals^2))

Hg_NP_pred_list[["ML"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_NP_LMER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = Hg_NP_LMER_test$VALUE_LOG,
                                      PRED = Hg_NP_LMER_test$fit)

## Compare with INLA 
Hg_INLA_NP_MASS_predict

Hg_NP_INLA_test <- merge(Hg_NP_ind, Hg_INLA_NP_MASS_predict[, c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")] )
Hg_NP_INLA_test$Hg_INLA_pred_LOG <- log(Hg_NP_INLA_test$`0.5quant`)

Hg_NP_INLA_test_lm <- lm(Hg_NP_INLA_test$VALUE_LOG ~ Hg_NP_INLA_test$Hg_INLA_pred_LOG)
summary(Hg_NP_INLA_test_lm)
sqrt(mean(Hg_NP_INLA_test_lm$residuals^2))
nrow(Hg_NP_INLA_test)

Hg_NP_pred_list[["AB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_NP_INLA_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = Hg_NP_INLA_test$VALUE_LOG,
                                      PRED = Hg_NP_INLA_test$Hg_INLA_pred_LOG)

## Compare with STAN 

Hg_STAN_NP_MASS_predictions

Hg_NP_STAN_test <- merge(Hg_NP_ind, Hg_STAN_NP_MASS_predictions[, c("WATERBODY_CODE_SAMPLE_YEAR", "Hg_STAN_0.5quant")] )
Hg_NP_STAN_test$Hg_STAN_pred_LOG <- log(Hg_NP_STAN_test$Hg_STAN_0.5quant)

Hg_NP_STAN_test_lm <- lm(Hg_NP_STAN_test$VALUE_LOG ~ Hg_NP_STAN_test$Hg_STAN_pred_LOG)
summary(Hg_NP_STAN_test_lm)
sqrt(mean(Hg_NP_STAN_test_lm$residuals^2))
nrow(Hg_NP_STAN_test)

Hg_NP_pred_list[["MB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_NP_STAN_test$WATERBODY_CODE_SAMPLE_YEAR,
                                      OBS = Hg_NP_STAN_test$VALUE_LOG,
                                      PRED = Hg_NP_STAN_test$Hg_STAN_pred_LOG)

par(mfrow = c(2,4))
with(Hg_NP_pred_list[["SER"]], plot(OBS ~ PRED))
with(Hg_NP_pred_list[["ML"]], plot(OBS ~ PRED))
with(Hg_NP_pred_list[["AB"]], plot(OBS ~ PRED))
with(Hg_NP_pred_list[["AB"]], plot(OBS ~ PRED))

## WE ----

Hg_WE_ind <- which(Hg_WE$WEIGHT_GRAM > 950 & Hg_WE$WEIGHT_GRAM < 1050)
Hg_WE_ind <- Hg_WE[Hg_WE_ind,]
nrow(Hg_WE_ind) # 773
nrow(Hg_WE) # 13870

Hg_WE_pt_SER <- subset(hg_sub_preds_mass, (WATERBODY_CODE_SAMPLE_YEAR %in% Hg_WE_ind$WATERBODY_CODE_SAMPLE_YEAR) &
                         SPECIES_CODE == "334")

Hg_WE_SER_test <- merge(Hg_WE_ind, Hg_WE_pt_SER[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])
nrow(Hg_WE_SER_test)

Hg_WE_SER_test$fit_LOG <- log(as.numeric(Hg_WE_SER_test$fit))
plot(Hg_WE_SER_test$fit_LOG, Hg_WE_SER_test$VALUE_LOG)

Hg_WE_SER_test_lm <- lm(Hg_WE_SER_test$VALUE_LOG ~ Hg_WE_SER_test$fit_LOG)
summary(Hg_WE_SER_test_lm)
sqrt(mean(Hg_WE_SER_test_lm$residuals^2))

Hg_WE_pred_list[["SER"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_WE_SER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                       OBS = Hg_WE_SER_test$VALUE_LOG,
                                       PRED = Hg_WE_SER_test$fit_LOG)

## Compare with LMER

Hg_WE_LMER_test <- merge(Hg_WE_ind, cbind(Hg_LMER_WE_MASS_predict[, "WATERBODY_CODE_SAMPLE_YEAR"], Hg_LMER_WE_MASS_boot_pred_res_sum))
nrow(Hg_WE_LMER_test)

plot(Hg_WE_LMER_test$fit, Hg_WE_LMER_test$VALUE_LOG)
summary(lm(Hg_WE_LMER_test$VALUE_LOG ~ Hg_WE_LMER_test$fit)) # 0.75

Hg_WE_LMER_test_lm <- lm(Hg_WE_LMER_test$VALUE_LOG ~ Hg_WE_LMER_test$Hg_LMER_pred_LOG)
summary(Hg_WE_LMER_test_lm)
sqrt(mean(Hg_WE_LMER_test_lm$residuals^2))

Hg_WE_pred_list[["ML"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_WE_LMER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = Hg_WE_LMER_test$VALUE_LOG,
                                      PRED = Hg_WE_LMER_test$fit)

## Compare with INLA 
Hg_INLA_WE_MASS_predict

Hg_WE_INLA_test <- merge(Hg_WE_ind, Hg_INLA_WE_MASS_predict[, c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")] )
Hg_WE_INLA_test$Hg_INLA_pred_LOG <- log(Hg_WE_INLA_test$`0.5quant`)

Hg_WE_INLA_test_lm <- lm(Hg_WE_INLA_test$VALUE_LOG ~ Hg_WE_INLA_test$Hg_INLA_pred_LOG)
summary(Hg_WE_INLA_test_lm)
sqrt(mean(Hg_WE_INLA_test_lm$residuals^2))
nrow(Hg_WE_INLA_test)

Hg_WE_pred_list[["AB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_WE_INLA_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = Hg_WE_INLA_test$VALUE_LOG,
                                      PRED = Hg_WE_INLA_test$Hg_INLA_pred_LOG)

## Compare with STAN 

Hg_STAN_WE_MASS_predictions

Hg_WE_STAN_test <- merge(Hg_WE_ind, Hg_STAN_WE_MASS_predictions[, c("WATERBODY_CODE_SAMPLE_YEAR", "Hg_STAN_0.5quant")] )
Hg_WE_STAN_test$Hg_STAN_pred_LOG <- log(Hg_WE_STAN_test$Hg_STAN_0.5quant)

Hg_WE_STAN_test_lm <- lm(Hg_WE_STAN_test$VALUE_LOG ~ Hg_WE_STAN_test$Hg_STAN_pred_LOG)
summary(Hg_WE_STAN_test_lm)
sqrt(mean(Hg_WE_STAN_test_lm$residuals^2))
nrow(Hg_WE_STAN_test)

Hg_WE_pred_list[["MB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_WE_STAN_test$WATERBODY_CODE_SAMPLE_YEAR,
                                      OBS = Hg_WE_STAN_test$VALUE_LOG,
                                      PRED = Hg_WE_STAN_test$Hg_STAN_pred_LOG)

par(mfrow = c(2,4))
with(Hg_WE_pred_list[["SER"]], plot(OBS ~ PRED))
with(Hg_WE_pred_list[["ML"]], plot(OBS ~ PRED))
with(Hg_WE_pred_list[["AB"]], plot(OBS ~ PRED))
with(Hg_WE_pred_list[["AB"]], plot(OBS ~ PRED))


## As ----

As_sub_preds_mass$WATERBODY_CODE_SAMPLE_YEAR <- paste(As_sub_preds_mass$WATERBODY_CODE, 
                                                      As_sub_preds_mass$SAMPLE_YEAR, sep = "_")

As_LT_pred_list <- list()
As_NP_pred_list <- list()
As_WE_pred_list <- list()

## LT ----

As_LT_ind <- which(As_LT$WEIGHT_GRAM > 950 & As_LT$WEIGHT_GRAM < 1050)
As_LT_ind <- As_LT[As_LT_ind,]
nrow(As_LT_ind) # 773
nrow(As_LT) # 13870

(As_LT_pt_SER <- subset(As_sub_preds_mass, (WATERBODY_CODE_SAMPLE_YEAR %in% As_LT_ind$WATERBODY_CODE_SAMPLE_YEAR) &
                          SPECIES_CODE == "LT"))

As_LT_SER_test <- merge(As_LT_ind, As_LT_pt_SER[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])
nrow(As_LT_SER_test)

As_LT_SER_test$fit_LOG <- log(as.numeric(As_LT_SER_test$fit))
plot(As_LT_SER_test$fit_LOG, As_LT_SER_test$VALUE_LOG)

As_LT_SER_test_lm <- lm(As_LT_SER_test$VALUE_LOG ~ As_LT_SER_test$fit_LOG)
summary(As_LT_SER_test_lm)
sqrt(mean(As_LT_SER_test_lm$residuals^2))

As_LT_pred_list[["SER"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_LT_SER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                       OBS = As_LT_SER_test$VALUE_LOG,
                                       PRED = As_LT_SER_test$fit_LOG)

## Compare with LMER

As_LT_LMER_test <- merge(As_LT_ind, 
                         data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_LMER_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR,
                                    As_LMER_LT_MASS_boot_pred_res_sum))

plot(As_LT_LMER_test$fit, As_LT_LMER_test$VALUE_LOG)
summary(lm(As_LT_LMER_test$VALUE_LOG ~ As_LT_LMER_test$fit)) # 0.75

As_LT_LMER_test_lm <- lm(As_LT_LMER_test$VALUE_LOG ~ As_LT_LMER_test$fit)
summary(As_LT_LMER_test_lm)
sqrt(mean(As_LT_LMER_test_lm$residuals^2))

As_LT_pred_list[["ML"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_LT_LMER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = As_LT_LMER_test$VALUE_LOG,
                                      PRED = As_LT_LMER_test$fit)

## Compare with INLA 
As_INLA_LT_MASS_predict

As_LT_INLA_test <- merge(As_LT_ind, As_INLA_LT_MASS_predict[, c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")] )
As_LT_INLA_test$As_INLA_pred_LOG <- log(As_LT_INLA_test$`0.5quant`)

As_LT_INLA_test_lm <- lm(As_LT_INLA_test$VALUE_LOG ~ As_LT_INLA_test$As_INLA_pred_LOG)
summary(As_LT_INLA_test_lm)
sqrt(mean(As_LT_INLA_test_lm$residuals^2))
nrow(As_LT_INLA_test)

As_LT_pred_list[["AB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_LT_INLA_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = As_LT_INLA_test$VALUE_LOG,
                                      PRED = As_LT_INLA_test$As_INLA_pred_LOG)

## Compare with STAN 

As_LT_STAN_test <- merge(As_LT_ind, As_STAN_LT_MASS_predictions[, c("WATERBODY_CODE_SAMPLE_YEAR", "As_STAN_0.5quant")] )
As_LT_STAN_test$As_STAN_pred_LOG <- log(As_LT_STAN_test$As_STAN_0.5quant)

As_LT_STAN_test_lm <- lm(As_LT_STAN_test$VALUE_LOG ~ As_LT_STAN_test$As_STAN_pred_LOG)
summary(As_LT_STAN_test_lm)
sqrt(mean(As_LT_STAN_test_lm$residuals^2))
nrow(As_LT_STAN_test)

As_LT_pred_list[["MB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_LT_STAN_test$WATERBODY_CODE_SAMPLE_YEAR,
                                      OBS = As_LT_STAN_test$VALUE_LOG,
                                      PRED = As_LT_STAN_test$As_STAN_pred_LOG)

par(mfrow = c(2,4))
with(As_LT_pred_list[["SER"]], plot(OBS ~ PRED))
with(As_LT_pred_list[["ML"]], plot(OBS ~ PRED))
with(As_LT_pred_list[["AB"]], plot(OBS ~ PRED))
with(As_LT_pred_list[["AB"]], plot(OBS ~ PRED))

## NP ----

As_NP_ind <- which(As_NP$WEIGHT_GRAM > 950 & As_NP$WEIGHT_GRAM < 1050)
As_NP_ind <- As_NP[As_NP_ind,]
nrow(As_NP_ind) # 773
nrow(As_NP) # 13870

unique(As_NP_ind)

subset(As_sub_preds_mass, SPECIES_CODE == "NP")

As_sub_preds_mass

As_NP_pt_SER <- subset(As_sub_preds_mass, (WATERBODY_CODE_SAMPLE_YEAR %in% As_NP_ind$WATERBODY_CODE_SAMPLE_YEAR) &
                         SPECIES_CODE == "NP")
As_NP_SER_test <- merge(As_NP_ind, As_NP_pt_SER[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])

As_NP_SER_test$fit_LOG <- log(as.numeric(As_NP_SER_test$fit))
plot(As_NP_SER_test$fit_LOG, As_NP_SER_test$VALUE_LOG)

As_NP_SER_test_lm <- lm(As_NP_SER_test$VALUE_LOG ~ As_NP_SER_test$fit_LOG)
summary(As_NP_SER_test_lm)
sqrt(mean(As_NP_SER_test_lm$residuals^2))

As_NP_pred_list[["SER"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_NP_SER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                       OBS = As_NP_SER_test$VALUE_LOG,
                                       PRED = As_NP_SER_test$fit_LOG)

## Compare with LMER

As_NP_LMER_test <- merge(As_NP_ind, 
                         data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_LMER_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR,
                                    As_LMER_NP_MASS_boot_pred_res_sum))
nrow(As_NP_LMER_test)

plot(As_NP_LMER_test$fit, As_NP_LMER_test$VALUE_LOG)
summary(lm(As_NP_LMER_test$VALUE_LOG ~ As_NP_LMER_test$fit)) # 0.75

As_NP_LMER_test_lm <- lm(As_NP_LMER_test$VALUE_LOG ~ As_NP_LMER_test$As_LMER_pred_LOG)
summary(As_NP_LMER_test_lm)
sqrt(mean(As_NP_LMER_test_lm$residuals^2))

As_NP_pred_list[["ML"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_NP_LMER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = As_NP_LMER_test$VALUE_LOG,
                                      PRED = As_NP_LMER_test$fit)

## Compare with INLA 
As_INLA_NP_MASS_predict

As_NP_INLA_test <- merge(As_NP_ind, As_INLA_NP_MASS_predict[, c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")] )
As_NP_INLA_test$As_INLA_pred_LOG <- log(As_NP_INLA_test$`0.5quant`)

As_NP_INLA_test_lm <- lm(As_NP_INLA_test$VALUE_LOG ~ As_NP_INLA_test$As_INLA_pred_LOG)
summary(As_NP_INLA_test_lm)
sqrt(mean(As_NP_INLA_test_lm$residuals^2))
nrow(As_NP_INLA_test)

As_NP_pred_list[["AB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_NP_INLA_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = As_NP_INLA_test$VALUE_LOG,
                                      PRED = As_NP_INLA_test$As_INLA_pred_LOG)

## Compare with STAN 

As_STAN_NP_MASS_predictions

As_NP_STAN_test <- merge(As_NP_ind, As_STAN_NP_MASS_predictions[, c("WATERBODY_CODE_SAMPLE_YEAR", "As_STAN_0.5quant")] )
As_NP_STAN_test$As_STAN_pred_LOG <- log(As_NP_STAN_test$As_STAN_0.5quant)

As_NP_STAN_test_lm <- lm(As_NP_STAN_test$VALUE_LOG ~ As_NP_STAN_test$As_STAN_pred_LOG)
summary(As_NP_STAN_test_lm)
sqrt(mean(As_NP_STAN_test_lm$residuals^2))
nrow(As_NP_STAN_test)

As_NP_pred_list[["MB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_NP_STAN_test$WATERBODY_CODE_SAMPLE_YEAR,
                                      OBS = As_NP_STAN_test$VALUE_LOG,
                                      PRED = As_NP_STAN_test$As_STAN_pred_LOG)

par(mfrow = c(2,4))
with(As_NP_pred_list[["SER"]], plot(OBS ~ PRED))
with(As_NP_pred_list[["ML"]], plot(OBS ~ PRED))
with(As_NP_pred_list[["AB"]], plot(OBS ~ PRED))
with(As_NP_pred_list[["MB"]], plot(OBS ~ PRED))

head(As_NP_pred_list[["SER"]])
head(As_NP_pred_list[["ML"]])
head(As_NP_pred_list[["AB"]])


## WE ----

As_WE_ind <- which(As_WE$WEIGHT_GRAM > 950 & As_WE$WEIGHT_GRAM < 1050)
As_WE_ind <- As_WE[As_WE_ind,]
nrow(As_WE_ind) # 773
nrow(As_WE) # 13870

As_WE_pt_SER <- subset(As_sub_preds_mass, (WATERBODY_CODE_SAMPLE_YEAR %in% As_WE_ind$WATERBODY_CODE_SAMPLE_YEAR) &
                         SPECIES_CODE == "334")

As_WE_SER_test <- merge(As_WE_ind, As_WE_pt_SER[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])
nrow(As_WE_SER_test)

As_WE_SER_test$fit_LOG <- log(as.numeric(As_WE_SER_test$fit))
plot(As_WE_SER_test$fit_LOG, As_WE_SER_test$VALUE_LOG)

As_WE_SER_test_lm <- lm(As_WE_SER_test$VALUE_LOG ~ As_WE_SER_test$fit_LOG)
summary(As_WE_SER_test_lm)
sqrt(mean(As_WE_SER_test_lm$residuals^2))

As_WE_pred_list[["SER"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_WE_SER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                       OBS = As_WE_SER_test$VALUE_LOG,
                                       PRED = As_WE_SER_test$fit_LOG)

## Compare with LMER

As_WE_LMER_test <- merge(As_WE_ind, cbind(As_LMER_WE_MASS_predict[, "WATERBODY_CODE_SAMPLE_YEAR"], As_LMER_WE_MASS_boot_pred_res_sum))
nrow(As_WE_LMER_test)

plot(As_WE_LMER_test$fit, As_WE_LMER_test$VALUE_LOG)
summary(lm(As_WE_LMER_test$VALUE_LOG ~ As_WE_LMER_test$fit)) # 0.75

As_WE_LMER_test_lm <- lm(As_WE_LMER_test$VALUE_LOG ~ As_WE_LMER_test$As_LMER_pred_LOG)
summary(As_WE_LMER_test_lm)
sqrt(mean(As_WE_LMER_test_lm$residuals^2))

As_WE_pred_list[["ML"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_WE_LMER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = As_WE_LMER_test$VALUE_LOG,
                                      PRED = As_WE_LMER_test$fit)

## Compare with INLA 
As_INLA_WE_MASS_predict

As_WE_INLA_test <- merge(As_WE_ind, As_INLA_WE_MASS_predict[, c("WATERBODY_CODE_SAMPLE_YEAR", "0.5quant")] )
As_WE_INLA_test$As_INLA_pred_LOG <- log(As_WE_INLA_test$`0.5quant`)

As_WE_INLA_test_lm <- lm(As_WE_INLA_test$VALUE_LOG ~ As_WE_INLA_test$As_INLA_pred_LOG)
summary(As_WE_INLA_test_lm)
sqrt(mean(As_WE_INLA_test_lm$residuals^2))
nrow(As_WE_INLA_test)

As_WE_pred_list[["AB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_WE_INLA_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                      OBS = As_WE_INLA_test$VALUE_LOG,
                                      PRED = As_WE_INLA_test$As_INLA_pred_LOG)

## Compare with STAN 

As_STAN_WE_MASS_predictions

As_WE_STAN_test <- merge(As_WE_ind, As_STAN_WE_MASS_predictions[, c("WATERBODY_CODE_SAMPLE_YEAR", "As_STAN_0.5quant")] )
As_WE_STAN_test$As_STAN_pred_LOG <- log(As_WE_STAN_test$As_STAN_0.5quant)

As_WE_STAN_test_lm <- lm(As_WE_STAN_test$VALUE_LOG ~ As_WE_STAN_test$As_STAN_pred_LOG)
summary(As_WE_STAN_test_lm)
sqrt(mean(As_WE_STAN_test_lm$residuals^2))
nrow(As_WE_STAN_test)

As_WE_pred_list[["MB"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = As_WE_STAN_test$WATERBODY_CODE_SAMPLE_YEAR,
                                      OBS = As_WE_STAN_test$VALUE_LOG,
                                      PRED = As_WE_STAN_test$As_STAN_pred_LOG)

par(mfrow = c(2,4))
with(As_WE_pred_list[["SER"]], plot(OBS ~ PRED))
with(As_WE_pred_list[["ML"]], plot(OBS ~ PRED))
with(As_WE_pred_list[["AB"]], plot(OBS ~ PRED))
with(As_WE_pred_list[["AB"]], plot(OBS ~ PRED))












































## Some investigative plots
lapply(Hg_ML, function(xx){
  
  mod_resid <- (resid(xx$model, type = "pearson"))
  mod_fitted <- fitted(xx$model)
  ranef_wb <- (ranef(xx$model)$WATERBODY_CODE$`(Intercept)`)
  ranef_wbsy <- (ranef(xx$model)$WATERBODY_CODE_SAMPLE_YEAR$`(Intercept)`)
  ranef_slope <- (ranef(xx$model)$WATERBODY_CODE$WEIGHT_GRAM_LOG) 
  
  p1 <- ggplot() + geom_histogram(aes(mod_resid), bins = 100) + 
    theme_classic() + xlab("Residuals")
  p2 <- ggplot() + geom_point(aes(x = mod_fitted, y = mod_resid)) + 
    theme_classic() + xlab("Fitted") + ylab("Residuals")
  p3 <- ggplot() + geom_histogram(aes(ranef_wb), bins = 100) + 
    theme_classic() + xlab("WATERBODY")
  p4 <- ggplot() + geom_histogram(aes(ranef_wbsy), bins = 100) + 
    theme_classic() + xlab("WATERBODY:SAMPLE YEAR")
  p5 <- ggplot() + geom_histogram(aes(ranef_slope), bins = 100) + 
    xlab("SLOPE") + theme_classic()
  
  
  (p1 + p2) / (p3 + p4 + p5)
})




## ***********************
## Global RMSE and R2 ----
## ***********************  









## Global RMSE and R2 ----