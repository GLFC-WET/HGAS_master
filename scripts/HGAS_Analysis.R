## HGAS_Analysis
## Author(s): Brian Kielstra
## Originated: 2021-10-14
##
##
## Investigates best methods for estimating lake-level Hg and As concentrations
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

## **************************************************
## STANDARDIZATION APPROACHES & DATA EXPLORATION ----
## **************************************************

## If standardizing Hg concentration by fish length/mass, what is the best way?
## There are a few options tried here:
## 1) Sampling event-based regressions (lake-year-species) and
##    predict to common length or mass by species
## 2) Global mixed model

## An advange to the sampling event regressions are that you're honing in
## on an individual lake and might get more accurate predictions for a specific
## lake. A disadvantage is that the fit might not be very good on a sampling
## event basis or you need to throw out an entire regression because there
## were not enough data to fit the relationships (here, at least 5 fish were needed)

## An advantage to the global mixed model is that one can borrow from
## other observations to produce a prediction even if only one fish were in the
## lake.

## 1) HG ----

## Read the saved file
hg <- readxl::read_excel("./data/Ontario Inland 3 Species Hg 2008-2020-12-16.xlsx")
nrow(hg) #39319

sum( (hg$WEIGHT_GRAM > 900 & hg$WEIGHT_GRAM < 1100) , na.rm = TRUE) # 3843 fish

## Log transform
hg$LENGTH_CM_LOG <- log(hg$LENGTH_CM) 
hg$WEIGHT_GRAM_LOG <- log(hg$WEIGHT_GRAM)
hg$VALUE_LOG <- log(hg$VALUE)

hg$WATERBODY_CODE_SAMPLE_YEAR <- paste(hg$WATERBODY_CODE, hg$SAMPLE_YEAR, sep = "_")

## Subset by "SKINLESS, BONELESS FILLET" OR "SKINLESS, BONELESS FILLET PLUG"
hg_sub <- subset(hg, PORTION_TYPE_DESC %in% 
                         c("SKINLESS, BONELESS FILLET (STANDARD MOE DORSAL FILLET)", 
                           "SKINLESS, BONELESS FILLET PLUG (SKIN-OFF)"))

## Subset and remove those lengths and weights that are NA
hg_sub <- subset(hg_sub, !is.na(LENGTH_CM_LOG))
hg_sub <- subset(hg_sub, !is.na(WEIGHT_GRAM_LOG))

length(unique(hg_sub$LOCATION_NAME)) # 1197 lakes
length(unique(hg_sub$WATERBODY_CODE_SAMPLE_YEAR)) # 1915 sampling events

unique(hg_sub[, c("SPECIES_CODE", "SPECIES_NAME")])

nrow(hg_sub)

## 1.x) Data exploration for VALUE (Hg), WEIGHT_GRAM, and LENGTH_CM ----
hist(hg_sub$VALUE); summary(hg_sub$VALUE)
dotchart(hg_sub$VALUE) # majority of data < 1 but some outliers 
dotchart(log(hg_sub$VALUE)) # few outliers and majority of data b/w -1 and 0 as expected 
hist(log(hg_sub$VALUE))

hist(hg_sub$WEIGHT_GRAM); summary(hg_sub$WEIGHT_GRAM)
dotchart(hg_sub$WEIGHT_GRAM) # majority of data < 2000 but some outliers 
dotchart(log(hg_sub$WEIGHT_GRAM)) # few outliers and majority of data b/w -1 and 0 as expected 
hist(log(hg_sub$WEIGHT_GRAM))

hist(hg_sub$LENGTH_CM); summary(hg_sub$LENGTH_CM)
dotchart(hg_sub$LENGTH_CM) # majority of data < 2000 but some outliers 
dotchart(log(hg_sub$LENGTH_CM)) # few outliers and majority of data b/w -1 and 0 as expected 
hist(log(hg_sub$LENGTH_CM))

## 1.x) All species w/ dataset-based median length, mass ----

median_LENGTH_CM <- sapply(unique(hg_sub$SPECIES_CODE), function(x) {
  hg_sub_sub <- subset(hg_sub, SPECIES_CODE == x)
  res <- median(hg_sub_sub$LENGTH_CM, na.rm = T)
  return(res)
}, USE.NAMES = TRUE)

hist(hg_sub$LENGTH_CM, breaks = 100); abline(v = median_LENGTH_CM, lty = 2, col = "red")

median_WEIGHT_GRAM <- sapply(unique(hg_sub$SPECIES_CODE), function(x) {
  hg_sub_sub <- subset(hg_sub, SPECIES_CODE == x)
  res <- median(hg_sub_sub$WEIGHT_GRAM, na.rm = TRUE)
  return(res)
}, USE.NAMES = TRUE)

## Get median size for all species and construct prediction dataframe
median_frame <- data.frame(
  SPECIES_CODE = unique(hg_sub$SPECIES_CODE),
  LENGTH_CM = median_LENGTH_CM,
  WEIGHT_GRAM = median_WEIGHT_GRAM,
  LENGTH_CM_LOG = log(median_LENGTH_CM), 
  WEIGHT_GRAM_LOG = log(median_WEIGHT_GRAM), 
  stringsAsFactors = FALSE
)

Hg_LT <- subset(hg_sub, SPECIES_NAME == "Lake Trout")
Hg_NP <- subset(hg_sub, SPECIES_NAME == "Northern Pike")
Hg_WE <- subset(hg_sub, SPECIES_NAME == "Walleye")

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

## 1.x) SER - sampling event regressions ----

 rm(x, y, z)
 x = unique(hg$WATERBODY_CODE)[1]
 y = unique(sub_waterbody$SAMPLE_YEAR)[1]
 z = unique(sub_year$SPECIES_CODE)[1]

## Standardize by mass, >4 individuals needed and using Rob Mackereth's standardized values
hg_sub_preds_mass <- lapply(unique(hg_sub$WATERBODY_CODE), function(x) {
  
  ## Subsets the waterbody of interest
  sub_waterbody <- subset(hg_sub, WATERBODY_CODE == x)
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$SAMPLE_YEAR), function(y) {
    
    ## Subsets a year within a given waterbody
    sub_year <- subset(sub_waterbody, SAMPLE_YEAR == y)
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$SPECIES_CODE), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, SPECIES_CODE == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      ## the sum evaluation statement is meant to make sure that all
      ## masses are not NAs which was a case in these data
      if (nrow(sub_species) > 4 & !sum(sub_species$WEIGHT_GRAM, na.rm = T) == 0) {
        
        ## log mercury concentrations and masss (converted to mm)
        (logm <- log(sub_species$VALUE))
        (logl <- log(sub_species$WEIGHT_GRAM))
        
        percentile <- ecdf(sub_species$WEIGHT_GRAM)
        
        ## regression relationship
        rel <- lm(logm ~ logl)
        
        #plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        #abline(rel)
        
        ## Generating context and regression summary statistics
        (WATERBODY_CODE <- x) # waterbody
        (SAMPLE_YEAR <- y) # year
        (SPECIES_CODE <- z) # species
        (n <- length(logl)) # number of individuals
        (int <- formatC(rel$coefficients[1], digits = 4, format = "f")) # estimated intercept
        (slp <- formatC(rel$coefficients[2], digits = 4, format = "f")) # estimated slope
        (int_confint <- paste(formatC(confint(rel)[1, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated intercept confidence interval
        (slp_confint <- paste(formatC(confint(rel)[2, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated slope confidence interval
        (r2 <- formatC(summary(rel)$r.squared, digits = 4, format = "f")) # adjusted R2
        (RMSE <- sqrt(mean(rel$residuals^2)))
        
        ## Modify the predicted standardized length based on species
        ## (pred_modify <- median_frame[median_frame$SPECIES_CODE == z, "WEIGHT_GRAM"])
        pred_modify <- 1000
        
        (target_size_percentile <- percentile(pred_modify))
        
        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
                             newdata = data.frame(logl = log(pred_modify)),
                             interval = "confidence"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(WATERBODY_CODE, SAMPLE_YEAR, SPECIES_CODE, n,
                             int, int_confint, slp, slp_confint, r2, RMSE,
                             target_size_percentile, pred,
                             stringsAsFactors = FALSE
        ))
      } else {
        frame <- NULL
      }
      return(frame)
    })
    
    ## rbind species-specific regression relationships and predictions
    sub_species_frame <- do.call(rbind, sub_species_output)
    return(sub_species_frame)
  })
  ## rbind years for a waterbody
  sub_year_frame <- do.call(rbind, sub_year_output)
  sub_year_frame
  return(sub_year_frame)
})
hg_sub_preds_mass <- hg_sub_preds_mass[!sapply(hg_sub_preds_mass, is.null)]
hg_sub_preds_mass <- dplyr::bind_rows(hg_sub_preds_mass)

## 1.x) lmer - mixed model, maximum likelihood ----

## Some functions for easier processing
LMER_boot_est <- function(.) {
  c(beta=fixef(.), 
    as.data.frame(VarCorr(.))$sdcor[c(1:3,5)])
}

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

LMER_boot_pred_summary <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

LMER_boot_initiate <- function(varlist, nclust){
  closeAllConnections()
  clust <- parallel::makeCluster(nclust)
  parallel::clusterEvalQ(clust, library("lme4"))
  parallel::clusterExport(cl = clust, varlist = varlist)
  showConnections()
  clust
}

## 1.x.x) Lake Trout ----

Hg_LMER_LT_MASS <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                          (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                          (1|WATERBODY_CODE_SAMPLE_YEAR), data = Hg_LT, 
                        control = lmerControl(optimizer = "Nelder_Mead")) 
summary(Hg_LMER_LT_MASS)
summary(lm(Hg_LT$VALUE_LOG ~ fitted(Hg_LMER_LT_MASS)))

hist(resid(Hg_LMER_LT_MASS), breaks = 100) # seems good 
plot(resid(Hg_LMER_LT_MASS) ~ fitted(Hg_LMER_LT_MASS)) # seems reasonable 

## Generate predictions for each waterbody based on the median size of LT in whole dataset 
Hg_LMER_LT_MASS_predict <- unique(Hg_LT[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
Hg_LMER_LT_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

## Bootstrap confidence intervals from bootMer
clust <- LMER_boot_initiate(varlist = "Hg_LMER_LT_MASS_predict")
system.time(Hg_LMER_LT_MASS_boot_est <- lme4::bootMer(Hg_LMER_LT_MASS, LMER_boot_est, 
                                                      nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                      parallel = "snow", 
                                                      cl = clust, 
                                                      ncpus = 4)) #17s for 100

LMER_boot_est_summary(Hg_LMER_LT_MASS_boot_est) ## 889s

Hg_LMER_LT_MASS_boot_pred <- function(., newdata) {
  predict(., newdata=Hg_LMER_LT_MASS_predict)
}

clust <- LMER_boot_initiate(varlist = "Hg_LMER_LT_MASS_predict")
system.time(Hg_LMER_LT_MASS_boot_pred_res <- lme4::bootMer(Hg_LMER_LT_MASS, Hg_LMER_LT_MASS_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 

(Hg_LMER_LT_MASS_boot_pred_res_sum <- LMER_boot_pred_summary(Hg_LMER_LT_MASS_boot_pred_res))

## Results 
LMER_boot_est_summary(Hg_LMER_LT_MASS_boot_est_res)
head(Hg_LMER_LT_MASS_boot_pred_res_sum)

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/glmmHg_LMER/vignettes/model_evaluation.pdf

## Easy to generate effects plots with glmmHg_LMER 
ae_MASS <- allEffects(Hg_LMER_LT_MASS)

par(mfrow = c(1,2))
plot(ae_MASS$WEIGHT_GRAM_LOG$data$VALUE_LOG ~ ae_MASS$WEIGHT_GRAM_LOG$data$WEIGHT_GRAM_LOG, main = "MASS")  
lines(ae_MASS$WEIGHT_GRAM_LOG$fit ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$upper ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, lty = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$lower ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, lty = 2)

(pp <- sjPlot::plot_model(Hg_LMER_LT_MASS,type=c("pred"),
                          terms=c("WEIGHT_GRAM_LOG","WATERBODY_CODE"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))

## 1.x.x) Northern Pike ----

Hg_LMER_NP_MASS <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                          (1|WATERBODY_CODE_SAMPLE_YEAR), data = Hg_NP,
                        control = lmerControl(optimizer = "Nelder_Mead")) 
summary(Hg_LMER_NP_MASS)
summary(lm(Hg_NP$VALUE_LOG ~ fitted(Hg_LMER_NP_MASS)))

hist(resid(Hg_LMER_NP_MASS), breaks = 100) # seems good 
plot(resid(Hg_LMER_NP_MASS) ~ fitted(Hg_LMER_NP_MASS)) # seems reasonable 

## Generate predictions for each waterbody based on the median size of NP in whole dataset 
Hg_LMER_NP_MASS_predict <- unique(Hg_NP[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
Hg_LMER_NP_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

## Bootstrap confidence intervals from bootMer
clust <- LMER_boot_initiate(varlist = "Hg_LMER_NP_MASS_predict")
system.time(Hg_LMER_NP_MASS_boot_est <- lme4::bootMer(Hg_LMER_NP_MASS, LMER_boot_est, 
                                                      nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                      parallel = "snow", 
                                                      cl = clust, 
                                                      ncpus = 4)) #17s for 100
LMER_boot_est_summary(Hg_LMER_NP_MASS_boot_est) ## 889s

Hg_LMER_NP_MASS_boot_pred <- function(., newdata) {
  predict(., newdata=Hg_LMER_NP_MASS_predict)
}

clust <- LMER_boot_initiate(varlist = "Hg_LMER_NP_MASS_predict")
system.time(Hg_LMER_NP_MASS_boot_pred_res <- lme4::bootMer(Hg_LMER_NP_MASS, Hg_LMER_NP_MASS_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 

(Hg_LMER_NP_MASS_boot_pred_res_sum <- LMER_boot_pred_summary(Hg_LMER_NP_MASS_boot_pred_res))

## ResuNPs 
LMER_boot_est_summary(Hg_LMER_NP_MASS_boot_est_res)
head(Hg_LMER_NP_MASS_boot_pred_res_sum)

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/glmmHg_LMER/vignettes/model_evaluation.pdf

## Easy to generate effects plots with glmmHg_LMER 
ae_MASS <- allEffects(Hg_LMER_NP_MASS)

par(mfrow = c(1,2))
plot(ae_MASS$WEIGHT_GRAM_LOG$data$VALUE_LOG ~ ae_MASS$WEIGHT_GRAM_LOG$data$WEIGHT_GRAM_LOG, main = "MASS")  
lines(ae_MASS$WEIGHT_GRAM_LOG$fit ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$upper ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, NPy = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$lower ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, NPy = 2)

(pp <- sjPlot::plot_model(Hg_LMER_NP_MASS,type=c("pred"),
                          terms=c("WEIGHT_GRAM_LOG","WATERBODY_CODE"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))

## 1.x.x) Walleye ----

Hg_LMER_WE_MASS <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                          (1|WATERBODY_CODE_SAMPLE_YEAR), data = Hg_WE, 
                        control = lmerControl(optimizer = "Nelder_Mead")) 
summary(Hg_LMER_WE_MASS)
summary(lm(Hg_WE$VALUE_LOG ~ fitted(Hg_LMER_WE_MASS)))

hist(resid(Hg_LMER_WE_MASS), breaks = 100) # seems good 
plot(resid(Hg_LMER_WE_MASS) ~ fitted(Hg_LMER_WE_MASS)) # seems reasonable 

## Generate predictions for each waterbody based on the median size of WE in whole dataset 
Hg_LMER_WE_MASS_predict <- unique(Hg_WE[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
Hg_LMER_WE_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

## Bootstrap confidence intervals from bootMer
clust <- LMER_boot_initiate(varlist = "Hg_LMER_WE_MASS_predict")
system.time(Hg_LMER_WE_MASS_boot_est <- lme4::bootMer(Hg_LMER_WE_MASS, LMER_boot_est, 
                                                      nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                      parallel = "snow", 
                                                      cl = clust, 
                                                      ncpus = 4)) #17s for 100
LMER_boot_est_summary(Hg_LMER_WE_MASS_boot_est) ## 889s

Hg_LMER_WE_MASS_boot_pred <- function(., newdata) {
  predict(., newdata=Hg_LMER_WE_MASS_predict)
}

clust <- LMER_boot_initiate(varlist = "Hg_LMER_WE_MASS_predict")
system.time(Hg_LMER_WE_MASS_boot_pred_res <- lme4::bootMer(Hg_LMER_WE_MASS, Hg_LMER_WE_MASS_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 

(Hg_LMER_WE_MASS_boot_pred_res_sum <- LMER_boot_pred_summary(Hg_LMER_WE_MASS_boot_pred_res))

## ResuWEs 
LMER_boot_est_summary(Hg_LMER_WE_MASS_boot_est)
head(Hg_LMER_WE_MASS_boot_pred_res_sum)

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/glmmHg_LMER/vignettes/model_evaluation.pdf

## Easy to generate effects plots with glmmHg_LMER 
ae_MASS <- allEffects(Hg_LMER_WE_MASS)

par(mfrow = c(1,2))
plot(ae_MASS$WEIGHT_GRAM_LOG$data$VALUE_LOG ~ ae_MASS$WEIGHT_GRAM_LOG$data$WEIGHT_GRAM_LOG, main = "MASS")  
lines(ae_MASS$WEIGHT_GRAM_LOG$fit ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$upper ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, WEy = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$lower ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, WEy = 2)

(pp <- sjPlot::plot_model(Hg_LMER_WE_MASS,type=c("pred"),
                          terms=c("WEIGHT_GRAM_LOG","WATERBODY_CODE"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))

save.image(file = "./out_workspaces/HGAS_scratch.RData")

## 1.5) INLA - mixed model, Bayesian, fast ----

## Set up precision --> standard deviation formula; Bayesian models use precision (tau) where sd = 1/sqrt(tau) 
MySqrt <- function(x) {
  1 / sqrt(x)
}

## 1.x.x) Lake Trout ---- 

Hg_LT$WATERBODY_CODE <- factor(Hg_LT$WATERBODY_CODE) # for inla model
Hg_LT$WATERBODY_CODE1 <- as.integer(Hg_LT$WATERBODY_CODE) # for inla model
Hg_LT$WATERBODY_CODE2 <- Hg_LT$WATERBODY_CODE1 + max(Hg_LT$WATERBODY_CODE1) # for inla model
Hg_LT_n_waterbody <- dplyr::n_distinct(Hg_LT$WATERBODY_CODE) # for inla model

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
Hg_INLA_LT_MASS <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                       f(WATERBODY_CODE1, n = 2 * Hg_LT$n_waterbody, model = "iid2d") +   
                       f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                       f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                     data = Hg_LT, 
                     control.predictor = list(
                       compute = TRUE, 
                       quantiles = c(0.025, 0.5, 0.975)
                     ),
                     control.compute = list(
                       cpo = TRUE
                     )
)
summary(Hg_INLA_LT_MASS)

## There is a lot of information in the Hg_INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
Hg_INLA_LT_MASS_fits <- data.frame(WATERBODY_CODE = Hg_LT$WATERBODY_CODE, 
                                VALUE_LOG = Hg_LT$VALUE_LOG)

## Alternative comparison to fitted and residuals
Hg_INLA_LT_MASS_fits$inla_posterior_q50 <- Hg_INLA_LT_MASS$summary.fitted.values[, "0.5quant"]
Hg_INLA_LT_MASS_fits$inla_posterior_q2p5 <- Hg_INLA_LT_MASS$summary.fitted.values[, "0.025quant"]
Hg_INLA_LT_MASS_fits$inla_posterior_q97p5 <- Hg_INLA_LT_MASS$summary.fitted.values[, "0.975quant"]
Hg_INLA_LT_MASS_fits$resid_inla <- Hg_INLA_LT_MASS_fits$VALUE_LOG - Hg_INLA_LT_MASS_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = Hg_INLA_LT_MASS_fits)
j <- order(Hg_INLA_LT_MASS_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = Hg_INLA_LT_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ Hg_INLA_LT_MASS_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(Hg_INLA_LT_MASS_fits$VALUE_LOG ~ Hg_INLA_LT_MASS_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(Hg_INLA_LT_MASS_fits$VALUE_LOG ~ Hg_INLA_LT_MASS_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(Hg_INLA_LT_MASS_fits$resid_inla, breaks = 100) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
Hg_INLA_LT_MASS_fixed <- data.frame(
  ID = rownames(Hg_INLA_LT_MASS$summary.fixed),
  Hg_INLA_LT_MASS$summary.fixed, stringsAsFactors = FALSE
)
names(Hg_INLA_LT_MASS_fixed) <- c("ID", names(Hg_INLA_LT_MASS$summary.fixed))
Hg_INLA_LT_MASS_fixed$Type <- "Fixed"
head(Hg_INLA_LT_MASS_fixed)

summary(Hg_INLA_LT_MASS)

1/sqrt(Hg_INLA_LT_MASS$summary.hyperpar)

Hg_INLA_LT_MASS_random_intercept <- Hg_INLA_LT_MASS$summary.random$WATERBODY_CODE1
Hg_INLA_LT_MASS_random_intercept$ID <- as.character(Hg_INLA_LT_MASS_random_intercept$ID)
Hg_INLA_LT_MASS_random_intercept$WATERBODY_CODE1 <- Hg_INLA_LT_MASS_random_intercept$ID
Hg_INLA_LT_MASS_random_intercept <- merge(Hg_INLA_LT_MASS_random_intercept, LT[,c("WATERBODY_CODE1", "WATERBODY_CODE")], no.dups = TRUE)
Hg_INLA_LT_MASS_random_intercept <- Hg_INLA_LT_MASS_random_intercept[!duplicated(Hg_INLA_LT_MASS_random_intercept),]
Hg_INLA_LT_MASS_random_intercept$Type <- "Random Intercept - Waterbody"
head(Hg_INLA_LT_MASS_random_intercept)

Hg_INLA_LT_MASS_random_slope <- Hg_INLA_LT_MASS$summary.random$WATERBODY_CODE2
Hg_INLA_LT_MASS_random_slope$ID <- as.character(Hg_INLA_LT_MASS_random_slope$ID)
Hg_INLA_LT_MASS_random_slope$WATERBODY_CODE2 <- Hg_INLA_LT_MASS_random_slope$ID
Hg_INLA_LT_MASS_random_slope <- merge(Hg_INLA_LT_MASS_random_slope, LT[,c("WATERBODY_CODE2", "WATERBODY_CODE")])
Hg_INLA_LT_MASS_random_slope <- Hg_INLA_LT_MASS_random_slope[!duplicated(Hg_INLA_LT_MASS_random_slope),]
Hg_INLA_LT_MASS_random_slope$Type <- "Random Slope - Waterbody"
head(Hg_INLA_LT_MASS_random_slope)

Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_LT_MASS$summary.random$WATERBODY_CODE_SAMPLE_YEAR
Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$ID
Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR[!duplicated(Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR),]
Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$Type <- "Random Intercept - Waterbody Sampling Year "
head(Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR)

Hg_INLA_LT_MASS_summary <- list(
  Hg_INLA_LT_MASS_fixed,
  Hg_INLA_LT_MASS_random_intercept,
  Hg_INLA_LT_MASS_random_slope,
  Hg_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR
)

Hg_INLA_LT_MASS_summary <- dplyr::bind_rows(Hg_INLA_LT_MASS_summary)
head(Hg_INLA_LT_MASS_summary)
tail(Hg_INLA_LT_MASS_summary)

## Variance parameters, difficult to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for the Gaussian observations`
Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`

tau_WATERBODY_CODE1 <- Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`
tau_residual <- Hg_INLA_LT_MASS$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_WATERBODY_CODE_SAMPLE_YEAR <- inla.emarginal(MySqrt, tau_WATERBODY_CODE_SAMPLE_YEAR))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(Hg_LT_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)^2)

summary(Hg_INLA_LT_MASS)

## TMB vs. Hg_INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)
summary(Hg_TMB_LT_MASS)

## WATERBODY_CODE absolute variation >> WATERBODY_CODE VALUE_LOG~LENGTH_CM relation > WATERBODY_SAMPLE_YEAR > RESIDUAL 

## Now predictions 

## Predictions in Hg_INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the VALUE_LOG values are set to NA
## Then you extract the posterior information 

## Take TMB predictions but fit the Waterbody Codes to the dataframe 
head(Hg_LMER_LT_MASS_predict) 

Hg_INLA_LT_MASS_predict <- merge(Hg_LMER_LT_MASS_predict, Hg_LT[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
Hg_INLA_LT_MASS_predict <- Hg_INLA_LT_MASS_predict[!duplicated(Hg_INLA_LT_MASS_predict),]
Hg_INLA_LT_MASS_predict$VALUE_LOG <- NA

nrow(Hg_LMER_LT_MASS_predict)
nrow(Hg_INLA_LT_MASS_predict)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(Hg_TMB_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, Hg_INLA_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

Hg_INLA_LT_MASS_prediction <- Hg_LT[names(Hg_INLA_LT_MASS_predict)]

Hg_INLA_LT_MASS_prediction <- rbind(Hg_INLA_LT_MASS_prediction, 
                                 Hg_INLA_LT_MASS_predict)

head(Hg_INLA_LT_MASS_prediction)
tail(Hg_INLA_LT_MASS_prediction)

nrow(Hg_INLA_LT_MASS_prediction) - nrow(Hg_LT)

Hg_INLA_LT_MASS_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                        f(WATERBODY_CODE1, n = 2 * Hg_LT$n_waterbody, model = "iid2d") +   
                                        f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                        f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                      data = Hg_INLA_LT_MASS_prediction, 
                                      control.predictor = list(
                                        compute = TRUE, 
                                        quantiles = c(0.025, 0.5, 0.975)
                                      ),
                                      control.compute = list(
                                        cpo = TRUE
                                      )
)

Hg_INLA_LT_MASS_prediction_posteriors <- Hg_INLA_LT_MASS_prediction_model$summary.fitted.values[
  (nrow(Hg_LT) + 1):nrow(Hg_INLA_LT_MASS_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(Hg_INLA_LT_MASS_prediction_posteriors)

Hg_INLA_LT_MASS_prediction_posteriors <- exp(Hg_INLA_LT_MASS_prediction_posteriors)

Hg_INLA_LT_MASS_predict <- cbind(Hg_INLA_LT_MASS_predict, Hg_INLA_LT_MASS_prediction_posteriors)
Hg_INLA_LT_MASS_predict <- Hg_INLA_LT_MASS_predict[!duplicated(Hg_INLA_LT_MASS_predict),]

## Does this condition still hold true for comparison?
identical(Hg_TMB_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, Hg_INLA_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

## Produces very similar results but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(exp(Hg_LMER_LT_MASS_boot_pred_res_sum$fit)  ~ Hg_INLA_LT_MASS_predict$`0.5quant`)
abline(0,1)
head(Hg_TMB_LT_MASS_predict)
head(Hg_INLA_LT_MASS_predict)

## 1.x.x) Northern Pike ---- 

Hg_NP$WATERBODY_CODE <- factor(Hg_NP$WATERBODY_CODE) # for inla model
Hg_NP$WATERBODY_CODE1 <- as.integer(Hg_NP$WATERBODY_CODE) # for inla model
Hg_NP$WATERBODY_CODE2 <- Hg_NP$WATERBODY_CODE1 + max(Hg_NP$WATERBODY_CODE1) # for inla model
Hg_NP$n_waterbody <- dplyr::n_distinct(Hg_NP$WATERBODY_CODE) # for inla model

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
Hg_INLA_NP_MASS <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                          f(WATERBODY_CODE1, n = 2 * Hg_NP$n_waterbody, model = "iid2d") +   
                          f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                          f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                        data = Hg_NP, 
                        control.predictor = list(
                          compute = TRUE, 
                          quantiles = c(0.025, 0.5, 0.975)
                        ),
                        control.compute = list(
                          cpo = TRUE
                        )
)
summary(Hg_INLA_NP_MASS)

## There is a lot of information in the Hg_INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
Hg_INLA_NP_MASS_fits <- data.frame(WATERBODY_CODE = Hg_NP$WATERBODY_CODE, 
                                   VALUE_LOG = Hg_NP$VALUE_LOG)

## ANPernative comparison to fitted and residuals
Hg_INLA_NP_MASS_fits$inla_posterior_q50 <- Hg_INLA_NP_MASS$summary.fitted.values[, "0.5quant"]
Hg_INLA_NP_MASS_fits$inla_posterior_q2p5 <- Hg_INLA_NP_MASS$summary.fitted.values[, "0.025quant"]
Hg_INLA_NP_MASS_fits$inla_posterior_q97p5 <- Hg_INLA_NP_MASS$summary.fitted.values[, "0.975quant"]
Hg_INLA_NP_MASS_fits$resid_inla <- Hg_INLA_NP_MASS_fits$VALUE_LOG - Hg_INLA_NP_MASS_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = Hg_INLA_NP_MASS_fits)
j <- order(Hg_INLA_NP_MASS_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = Hg_INLA_NP_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ Hg_INLA_NP_MASS_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(Hg_INLA_NP_MASS_fits$VALUE_LOG ~ Hg_INLA_NP_MASS_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(Hg_INLA_NP_MASS_fits$VALUE_LOG ~ Hg_INLA_NP_MASS_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(Hg_INLA_NP_MASS_fits$resid_inla, breaks = 100) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
Hg_INLA_NP_MASS_fixed <- data.frame(
  ID = rownames(Hg_INLA_NP_MASS$summary.fixed),
  Hg_INLA_NP_MASS$summary.fixed, stringsAsFactors = FALSE
)
names(Hg_INLA_NP_MASS_fixed) <- c("ID", names(Hg_INLA_NP_MASS$summary.fixed))
Hg_INLA_NP_MASS_fixed$Type <- "Fixed"
head(Hg_INLA_NP_MASS_fixed)

summary(Hg_INLA_NP_MASS)

1/sqrt(Hg_INLA_NP_MASS$summary.hyperpar)

Hg_INLA_NP_MASS_random_intercept <- Hg_INLA_NP_MASS$summary.random$WATERBODY_CODE1
Hg_INLA_NP_MASS_random_intercept$ID <- as.character(Hg_INLA_NP_MASS_random_intercept$ID)
Hg_INLA_NP_MASS_random_intercept$WATERBODY_CODE1 <- Hg_INLA_NP_MASS_random_intercept$ID
Hg_INLA_NP_MASS_random_intercept <- merge(Hg_INLA_NP_MASS_random_intercept, NP[,c("WATERBODY_CODE1", "WATERBODY_CODE")], no.dups = TRUE)
Hg_INLA_NP_MASS_random_intercept <- Hg_INLA_NP_MASS_random_intercept[!duplicated(Hg_INLA_NP_MASS_random_intercept),]
Hg_INLA_NP_MASS_random_intercept$Type <- "Random Intercept - Waterbody"
head(Hg_INLA_NP_MASS_random_intercept)

Hg_INLA_NP_MASS_random_slope <- Hg_INLA_NP_MASS$summary.random$WATERBODY_CODE2
Hg_INLA_NP_MASS_random_slope$ID <- as.character(Hg_INLA_NP_MASS_random_slope$ID)
Hg_INLA_NP_MASS_random_slope$WATERBODY_CODE2 <- Hg_INLA_NP_MASS_random_slope$ID
Hg_INLA_NP_MASS_random_slope <- merge(Hg_INLA_NP_MASS_random_slope, NP[,c("WATERBODY_CODE2", "WATERBODY_CODE")])
Hg_INLA_NP_MASS_random_slope <- Hg_INLA_NP_MASS_random_slope[!duplicated(Hg_INLA_NP_MASS_random_slope),]
Hg_INLA_NP_MASS_random_slope$Type <- "Random Slope - Waterbody"
head(Hg_INLA_NP_MASS_random_slope)

Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_NP_MASS$summary.random$WATERBODY_CODE_SAMPLE_YEAR
Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$ID
Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR[!duplicated(Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR),]
Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$Type <- "Random Intercept - Waterbody Sampling Year "
head(Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR)

Hg_INLA_NP_MASS_summary <- list(
  Hg_INLA_NP_MASS_fixed,
  Hg_INLA_NP_MASS_random_intercept,
  Hg_INLA_NP_MASS_random_slope,
  Hg_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR
)

Hg_INLA_NP_MASS_summary <- dplyr::bind_rows(Hg_INLA_NP_MASS_summary)
head(Hg_INLA_NP_MASS_summary)
tail(Hg_INLA_NP_MASS_summary)

## Variance parameters, difficuNP to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for the Gaussian observations`
Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`

tau_WATERBODY_CODE1 <- Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`
tau_residual <- Hg_INLA_NP_MASS$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_WATERBODY_CODE_SAMPLE_YEAR <- inla.emarginal(MySqrt, tau_WATERBODY_CODE_SAMPLE_YEAR))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(Hg_NP_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)^2)

## TMB vs. Hg_INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)

## WATERBODY_CODE absolute variation >> WATERBODY_CODE VALUE_LOG~LENGTH_CM relation > WATERBODY_SAMPLE_YEAR > RESIDUAL 

## Now predictions 

## Predictions in Hg_INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the VALUE_LOG values are set to NA
## Then you extract the posterior information 

## Take TMB predictions but fit the Waterbody Codes to the dataframe 
head(Hg_LMER_NP_MASS_predict) 

Hg_INLA_NP_MASS_predict <- merge(Hg_LMER_NP_MASS_predict, Hg_NP[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
Hg_INLA_NP_MASS_predict <- Hg_INLA_NP_MASS_predict[!duplicated(Hg_INLA_NP_MASS_predict),]
Hg_INLA_NP_MASS_predict$VALUE_LOG <- NA

nrow(Hg_LMER_NP_MASS_predict)
nrow(Hg_INLA_NP_MASS_predict)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(Hg_TMB_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, Hg_INLA_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

Hg_INLA_NP_MASS_prediction <- Hg_NP[names(Hg_INLA_NP_MASS_predict)]

Hg_INLA_NP_MASS_prediction <- rbind(Hg_INLA_NP_MASS_prediction, 
                                    Hg_INLA_NP_MASS_predict)

head(Hg_INLA_NP_MASS_prediction)
tail(Hg_INLA_NP_MASS_prediction)

nrow(Hg_INLA_NP_MASS_prediction) - nrow(Hg_NP)

Hg_INLA_NP_MASS_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * Hg_NP$n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = Hg_INLA_NP_MASS_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
)

log(1000)

Hg_INLA_NP_MASS_prediction_posteriors <- Hg_INLA_NP_MASS_prediction_model$summary.fitted.values[
  (nrow(Hg_NP) + 1):nrow(Hg_INLA_NP_MASS_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(Hg_INLA_NP_MASS_prediction_posteriors)

Hg_INLA_NP_MASS_prediction_posteriors <- exp(Hg_INLA_NP_MASS_prediction_posteriors)

Hg_INLA_NP_MASS_predict <- cbind(Hg_INLA_NP_MASS_predict, Hg_INLA_NP_MASS_prediction_posteriors)
Hg_INLA_NP_MASS_predict <- Hg_INLA_NP_MASS_predict[!duplicated(Hg_INLA_NP_MASS_predict),]

## Does this condition still hold true for comparison?
identical(Hg_TMB_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, Hg_INLA_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

## Produces very similar resuNPs but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(exp(Hg_LMER_NP_MASS_boot_pred_res_sum$fit)  ~ Hg_INLA_NP_MASS_predict$`0.5quant`)
abline(0,1)
head(Hg_TMB_NP_MASS_predict)
head(Hg_INLA_NP_MASS_predict)


## 1.x.x) Walleye ---- 

Hg_WE$WATERBODY_CODE <- factor(Hg_WE$WATERBODY_CODE) # for inla model
Hg_WE$WATERBODY_CODE1 <- as.integer(Hg_WE$WATERBODY_CODE) # for inla model
Hg_WE$WATERBODY_CODE2 <- Hg_WE$WATERBODY_CODE1 + max(Hg_WE$WATERBODY_CODE1) # for inla model
Hg_WE$n_waterbody <- dplyr::n_distinct(Hg_WE$WATERBODY_CODE) # for inla model

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
Hg_INLA_WE_MASS <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                          f(WATERBODY_CODE1, n = 2 * Hg_WE$n_waterbody, model = "iid2d") +   
                          f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                          f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                        data = Hg_WE, 
                        control.predictor = list(
                          compute = TRUE, 
                          quantiles = c(0.025, 0.5, 0.975)
                        ),
                        control.compute = list(
                          cpo = TRUE
                        )
)
summary(Hg_INLA_WE_MASS)

## There is a lot of information in the Hg_INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
Hg_INLA_WE_MASS_fits <- data.frame(WATERBODY_CODE = Hg_WE$WATERBODY_CODE, 
                                   VALUE_LOG = Hg_WE$VALUE_LOG)

## AWEernative comparison to fitted and residuals
Hg_INLA_WE_MASS_fits$inla_posterior_q50 <- Hg_INLA_WE_MASS$summary.fitted.values[, "0.5quant"]
Hg_INLA_WE_MASS_fits$inla_posterior_q2p5 <- Hg_INLA_WE_MASS$summary.fitted.values[, "0.025quant"]
Hg_INLA_WE_MASS_fits$inla_posterior_q97p5 <- Hg_INLA_WE_MASS$summary.fitted.values[, "0.975quant"]
Hg_INLA_WE_MASS_fits$resid_inla <- Hg_INLA_WE_MASS_fits$VALUE_LOG - Hg_INLA_WE_MASS_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = Hg_INLA_WE_MASS_fits)
j <- order(Hg_INLA_WE_MASS_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = Hg_INLA_WE_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ Hg_INLA_WE_MASS_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(Hg_INLA_WE_MASS_fits$VALUE_LOG ~ Hg_INLA_WE_MASS_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(Hg_INLA_WE_MASS_fits$VALUE_LOG ~ Hg_INLA_WE_MASS_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(Hg_INLA_WE_MASS_fits$resid_inla, breaks = 100) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
Hg_INLA_WE_MASS_fixed <- data.frame(
  ID = rownames(Hg_INLA_WE_MASS$summary.fixed),
  Hg_INLA_WE_MASS$summary.fixed, stringsAsFactors = FALSE
)
names(Hg_INLA_WE_MASS_fixed) <- c("ID", names(Hg_INLA_WE_MASS$summary.fixed))
Hg_INLA_WE_MASS_fixed$Type <- "Fixed"
head(Hg_INLA_WE_MASS_fixed)

summary(Hg_INLA_WE_MASS)

1/sqrt(Hg_INLA_WE_MASS$summary.hyperpar)

Hg_INLA_WE_MASS_random_intercept <- Hg_INLA_WE_MASS$summary.random$WATERBODY_CODE1
Hg_INLA_WE_MASS_random_intercept$ID <- as.character(Hg_INLA_WE_MASS_random_intercept$ID)
Hg_INLA_WE_MASS_random_intercept$WATERBODY_CODE1 <- Hg_INLA_WE_MASS_random_intercept$ID
Hg_INLA_WE_MASS_random_intercept <- merge(Hg_INLA_WE_MASS_random_intercept, WE[,c("WATERBODY_CODE1", "WATERBODY_CODE")], no.dups = TRUE)
Hg_INLA_WE_MASS_random_intercept <- Hg_INLA_WE_MASS_random_intercept[!duplicated(Hg_INLA_WE_MASS_random_intercept),]
Hg_INLA_WE_MASS_random_intercept$Type <- "Random Intercept - Waterbody"
head(Hg_INLA_WE_MASS_random_intercept)

Hg_INLA_WE_MASS_random_slope <- Hg_INLA_WE_MASS$summary.random$WATERBODY_CODE2
Hg_INLA_WE_MASS_random_slope$ID <- as.character(Hg_INLA_WE_MASS_random_slope$ID)
Hg_INLA_WE_MASS_random_slope$WATERBODY_CODE2 <- Hg_INLA_WE_MASS_random_slope$ID
Hg_INLA_WE_MASS_random_slope <- merge(Hg_INLA_WE_MASS_random_slope, WE[,c("WATERBODY_CODE2", "WATERBODY_CODE")])
Hg_INLA_WE_MASS_random_slope <- Hg_INLA_WE_MASS_random_slope[!duplicated(Hg_INLA_WE_MASS_random_slope),]
Hg_INLA_WE_MASS_random_slope$Type <- "Random Slope - Waterbody"
head(Hg_INLA_WE_MASS_random_slope)

Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_WE_MASS$summary.random$WATERBODY_CODE_SAMPLE_YEAR
Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$ID
Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR[!duplicated(Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR),]
Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$Type <- "Random Intercept - Waterbody Sampling Year "
head(Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR)

Hg_INLA_WE_MASS_summary <- list(
  Hg_INLA_WE_MASS_fixed,
  Hg_INLA_WE_MASS_random_intercept,
  Hg_INLA_WE_MASS_random_slope,
  Hg_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR
)

Hg_INLA_WE_MASS_summary <- dplyr::bind_rows(Hg_INLA_WE_MASS_summary)
head(Hg_INLA_WE_MASS_summary)
tail(Hg_INLA_WE_MASS_summary)

## Variance parameters, difficuWE to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for the Gaussian observations`
Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`

tau_WATERBODY_CODE1 <- Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_WATERBODY_CODE_SAMPLE_YEAR <- Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`
tau_residual <- Hg_INLA_WE_MASS$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_WATERBODY_CODE_SAMPLE_YEAR <- inla.emarginal(MySqrt, tau_WATERBODY_CODE_SAMPLE_YEAR))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(Hg_WE_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)^2)

## TMB vs. Hg_INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)
summary(Hg_TMB_WE_MASS)

## WATERBODY_CODE absolute variation >> WATERBODY_CODE VALUE_LOG~LENGTH_CM relation > WATERBODY_SAMPLE_YEAR > RESIDUAL 

## Now predictions 

## Predictions in Hg_INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the VALUE_LOG values are set to NA
## Then you extract the posterior information 

## Take TMB predictions but fit the Waterbody Codes to the dataframe 
head(Hg_LMER_WE_MASS_predict) 

Hg_INLA_WE_MASS_predict <- merge(Hg_LMER_WE_MASS_predict, Hg_WE[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
Hg_INLA_WE_MASS_predict <- Hg_INLA_WE_MASS_predict[!duplicated(Hg_INLA_WE_MASS_predict),]
Hg_INLA_WE_MASS_predict$VALUE_LOG <- NA

nrow(Hg_LMER_WE_MASS_predict)
nrow(Hg_INLA_WE_MASS_predict)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(Hg_TMB_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, Hg_INLA_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

Hg_INLA_WE_MASS_prediction <- Hg_WE[names(Hg_INLA_WE_MASS_predict)]

Hg_INLA_WE_MASS_prediction <- rbind(Hg_INLA_WE_MASS_prediction, 
                                    Hg_INLA_WE_MASS_predict)

head(Hg_INLA_WE_MASS_prediction)
tail(Hg_INLA_WE_MASS_prediction)

nrow(Hg_INLA_WE_MASS_prediction) - nrow(Hg_WE)

Hg_INLA_WE_MASS_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * Hg_WE$n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = Hg_INLA_WE_MASS_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
)

Hg_INLA_WE_MASS_prediction_posteriors <- Hg_INLA_WE_MASS_prediction_model$summary.fitted.values[
  (nrow(Hg_WE) + 1):nrow(Hg_INLA_WE_MASS_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(Hg_INLA_WE_MASS_prediction_posteriors)

Hg_INLA_WE_MASS_prediction_posteriors <- exp(Hg_INLA_WE_MASS_prediction_posteriors)

Hg_INLA_WE_MASS_predict <- cbind(Hg_INLA_WE_MASS_predict, Hg_INLA_WE_MASS_prediction_posteriors)
Hg_INLA_WE_MASS_predict <- Hg_INLA_WE_MASS_predict[!duplicated(Hg_INLA_WE_MASS_predict),]

## Does this condition still hold true for comparison?
identical(Hg_TMB_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, Hg_INLA_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

## Produces very similar resuWEs but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(exp(Hg_LMER_WE_MASS_boot_pred_res_sum$fit)  ~ Hg_INLA_WE_MASS_predict$`0.5quant`)
abline(0,1)
head(Hg_TMB_WE_MASS_predict)
head(Hg_INLA_WE_MASS_predict)

## 1.6) RSTAN - mixed model, Bayesian slow ---- 

## 1.x.x) Lake Trout ---- 

system.time(Hg_STAN_LT_MASS <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                             (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                             (1|WATERBODY_CODE_SAMPLE_YEAR), data = LT, 
                           cores = 4, chains = 4))
saveRDS(Hg_STAN_LT_MASS, "./out_workspaces/Hg_STAN_LT_MASS.rds")
Hg_STAN_LT_MASS <- readRDS("./out_workspaces/STAN_LT_MASS.rds")

# shinystan
launch_shinystan(Hg_STAN_LT_MASS)

# Marginal r-squared (no random effects)
Hg_STAN_LT_MASS_rsq_marg <- bayes_R2(Hg_STAN_LT_MASS, re.form = NA)
median(Hg_STAN_LT_MASS_rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
Hg_STAN_LT_MASS_rsq_cond_lake <- bayes_R2(Hg_STAN_LT_MASS, re.form = ~ (WEIGHT_GRAM_LOG|WATERBODY_CODE) )
median(Hg_STAN_LT_MASS_rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
Hg_STAN_LT_MASS_fits <- data.frame(
  VALUE_LOG = Hg_LT$VALUE_LOG,
  Hg_STAN_FIT = fitted(Hg_STAN_LT_MASS),
  Hg_STAN_RESID = residuals(Hg_STAN_LT_MASS)
)

lw1_stan <- loess(Hg_STAN_RESID ~ Hg_STAN_FIT, data = Hg_STAN_LT_MASS_fits)
j <- order(Hg_STAN_LT_MASS_fits$Hg_STAN_FIT)
plot(Hg_STAN_RESID ~ Hg_STAN_FIT, data = Hg_STAN_LT_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ Hg_STAN_LT_MASS_fits$Hg_STAN_FIT[j], lwd = 2)

## Histogram of residuals
hist(Hg_STAN_LT_MASS_fits$Hg_STAN_RESID, breaks = 100)

## Alternative comparison to fitted and residuals
Hg_STAN_LT_MASS_posterior_full <- rstanarm::posterior_predict(Hg_STAN_LT_MASS)
Hg_STAN_LT_MASS_posterior_est_full <- apply(Hg_STAN_LT_MASS_posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

Hg_STAN_LT_MASS_fits$stan_posterior_q50 <- Hg_STAN_LT_MASS_posterior_est_full[2, ]
Hg_STAN_LT_MASS_fits$stan_posterior_q2p5 <- Hg_STAN_LT_MASS_posterior_est_full[1, ]
Hg_STAN_LT_MASS_fits$stan_posterior_q97p5 <- Hg_STAN_LT_MASS_posterior_est_full[3, ]

## Parameter estimates from posterior distribution
print(Hg_STAN_LT_MASS, detail = TRUE)

summary(Hg_STAN_LT_MASS, pars = "WATERBODY_SAMPLE_CODE")

Hg_STAN_LT_MASS_sims <- as.data.frame(as.matrix(Hg_STAN_LT_MASS))

head(Hg_STAN_LT_MASS_sims[,grep("Sigma", names(Hg_STAN_LT_MASS_sims))])

a_quant <- apply(
  X = Hg_STAN_LT_MASS_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

Hg_STAN_LT_MASS_summary <- a_quant

## Variance parameters
Hg_STAN_LT_RE <- Hg_STAN_LT_MASS_sims[,grep("Sigma", names(Hg_STAN_LT_MASS_sims), ignore.case = TRUE)]

Hg_STAN_LT_RE <- apply(
  X = Hg_STAN_LT_RE,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.50, 0.025, 0.975)
)

sqrt(Hg_STAN_LT_RE)

## Predictions 

Hg_STAN_LT_MASS_predict <- unique(Hg_LT[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
Hg_STAN_LT_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

Hg_STAN_LT_MASS_predictions <- rstanarm::posterior_predict(Hg_STAN_LT_MASS,
                                             newdata = Hg_STAN_LT_MASS_predict
)

Hg_STAN_LT_MASS_predictions <- apply(Hg_STAN_LT_MASS_predictions,
                           MARGIN = 2,
                           function(x) {
                             quantile(x, probs = c(0.025, 0.5, 0.975))
                           }
)

Hg_STAN_LT_MASS_predictions <- exp(t(Hg_STAN_LT_MASS_predictions))
colnames(Hg_STAN_LT_MASS_predictions) <- c("Hg_STAN_0.025quant", "Hg_STAN_0.5quant", "Hg_STAN_0.975quant")
Hg_STAN_LT_MASS_predictions <- cbind(Hg_STAN_LT_MASS_predict, Hg_STAN_LT_MASS_predictions)

par(mfrow = c(1,2))
plot(exp(Hg_LMER_LT_MASS_boot_pred_res_sum$fit), Hg_STAN_LT_MASS_predictions$Hg_STAN_0.5quant, xlab = "Maximum Likelihood Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))
plot((Hg_INLA_LT_MASS_prediction_posteriors$`0.5quant`), Hg_STAN_LT_MASS_predictions$Hg_STAN_0.5quant,  xlab = "Bayesian INLA Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))

## 1.x.x) Northern Pike ---- 

system.time(Hg_STAN_NP_MASS <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                             (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                             (1|WATERBODY_CODE_SAMPLE_YEAR), data = NP,
                           cores = 4, chains = 4))
saveRDS(Hg_STAN_NP_MASS, "./out_workspaces/Hg_STAN_NP_MASS.rds")
Hg_STAN_NP_MASS <- readRDS( "./out_workspaces/Hg_STAN_NP_MASS.rds")

# shinystan
launch_shinystan(Hg_STAN_NP_MASS)

# Marginal r-squared (no random effects)
Hg_STAN_NP_MASS_rsq_marg <- bayes_R2(Hg_STAN_NP_MASS, re.form = NA)
median(Hg_STAN_NP_MASS_rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
Hg_STAN_NP_MASS_rsq_cond_lake <- bayes_R2(Hg_STAN_NP_MASS, re.form = ~ (WEIGHT_GRAM_LOG|WATERBODY_CODE) )
median(Hg_STAN_NP_MASS_rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
Hg_STAN_NP_MASS_fits <- data.frame(
  VALUE_LOG = Hg_NP$VALUE_LOG,
  Hg_STAN_FIT = fitted(Hg_STAN_NP_MASS),
  Hg_STAN_RESID = residuals(Hg_STAN_NP_MASS)
)

lw1_stan <- loess(Hg_STAN_RESID ~ Hg_STAN_FIT, data = Hg_STAN_NP_MASS_fits)
j <- order(Hg_STAN_NP_MASS_fits$Hg_STAN_FIT)
plot(Hg_STAN_RESID ~ Hg_STAN_FIT, data = Hg_STAN_NP_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ Hg_STAN_NP_MASS_fits$Hg_STAN_FIT[j], lwd = 2)

## Histogram of residuals
hist(Hg_STAN_NP_MASS_fits$Hg_STAN_RESID, breaks = 100)

## Alternative comparison to fitted and residuals
Hg_STAN_NP_MASS_posterior_full <- rstanarm::posterior_predict(Hg_STAN_NP_MASS)
Hg_STAN_NP_MASS_posterior_est_full <- apply(Hg_STAN_NP_MASS_posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

Hg_STAN_NP_MASS_fits$stan_posterior_q50 <- Hg_STAN_NP_MASS_posterior_est_full[2, ]
Hg_STAN_NP_MASS_fits$stan_posterior_q2p5 <- Hg_STAN_NP_MASS_posterior_est_full[1, ]
Hg_STAN_NP_MASS_fits$stan_posterior_q97p5 <- Hg_STAN_NP_MASS_posterior_est_full[3, ]

## Parameter estimates from posterior distribution
Hg_STAN_NP_MASS_sims <- as.data.frame(as.matrix(Hg_STAN_NP_MASS))
names(Hg_STAN_NP_MASS_sims)

a_quant <- apply(
  X = Hg_STAN_NP_MASS_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

Hg_STAN_NP_MASS_summary <- a_quant

## Variance parameters
Hg_STAN_NP_RE <- Hg_STAN_NP_MASS_sims[,grep("Sigma", names(Hg_STAN_NP_MASS_sims), ignore.case = TRUE)]

Hg_STAN_NP_RE <- apply(
  X = Hg_STAN_NP_RE,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.50, 0.025, 0.975)
)

sqrt(Hg_STAN_NP_RE)

## Predictions 
Hg_STAN_NP_MASS_predict <- unique(Hg_NP[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
Hg_STAN_NP_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)


Hg_STAN_NP_MASS_predictions <- rstanarm::posterior_predict(Hg_STAN_NP_MASS,
                                                        newdata = Hg_STAN_NP_MASS_predict
)

Hg_STAN_NP_MASS_predictions <- apply(Hg_STAN_NP_MASS_predictions,
                                  MARGIN = 2,
                                  function(x) {
                                    quantile(x, probs = c(0.025, 0.5, 0.975))
                                  }
)

Hg_STAN_NP_MASS_predictions <- exp(t(Hg_STAN_NP_MASS_predictions))
colnames(Hg_STAN_NP_MASS_predictions) <- c("Hg_STAN_0.025quant", "Hg_STAN_0.5quant", "Hg_STAN_0.975quant")
Hg_STAN_NP_MASS_predictions <- cbind(Hg_STAN_NP_MASS_predict, Hg_STAN_NP_MASS_predictions)

par(mfrow = c(1,2))
plot(exp(Hg_LMER_NP_MASS_boot_pred_res_sum$fit), Hg_STAN_NP_MASS_predictions$Hg_STAN_0.5quant, xlab = "Maximum Likelihood Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))
plot((Hg_INLA_NP_MASS_prediction_posteriors$`0.5quant`), Hg_STAN_NP_MASS_predictions$Hg_STAN_0.5quant,  xlab = "Bayesian INLA Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))

## 1.x.x) Walleye ----

system.time(Hg_STAN_WE_MASS <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                             (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                             (1|WATERBODY_CODE_SAMPLE_YEAR), data = WE,
                           cores = 4, chains = 4))
saveRDS(Hg_STAN_WE_MASS, "./out_workspaces/Hg_STAN_WE_MASS.rds")
Hg_STAN_WE_MASS <- readRDS("./out_workspaces/Hg_STAN_WE_MASS.rds")

# shinystan
launch_shinystan(Hg_STAN_WE_MASS)

# Marginal r-squared (no random effects)
Hg_STAN_WE_MASS_rsq_marg <- bayes_R2(Hg_STAN_WE_MASS, re.form = NA)
median(Hg_STAN_WE_MASS_rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
Hg_STAN_WE_MASS_rsq_cond_lake <- bayes_R2(Hg_STAN_WE_MASS, re.form = ~ (WEIGHT_GRAM_LOG|WATERBODY_CODE) )
median(Hg_STAN_WE_MASS_rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
Hg_STAN_WE_MASS_fits <- data.frame(
  VALUE_LOG = Hg_WE$VALUE_LOG,
  Hg_STAN_FIT = fitted(Hg_STAN_WE_MASS),
  Hg_STAN_RESID = residuals(Hg_STAN_WE_MASS)
)

lw1_stan <- loess(Hg_STAN_RESID ~ Hg_STAN_FIT, data = Hg_STAN_WE_MASS_fits)
j <- order(Hg_STAN_WE_MASS_fits$Hg_STAN_FIT)
plot(Hg_STAN_RESID ~ Hg_STAN_FIT, data = Hg_STAN_WE_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ Hg_STAN_WE_MASS_fits$Hg_STAN_FIT[j], lwd = 2)

## Histogram of residuals
hist(Hg_STAN_WE_MASS_fits$Hg_STAN_RESID, breaks = 100)

## Alternative comparison to fitted and residuals
Hg_STAN_WE_MASS_posterior_full <- rstanarm::posterior_predict(Hg_STAN_WE_MASS)
Hg_STAN_WE_MASS_posterior_est_full <- apply(Hg_STAN_WE_MASS_posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

Hg_STAN_WE_MASS_fits$stan_posterior_q50 <- Hg_STAN_WE_MASS_posterior_est_full[2, ]
Hg_STAN_WE_MASS_fits$stan_posterior_q2p5 <- Hg_STAN_WE_MASS_posterior_est_full[1, ]
Hg_STAN_WE_MASS_fits$stan_posterior_q97p5 <- Hg_STAN_WE_MASS_posterior_est_full[3, ]

## Parameter estimates from posterior distribution
Hg_STAN_WE_MASS_sims <- as.data.frame(as.matrix(Hg_STAN_WE_MASS))
names(Hg_STAN_WE_MASS_sims)

a_quant <- apply(
  X = Hg_STAN_WE_MASS_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

Hg_STAN_WE_MASS_summary <- a_quant

## Variance parameters
Hg_STAN_WE_RE <- Hg_STAN_WE_MASS_sims[,grep("Sigma", names(Hg_STAN_WE_MASS_sims), ignore.case = TRUE)]

Hg_STAN_WE_RE <- apply(
  X = Hg_STAN_WE_RE,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.50, 0.025, 0.975)
)

sqrt(Hg_STAN_WE_RE)

## Predictions 
Hg_STAN_WE_MASS_predict <- unique(Hg_WE[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
Hg_STAN_WE_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)


Hg_STAN_WE_MASS_predictions <- rstanarm::posterior_predict(Hg_STAN_WE_MASS,
                                                        newdata = Hg_STAN_WE_MASS_predict
)

Hg_STAN_WE_MASS_predictions <- apply(Hg_STAN_WE_MASS_predictions,
                                  MARGIN = 2,
                                  function(x) {
                                    quantile(x, probs = c(0.025, 0.5, 0.975))
                                  }
)

Hg_STAN_WE_MASS_predictions <- exp(t(Hg_STAN_WE_MASS_predictions))
colnames(Hg_STAN_WE_MASS_predictions) <- c("Hg_STAN_0.025quant", "Hg_STAN_0.5quant", "Hg_STAN_0.975quant")
Hg_STAN_WE_MASS_predictions <- cbind(Hg_STAN_WE_MASS_predict, Hg_STAN_WE_MASS_predictions)

par(mfrow = c(1,2))
plot(exp(Hg_LMER_WE_MASS_boot_pred_res_sum$fit), Hg_STAN_WE_MASS_predictions$Hg_STAN_0.5quant, xlab = "Maximum Likelihood Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))
plot((Hg_INLA_WE_MASS_prediction_posteriors$`0.5quant`), Hg_STAN_WE_MASS_predictions$Hg_STAN_0.5quant,  xlab = "Bayesian INLA Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))

## 1.7) Model comparison statistics ----

## 2) AS ----

## Read the saved file
As <- read.csv("./data/Fish_As_2021.12.01.csv")
As_sub <- subset(As, System_Type == "Lake")
As_sub <- subset(As_sub, Data_source == "MECP")

## Subset by "SKINLESS, BONELESS FILLET" OR "SKINLESS, BONELESS FILLET PLUG"
As_sub <- subset(As_sub, PORTION_TYPE_DESC %in% 
                   c("SKINLESS, BONELESS FILLET (STANDARD MOE DORSAL FILLET)", 
                     "SKINLESS, BONELESS FILLET PLUG (SKIN-OFF)"))

nrow(unique(As_sub[,c("Waterbody", "System_Type")])) #125

## Some data cleanup and variable creation to align with HG analysis
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
nrow(LT); length(with(LT, unique(WATERBODY_CODE))); 
median(table(LT$WATERBODY_CODE_SAMPLE_YEAR))
min(table(LT$WATERBODY_CODE_SAMPLE_YEAR))
max(table(LT$WATERBODY_CODE_SAMPLE_YEAR))

nrow(NP); length(with(NP, unique(WATERBODY_CODE))); 
median(table(NP$WATERBODY_CODE_SAMPLE_YEAR))
min(table(NP$WATERBODY_CODE_SAMPLE_YEAR))
max(table(NP$WATERBODY_CODE_SAMPLE_YEAR))

nrow(WE); length(with(WE, unique(WATERBODY_CODE))); 
median(table(WE$WATERBODY_CODE_SAMPLE_YEAR))
min(table(WE$WATERBODY_CODE_SAMPLE_YEAR))
max(table(WE$WATERBODY_CODE_SAMPLE_YEAR))

nrow(NP)
nrow(WE)

lk_years <- unique(As_sub[,c("WATERBODY_CODE", "SAMPLE_YEAR")])
(lk_years <- table((table(lk_years$WATERBODY_CODE))))

## 2.3) SER - sampling event regressions ----

## Standardize by mass, , >4 individuals needed and using Rob Mackereth's standardized values
x = As_sub$WATERBODY_CODE[1]
y = 2014
z = "WALL"

As_sub_preds_mass <- lapply(unique(As_sub$WATERBODY_CODE), function(x) {
  
  ## Subsets the waterbody of interest
  sub_waterbody <- subset(As_sub, WATERBODY_CODE == x)
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$SAMPLE_YEAR), function(y) {
    
    ## Subsets a year within a given waterbody
    sub_year <- subset(sub_waterbody, SAMPLE_YEAR == y)
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$SPECIES_CODE), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, SPECIES_CODE == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      ## the sum evaluation statement is meant to make sure that all
      ## masses are not NAs which was a case in these data
      if (nrow(sub_species) > 4 & !sum(sub_species$WEIGHT_GRAM, na.rm = T) == 0) {
        
        ## log mercury concentrations and masss (converted to mm)
        (logm <- log(sub_species$VALUE))
        (logl <- log(sub_species$WEIGHT_GRAM))
        
        percentile <- ecdf(sub_species$WEIGHT_GRAM)
        
        ## regression relationship
        rel <- lm(logm ~ logl)
        
        #plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        #abline(rel)
        
        ## Generating context and regression summary statistics
        (WATERBODY_CODE <- x) # waterbody
        (SAMPLE_YEAR <- y) # year
        (SPECIES_CODE <- z) # species
        (n <- length(logl)) # number of individuals
        (int <- formatC(rel$coefficients[1], digits = 4, format = "f")) # estimated intercept
        (slp <- formatC(rel$coefficients[2], digits = 4, format = "f")) # estimated slope
        (int_confint <- paste(formatC(confint(rel)[1, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated intercept confidence interval
        (slp_confint <- paste(formatC(confint(rel)[2, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated slope confidence interval
        (r2 <- formatC(summary(rel)$r.squared, digits = 4, format = "f")) # adjusted R2
        (RMSE <- sqrt(mean(rel$residuals^2)))
        
        ## Modify the predicted standardized length based on species
        (pred_modify <- 1000)
        
        (target_size_percentile <- percentile(pred_modify))
        
        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
                             newdata = data.frame(logl = log(pred_modify)),
                             interval = "confidence"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(WATERBODY_CODE, SAMPLE_YEAR, SPECIES_CODE, n,
                             int, int_confint, slp, slp_confint, r2, RMSE,
                             target_size_percentile, pred,
                             stringsAsFactors = FALSE
        ))
      } else {
        frame <- NULL
      }
      return(frame)
    })
    
    ## rbind species-specific regression relationships and predictions
    sub_species_frame <- do.call(rbind, sub_species_output)
    return(sub_species_frame)
  })
  ## rbind years for a waterbody
  sub_year_frame <- do.call(rbind, sub_year_output)
  sub_year_frame
  return(sub_year_frame)
})
As_sub_preds_mass <- As_sub_preds_mass[!sapply(As_sub_preds_mass, is.null)]
As_sub_preds_mass <- dplyr::bind_rows(As_sub_preds_mass)

## 2.4) lmer - mixed model, maximum likelihood ----

## 1.x.x) Lake Trout ----

As_LMER_LT_MASS <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_LT) 
summary(As_LMER_LT_MASS)

hist(resid(As_LMER_LT_MASS), breaks = 100) # seems good 
plot(resid(As_LMER_LT_MASS) ~ fitted(As_LMER_LT_MASS)) # seems reasonable; a few very high residuals 

## Generate predictions for each waterbody based on the median size of LT in whole dataset 

As_LMER_LT_MASS_predict <- unique(As_LT[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
As_LMER_LT_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

## Bootstrap confidence intervals from bootMer
clust <- LMER_boot_initiate(varlist = "As_LMER_LT_MASS_predict", nclust = 4)
system.time(As_LMER_LT_MASS_boot_est <- lme4::bootMer(As_LMER_LT_MASS, LMER_boot_est, 
                                                      nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                      parallel = "snow", 
                                                      cl = clust, 
                                                      ncpus = 4)) #17s for 100
LMER_boot_est_summary(As_LMER_LT_MASS_boot_est) ## 889s

As_LMER_LT_MASS_boot_pred <- function(., newdata) {
  predict(., newdata=As_LMER_LT_MASS_predict)
}

clust <- LMER_boot_initiate(varlist = "As_LMER_LT_MASS_predict", nclust = 4)
system.time(As_LMER_LT_MASS_boot_pred_res <- lme4::bootMer(As_LMER_LT_MASS, As_LMER_LT_MASS_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 

(As_LMER_LT_MASS_boot_pred_res_sum <- LMER_boot_pred_summary(As_LMER_LT_MASS_boot_pred_res))

## Results 
LMER_boot_est_summary(As_LMER_LT_MASS_boot_est)
head(As_LMER_LT_MASS_boot_pred_res_sum)

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/lmer/vignettes/model_evaluation.pdf

## Easy to generate effects plots with lmer 
ae_MASS <- allEffects(As_LMER_LT_MASS)

par(mfrow = c(1,2))
plot(ae_MASS$WEIGHT_GRAM_LOG$data$VALUE_LOG ~ ae_MASS$WEIGHT_GRAM_LOG$data$WEIGHT_GRAM_LOG, main = "MASS")  
lines(ae_MASS$WEIGHT_GRAM_LOG$fit ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$upper ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, lty = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$lower ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, lty = 2)

(pp <- sjPlot::plot_model(As_LMER_LT_MASS,type=c("pred"),
                          terms=c("WEIGHT_GRAM_LOG","WATERBODY_CODE"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))

## 1.x.x) Northern Pike ----

As_LMER_NP_MASS <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_NP) 
summary(As_LMER_NP_MASS)

hist(resid(As_LMER_NP_MASS), breaks = 100) # seems good 
plot(resid(As_LMER_NP_MASS) ~ fitted(As_LMER_NP_MASS)) # seems reasonable; a few very high residuals 

## Generate predictions for each waterbody based on the median size of NP in whole dataset 

As_LMER_NP_MASS_predict <- unique(As_NP[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
As_LMER_NP_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

## Bootstrap confidence intervals from bootMer
clust <- LMER_boot_initiate(varlist = "As_LMER_NP_MASS_predict", nclust = 4)
system.time(As_LMER_NP_MASS_boot_est <- lme4::bootMer(As_LMER_NP_MASS, LMER_boot_est, 
                                                      nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                      parallel = "snow", 
                                                      cl = clust, 
                                                      ncpus = 4)) #17s for 100
LMER_boot_est_summary(As_LMER_NP_MASS_boot_est) ## 889s

As_LMER_NP_MASS_boot_pred <- function(., newdata) {
  predict(., newdata=As_LMER_NP_MASS_predict)
}

clust <- LMER_boot_initiate(varlist = "As_LMER_NP_MASS_predict", nclust = 4)
system.time(As_LMER_NP_MASS_boot_pred_res <- lme4::bootMer(As_LMER_NP_MASS, As_LMER_NP_MASS_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 

(As_LMER_NP_MASS_boot_pred_res_sum <- LMER_boot_pred_summary(As_LMER_NP_MASS_boot_pred_res))

## ResuNPs 
LMER_boot_est_summary(As_LMER_NP_MASS_boot_est)
head(As_LMER_NP_MASS_boot_pred_res_sum)

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/lmer/vignettes/model_evaluation.pdf

## Easy to generate effects plots with lmer 
ae_MASS <- allEffects(As_LMER_NP_MASS)

par(mfrow = c(1,2))
plot(ae_MASS$WEIGHT_GRAM_LOG$data$VALUE_LOG ~ ae_MASS$WEIGHT_GRAM_LOG$data$WEIGHT_GRAM_LOG, main = "MASS")  
lines(ae_MASS$WEIGHT_GRAM_LOG$fit ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$upper ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, NPy = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$lower ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, NPy = 2)

(pp <- sjPlot::plot_model(As_LMER_NP_MASS,type=c("pred"),
                          terms=c("WEIGHT_GRAM_LOG","WATERBODY_CODE"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))

## 1.x.x) Walleye ----

As_LMER_WE_MASS <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + (WEIGHT_GRAM_LOG|WATERBODY_CODE) + (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_WE) 
summary(As_LMER_WE_MASS)

hist(resid(As_LMER_WE_MASS), breaks = 100) # seems good 
plot(resid(As_LMER_WE_MASS) ~ fitted(As_LMER_WE_MASS)) # seems reasonable; a few very high residuals 

## Generate predictions for each waterbody based on the median size of WE in whole dataset 

As_LMER_WE_MASS_predict <- unique(As_WE[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
As_LMER_WE_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

## Bootstrap confidence intervals from bootMer
clust <- LMER_boot_initiate(varlist = "As_LMER_WE_MASS_predict", nclust = 4)
system.time(As_LMER_WE_MASS_boot_est <- lme4::bootMer(As_LMER_WE_MASS, LMER_boot_est, 
                                                      nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                      parallel = "snow", 
                                                      cl = clust, 
                                                      ncpus = 4)) #17s for 100
LMER_boot_est_summary(As_LMER_WE_MASS_boot_est) ## 889s

As_LMER_WE_MASS_boot_pred <- function(., newdata) {
  predict(., newdata=As_LMER_WE_MASS_predict)
}

clust <- LMER_boot_initiate(varlist = "As_LMER_WE_MASS_predict", nclust = 4)
system.time(As_LMER_WE_MASS_boot_pred_res <- lme4::bootMer(As_LMER_WE_MASS, As_LMER_WE_MASS_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 

(As_LMER_WE_MASS_boot_pred_res_sum <- LMER_boot_pred_summary(As_LMER_WE_MASS_boot_pred_res))

## ResuWEs 
LMER_boot_est_summary(As_LMER_WE_MASS_boot_est)
head(As_LMER_WE_MASS_boot_pred_res_sum)

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/lmer/vignettes/model_evaluation.pdf

## Easy to generate effects plots with lmer 
ae_MASS <- allEffects(As_LMER_WE_MASS)

par(mfrow = c(1,2))
plot(ae_MASS$WEIGHT_GRAM_LOG$data$VALUE_LOG ~ ae_MASS$WEIGHT_GRAM_LOG$data$WEIGHT_GRAM_LOG, main = "MASS")  
lines(ae_MASS$WEIGHT_GRAM_LOG$fit ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$upper ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, WEy = 2)
lines(ae_MASS$WEIGHT_GRAM_LOG$lower ~ ae_MASS$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2, WEy = 2)

(pp <- sjPlot::plot_model(As_LMER_WE_MASS,type=c("pred"),
                          terms=c("WEIGHT_GRAM_LOG","WATERBODY_CODE"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))


## 2.5) INLA - mixed model, Bayesian, fast ----

## 1.x.x) Lake Trout ---- 

As_LT$WATERBODY_CODE <- factor(As_LT$WATERBODY_CODE) # for inla model
As_LT$WATERBODY_CODE1 <- as.integer(As_LT$WATERBODY_CODE) # for inla model
As_LT$WATERBODY_CODE2 <- As_LT$WATERBODY_CODE1 + max(As_LT$WATERBODY_CODE1) # for inla model
As_LT_n_waterbody <- dplyr::n_distinct(As_LT$WATERBODY_CODE) # for inla model

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
As_INLA_LT_MASS <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                       f(WATERBODY_CODE1, n = 2 * As_LT_n_waterbody, model = "iid2d") +   
                       f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                       f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                     data = As_LT, 
                     control.predictor = list(
                       compute = TRUE, 
                       quantiles = c(0.025, 0.5, 0.975)
                     ),
                     control.compute = list(
                       cpo = TRUE
                     )
)
summary(As_INLA_LT_MASS)

## There is a lot of information in the INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
As_INLA_LT_MASS_fits <- data.frame(WATERBODY_CODE = As_LT$WATERBODY_CODE, 
                                VALUE_LOG = As_LT$VALUE_LOG)

## Alternative comparison to fitted and residuals
As_INLA_LT_MASS_fits$inla_posterior_q50 <- As_INLA_LT_MASS$summary.fitted.values[, "0.5quant"]
As_INLA_LT_MASS_fits$inla_posterior_q2p5 <- As_INLA_LT_MASS$summary.fitted.values[, "0.025quant"]
As_INLA_LT_MASS_fits$inla_posterior_q97p5 <- As_INLA_LT_MASS$summary.fitted.values[, "0.975quant"]
As_INLA_LT_MASS_fits$resid_inla <- As_INLA_LT_MASS_fits$VALUE_LOG - As_INLA_LT_MASS_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = As_INLA_LT_MASS_fits)
j <- order(As_INLA_LT_MASS_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = As_INLA_LT_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ As_INLA_LT_MASS_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(As_INLA_LT_MASS_fits$VALUE_LOG ~ As_INLA_LT_MASS_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(As_INLA_LT_MASS_fits$VALUE_LOG ~ As_INLA_LT_MASS_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(As_INLA_LT_MASS_fits$resid_inla, breaks = 100) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
As_INLA_LT_MASS_fixed <- data.frame(
  ID = rownames(As_INLA_LT_MASS$summary.fixed),
  As_INLA_LT_MASS$summary.fixed, stringsAsFactors = FALSE
)
names(As_INLA_LT_MASS_fixed) <- c("ID", names(As_INLA_LT_MASS$summary.fixed))
As_INLA_LT_MASS_fixed$Type <- "Fixed"
head(As_INLA_LT_MASS_fixed)

As_INLA_LT_MASS_random_intercept <- As_INLA_LT_MASS$summary.random$WATERBODY_CODE1
As_INLA_LT_MASS_random_intercept$ID <- as.character(As_INLA_LT_MASS_random_intercept$ID)
As_INLA_LT_MASS_random_intercept$WATERBODY_CODE1 <- As_INLA_LT_MASS_random_intercept$ID
As_INLA_LT_MASS_random_intercept <- merge(As_INLA_LT_MASS_random_intercept, LT[,c("WATERBODY_CODE1", "WATERBODY_CODE")], no.dups = TRUE)
As_INLA_LT_MASS_random_intercept <- As_INLA_LT_MASS_random_intercept[!duplicated(As_INLA_LT_MASS_random_intercept),]
As_INLA_LT_MASS_random_intercept$Type <- "Random Intercept - Waterbody"
head(As_INLA_LT_MASS_random_intercept)

As_INLA_LT_MASS_random_slope <- As_INLA_LT_MASS$summary.random$WATERBODY_CODE2
As_INLA_LT_MASS_random_slope$ID <- as.character(As_INLA_LT_MASS_random_slope$ID)
As_INLA_LT_MASS_random_slope$WATERBODY_CODE2 <- As_INLA_LT_MASS_random_slope$ID
As_INLA_LT_MASS_random_slope <- merge(As_INLA_LT_MASS_random_slope, LT[,c("WATERBODY_CODE2", "WATERBODY_CODE")])
As_INLA_LT_MASS_random_slope <- As_INLA_LT_MASS_random_slope[!duplicated(As_INLA_LT_MASS_random_slope),]
As_INLA_LT_MASS_random_slope$Type <- "Random Slope - Waterbody"
head(As_INLA_LT_MASS_random_slope)

As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_LT_MASS$summary.random$WATERBODY_CODE_SAMPLE_YEAR
As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$ID
As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR[!duplicated(As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR),]
As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$Type <- "Random Intercept - Waterbody Sampling Year "
head(As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR)

As_INLA_LT_MASS_summary <- list(
  As_INLA_LT_MASS_fixed,
  As_INLA_LT_MASS_random_intercept,
  As_INLA_LT_MASS_random_slope,
  As_INLA_LT_MASS_random_WATERBODY_CODE_SAMPLE_YEAR
)

As_INLA_LT_MASS_summary <- dplyr::bind_rows(As_INLA_LT_MASS_summary)
head(As_INLA_LT_MASS_summary)
tail(As_INLA_LT_MASS_summary)

## Variance parameters, difficult to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
As_INLA_LT_MASS$marginals.hyperpar$`Precision for the Gaussian observations`
As_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
As_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
As_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`

tau_WATERBODY_CODE1 <- As_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- As_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_LT_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`
tau_residual <- As_INLA_LT_MASS$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_WATERBODY_CODE_SAMPLE_YEAR <- inla.emarginal(MySqrt, tau_WATERBODY_CODE_SAMPLE_YEAR))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(As_LT_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)^2)

## TMB vs. INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)
summary(As_TMB_LT_MASS)

## WATERBODY_CODE absolute variation >> WATERBODY_CODE VALUE_LOG~LENGTH_CM relation > WATERBODY_SAMPLE_YEAR > RESIDUAL 

## Now predictions 

## Predictions in INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the VALUE_LOG values are set to NA
## Then you extract the posterior information 

## Take LMER predictions but fit the Waterbody Codes to the dataframe 
head(As_LMER_LT_MASS_predict) 

As_INLA_LT_MASS_predict <- merge(As_LMER_LT_MASS_predict, As_LT[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
As_INLA_LT_MASS_predict <- As_INLA_LT_MASS_predict[!duplicated(As_INLA_LT_MASS_predict),]
As_INLA_LT_MASS_predict$VALUE_LOG <- NA

nrow(As_TMB_LT_MASS_predict)
nrow(As_INLA_LT_MASS_predict)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(As_TMB_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, As_INLA_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

As_INLA_LT_MASS_prediction <- As_LT[names(As_INLA_LT_MASS_predict)]

As_INLA_LT_MASS_prediction <- rbind(As_INLA_LT_MASS_prediction, 
                                 As_INLA_LT_MASS_predict)

head(As_INLA_LT_MASS_prediction)
tail(As_INLA_LT_MASS_prediction)

nrow(As_INLA_LT_MASS_prediction) - nrow(As_LT)

As_INLA_LT_MASS_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                        f(WATERBODY_CODE1, n = 2 * As_LT_n_waterbody, model = "iid2d") +   
                                        f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                        f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                      data = As_INLA_LT_MASS_prediction, 
                                      control.predictor = list(
                                        compute = TRUE, 
                                        quantiles = c(0.025, 0.5, 0.975)
                                      ),
                                      control.compute = list(
                                        cpo = TRUE
                                      )
)

As_INLA_LT_MASS_prediction_posteriors <- As_INLA_LT_MASS_prediction_model$summary.fitted.values[
  (nrow(As_LT) + 1):nrow(As_INLA_LT_MASS_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(As_INLA_LT_MASS_prediction_posteriors)

As_INLA_LT_MASS_prediction_posteriors <- exp(As_INLA_LT_MASS_prediction_posteriors)

As_INLA_LT_MASS_predict <- cbind(As_INLA_LT_MASS_predict, As_INLA_LT_MASS_prediction_posteriors)
As_INLA_LT_MASS_predict <- As_INLA_LT_MASS_predict[!duplicated(As_INLA_LT_MASS_predict),]

## Does this condition still hold true for comparison?
identical(As_TMB_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, As_INLA_LT_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

## Produces very similar results but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(exp(As_LMER_LT_MASS_boot_pred_res_sum$fit)  ~ As_INLA_LT_MASS_predict$`0.5quant`)
abline(0,1)
head(As_TMB_LT_MASS_predict)
head(As_INLA_LT_MASS_predict)

## 1.x.x) Northern Pike ---- 

As_NP$WATERBODY_CODE <- factor(As_NP$WATERBODY_CODE) # for inla model
As_NP$WATERBODY_CODE1 <- as.integer(As_NP$WATERBODY_CODE) # for inla model
As_NP$WATERBODY_CODE2 <- As_NP$WATERBODY_CODE1 + max(As_NP$WATERBODY_CODE1) # for inla model
As_NP_n_waterbody <- dplyr::n_distinct(As_NP$WATERBODY_CODE) # for inla model

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
As_INLA_NP_MASS <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                          f(WATERBODY_CODE1, n = 2 * As_NP_n_waterbody, model = "iid2d") +   
                          f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                          f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                        data = As_NP, 
                        control.predictor = list(
                          compute = TRUE, 
                          quantiles = c(0.025, 0.5, 0.975)
                        ),
                        control.compute = list(
                          cpo = TRUE
                        )
)
summary(As_INLA_NP_MASS)

## There is a lot of information in the INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
As_INLA_NP_MASS_fits <- data.frame(WATERBODY_CODE = As_NP$WATERBODY_CODE, 
                                   VALUE_LOG = As_NP$VALUE_LOG)

## ANPernative comparison to fitted and residuals
As_INLA_NP_MASS_fits$inla_posterior_q50 <- As_INLA_NP_MASS$summary.fitted.values[, "0.5quant"]
As_INLA_NP_MASS_fits$inla_posterior_q2p5 <- As_INLA_NP_MASS$summary.fitted.values[, "0.025quant"]
As_INLA_NP_MASS_fits$inla_posterior_q97p5 <- As_INLA_NP_MASS$summary.fitted.values[, "0.975quant"]
As_INLA_NP_MASS_fits$resid_inla <- As_INLA_NP_MASS_fits$VALUE_LOG - As_INLA_NP_MASS_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = As_INLA_NP_MASS_fits)
j <- order(As_INLA_NP_MASS_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = As_INLA_NP_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ As_INLA_NP_MASS_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(As_INLA_NP_MASS_fits$VALUE_LOG ~ As_INLA_NP_MASS_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(As_INLA_NP_MASS_fits$VALUE_LOG ~ As_INLA_NP_MASS_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(As_INLA_NP_MASS_fits$resid_inla, breaks = 100) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
As_INLA_NP_MASS_fixed <- data.frame(
  ID = rownames(As_INLA_NP_MASS$summary.fixed),
  As_INLA_NP_MASS$summary.fixed, stringsAsFactors = FALSE
)
names(As_INLA_NP_MASS_fixed) <- c("ID", names(As_INLA_NP_MASS$summary.fixed))
As_INLA_NP_MASS_fixed$Type <- "Fixed"
head(As_INLA_NP_MASS_fixed)

As_INLA_NP_MASS_random_intercept <- As_INLA_NP_MASS$summary.random$WATERBODY_CODE1
As_INLA_NP_MASS_random_intercept$ID <- as.character(As_INLA_NP_MASS_random_intercept$ID)
As_INLA_NP_MASS_random_intercept$WATERBODY_CODE1 <- As_INLA_NP_MASS_random_intercept$ID
As_INLA_NP_MASS_random_intercept <- merge(As_INLA_NP_MASS_random_intercept, NP[,c("WATERBODY_CODE1", "WATERBODY_CODE")], no.dups = TRUE)
As_INLA_NP_MASS_random_intercept <- As_INLA_NP_MASS_random_intercept[!duplicated(As_INLA_NP_MASS_random_intercept),]
As_INLA_NP_MASS_random_intercept$Type <- "Random Intercept - Waterbody"
head(As_INLA_NP_MASS_random_intercept)

As_INLA_NP_MASS_random_slope <- As_INLA_NP_MASS$summary.random$WATERBODY_CODE2
As_INLA_NP_MASS_random_slope$ID <- as.character(As_INLA_NP_MASS_random_slope$ID)
As_INLA_NP_MASS_random_slope$WATERBODY_CODE2 <- As_INLA_NP_MASS_random_slope$ID
As_INLA_NP_MASS_random_slope <- merge(As_INLA_NP_MASS_random_slope, NP[,c("WATERBODY_CODE2", "WATERBODY_CODE")])
As_INLA_NP_MASS_random_slope <- As_INLA_NP_MASS_random_slope[!duplicated(As_INLA_NP_MASS_random_slope),]
As_INLA_NP_MASS_random_slope$Type <- "Random Slope - Waterbody"
head(As_INLA_NP_MASS_random_slope)

As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_NP_MASS$summary.random$WATERBODY_CODE_SAMPLE_YEAR
As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$ID
As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR[!duplicated(As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR),]
As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$Type <- "Random Intercept - Waterbody Sampling Year "
head(As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR)

As_INLA_NP_MASS_summary <- list(
  As_INLA_NP_MASS_fixed,
  As_INLA_NP_MASS_random_intercept,
  As_INLA_NP_MASS_random_slope,
  As_INLA_NP_MASS_random_WATERBODY_CODE_SAMPLE_YEAR
)

As_INLA_NP_MASS_summary <- dplyr::bind_rows(As_INLA_NP_MASS_summary)
head(As_INLA_NP_MASS_summary)
tail(As_INLA_NP_MASS_summary)

## Variance parameters, difficuNP to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
As_INLA_NP_MASS$marginals.hyperpar$`Precision for the Gaussian observations`
As_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
As_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
As_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`

tau_WATERBODY_CODE1 <- As_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- As_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_NP_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`
tau_residual <- As_INLA_NP_MASS$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_WATERBODY_CODE_SAMPLE_YEAR <- inla.emarginal(MySqrt, tau_WATERBODY_CODE_SAMPLE_YEAR))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(As_NP_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)^2)

## TMB vs. INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)
summary(As_TMB_NP_MASS)

## WATERBODY_CODE absolute variation >> WATERBODY_CODE VALUE_LOG~LENGTH_CM relation > WATERBODY_SAMPLE_YEAR > RESIDUAL 

## Now predictions 

## Predictions in INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the VALUE_LOG values are set to NA
## Then you extract the posterior information 

## Take LMER predictions but fit the Waterbody Codes to the dataframe 
head(As_LMER_NP_MASS_predict) 

As_INLA_NP_MASS_predict <- merge(As_LMER_NP_MASS_predict, As_NP[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
As_INLA_NP_MASS_predict <- As_INLA_NP_MASS_predict[!duplicated(As_INLA_NP_MASS_predict),]
As_INLA_NP_MASS_predict$VALUE_LOG <- NA

nrow(As_TMB_NP_MASS_predict)
nrow(As_INLA_NP_MASS_predict)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(As_TMB_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, As_INLA_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

As_INLA_NP_MASS_prediction <- As_NP[names(As_INLA_NP_MASS_predict)]

As_INLA_NP_MASS_prediction <- rbind(As_INLA_NP_MASS_prediction, 
                                    As_INLA_NP_MASS_predict)

head(As_INLA_NP_MASS_prediction)
tail(As_INLA_NP_MASS_prediction)

nrow(As_INLA_NP_MASS_prediction) - nrow(As_NP)

As_INLA_NP_MASS_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * As_NP_n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = As_INLA_NP_MASS_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
)

As_INLA_NP_MASS_prediction_posteriors <- As_INLA_NP_MASS_prediction_model$summary.fitted.values[
  (nrow(As_NP) + 1):nrow(As_INLA_NP_MASS_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(As_INLA_NP_MASS_prediction_posteriors)

As_INLA_NP_MASS_prediction_posteriors <- exp(As_INLA_NP_MASS_prediction_posteriors)

As_INLA_NP_MASS_predict <- cbind(As_INLA_NP_MASS_predict, As_INLA_NP_MASS_prediction_posteriors)
As_INLA_NP_MASS_predict <- As_INLA_NP_MASS_predict[!duplicated(As_INLA_NP_MASS_predict),]

## Does this condition still hold true for comparison?
identical(As_TMB_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, As_INLA_NP_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

## Produces very similar resuNPs but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(exp(As_LMER_NP_MASS_boot_pred_res_sum$fit)  ~ As_INLA_NP_MASS_predict$`0.5quant`)
abline(0,1)
head(As_TMB_NP_MASS_predict)
head(As_INLA_NP_MASS_predict)

## 1.x.x) Walleye ---- 

As_WE$WATERBODY_CODE <- factor(As_WE$WATERBODY_CODE) # for inla model
As_WE$WATERBODY_CODE1 <- as.integer(As_WE$WATERBODY_CODE) # for inla model
As_WE$WATERBODY_CODE2 <- As_WE$WATERBODY_CODE1 + max(As_WE$WATERBODY_CODE1) # for inla model
As_WE_n_waterbody <- dplyr::n_distinct(As_WE$WATERBODY_CODE) # for inla model

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
As_INLA_WE_MASS <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                          f(WATERBODY_CODE1, n = 2 * As_WE_n_waterbody, model = "iid2d") +   
                          f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                          f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                        data = As_WE, 
                        control.predictor = list(
                          compute = TRUE, 
                          quantiles = c(0.025, 0.5, 0.975)
                        ),
                        control.compute = list(
                          cpo = TRUE
                        )
)
summary(As_INLA_WE_MASS)

## There is a lot of information in the INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
As_INLA_WE_MASS_fits <- data.frame(WATERBODY_CODE = As_WE$WATERBODY_CODE, 
                                   VALUE_LOG = As_WE$VALUE_LOG)

## AWEernative comparison to fitted and residuals
As_INLA_WE_MASS_fits$inla_posterior_q50 <- As_INLA_WE_MASS$summary.fitted.values[, "0.5quant"]
As_INLA_WE_MASS_fits$inla_posterior_q2p5 <- As_INLA_WE_MASS$summary.fitted.values[, "0.025quant"]
As_INLA_WE_MASS_fits$inla_posterior_q97p5 <- As_INLA_WE_MASS$summary.fitted.values[, "0.975quant"]
As_INLA_WE_MASS_fits$resid_inla <- As_INLA_WE_MASS_fits$VALUE_LOG - As_INLA_WE_MASS_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = As_INLA_WE_MASS_fits)
j <- order(As_INLA_WE_MASS_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = As_INLA_WE_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ As_INLA_WE_MASS_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(As_INLA_WE_MASS_fits$VALUE_LOG ~ As_INLA_WE_MASS_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(As_INLA_WE_MASS_fits$VALUE_LOG ~ As_INLA_WE_MASS_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(As_INLA_WE_MASS_fits$resid_inla, breaks = 100) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
As_INLA_WE_MASS_fixed <- data.frame(
  ID = rownames(As_INLA_WE_MASS$summary.fixed),
  As_INLA_WE_MASS$summary.fixed, stringsAsFactors = FALSE
)
names(As_INLA_WE_MASS_fixed) <- c("ID", names(As_INLA_WE_MASS$summary.fixed))
As_INLA_WE_MASS_fixed$Type <- "Fixed"
head(As_INLA_WE_MASS_fixed)

As_INLA_WE_MASS_random_intercept <- As_INLA_WE_MASS$summary.random$WATERBODY_CODE1
As_INLA_WE_MASS_random_intercept$ID <- as.character(As_INLA_WE_MASS_random_intercept$ID)
As_INLA_WE_MASS_random_intercept$WATERBODY_CODE1 <- As_INLA_WE_MASS_random_intercept$ID
As_INLA_WE_MASS_random_intercept <- merge(As_INLA_WE_MASS_random_intercept, WE[,c("WATERBODY_CODE1", "WATERBODY_CODE")], no.dups = TRUE)
As_INLA_WE_MASS_random_intercept <- As_INLA_WE_MASS_random_intercept[!duplicated(As_INLA_WE_MASS_random_intercept),]
As_INLA_WE_MASS_random_intercept$Type <- "Random Intercept - Waterbody"
head(As_INLA_WE_MASS_random_intercept)

As_INLA_WE_MASS_random_slope <- As_INLA_WE_MASS$summary.random$WATERBODY_CODE2
As_INLA_WE_MASS_random_slope$ID <- as.character(As_INLA_WE_MASS_random_slope$ID)
As_INLA_WE_MASS_random_slope$WATERBODY_CODE2 <- As_INLA_WE_MASS_random_slope$ID
As_INLA_WE_MASS_random_slope <- merge(As_INLA_WE_MASS_random_slope, WE[,c("WATERBODY_CODE2", "WATERBODY_CODE")])
As_INLA_WE_MASS_random_slope <- As_INLA_WE_MASS_random_slope[!duplicated(As_INLA_WE_MASS_random_slope),]
As_INLA_WE_MASS_random_slope$Type <- "Random Slope - Waterbody"
head(As_INLA_WE_MASS_random_slope)

As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_WE_MASS$summary.random$WATERBODY_CODE_SAMPLE_YEAR
As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$ID
As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR[!duplicated(As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR),]
As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR$Type <- "Random Intercept - Waterbody Sampling Year "
head(As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR)

As_INLA_WE_MASS_summary <- list(
  As_INLA_WE_MASS_fixed,
  As_INLA_WE_MASS_random_intercept,
  As_INLA_WE_MASS_random_slope,
  As_INLA_WE_MASS_random_WATERBODY_CODE_SAMPLE_YEAR
)

As_INLA_WE_MASS_summary <- dplyr::bind_rows(As_INLA_WE_MASS_summary)
head(As_INLA_WE_MASS_summary)
tail(As_INLA_WE_MASS_summary)

## Variance parameters, difficuWE to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
As_INLA_WE_MASS$marginals.hyperpar$`Precision for the Gaussian observations`
As_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
As_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
As_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`

tau_WATERBODY_CODE1 <- As_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- As_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_WATERBODY_CODE_SAMPLE_YEAR <- As_INLA_WE_MASS$marginals.hyperpar$`Precision for WATERBODY_CODE_SAMPLE_YEAR`
tau_residual <- As_INLA_WE_MASS$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_WATERBODY_CODE_SAMPLE_YEAR <- inla.emarginal(MySqrt, tau_WATERBODY_CODE_SAMPLE_YEAR))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(As_WE_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)^2)

## TMB vs. INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_WATERBODY_CODE_SAMPLE_YEAR)
summary(As_TMB_WE_MASS)

## WATERBODY_CODE absolute variation >> WATERBODY_CODE VALUE_LOG~LENGTH_CM relation > WATERBODY_SAMPLE_YEAR > RESIDUAL 

## Now predictions 

## Predictions in INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the VALUE_LOG values are set to NA
## Then you extract the posterior information 

## Take LMER predictions but fit the Waterbody Codes to the dataframe 
head(As_LMER_WE_MASS_predict) 

As_INLA_WE_MASS_predict <- merge(As_LMER_WE_MASS_predict, As_WE[,c("WATERBODY_CODE_SAMPLE_YEAR", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
As_INLA_WE_MASS_predict <- As_INLA_WE_MASS_predict[!duplicated(As_INLA_WE_MASS_predict),]
As_INLA_WE_MASS_predict$VALUE_LOG <- NA

nrow(As_TMB_WE_MASS_predict)
nrow(As_INLA_WE_MASS_predict)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(As_TMB_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, As_INLA_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

As_INLA_WE_MASS_prediction <- As_WE[names(As_INLA_WE_MASS_predict)]

As_INLA_WE_MASS_prediction <- rbind(As_INLA_WE_MASS_prediction, 
                                    As_INLA_WE_MASS_predict)

head(As_INLA_WE_MASS_prediction)
tail(As_INLA_WE_MASS_prediction)

nrow(As_INLA_WE_MASS_prediction) - nrow(As_WE)

As_INLA_WE_MASS_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * As_WE_n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = As_INLA_WE_MASS_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
)

As_INLA_WE_MASS_prediction_posteriors <- As_INLA_WE_MASS_prediction_model$summary.fitted.values[
  (nrow(As_WE) + 1):nrow(As_INLA_WE_MASS_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(As_INLA_WE_MASS_prediction_posteriors)

As_INLA_WE_MASS_prediction_posteriors <- exp(As_INLA_WE_MASS_prediction_posteriors)

As_INLA_WE_MASS_predict <- cbind(As_INLA_WE_MASS_predict, As_INLA_WE_MASS_prediction_posteriors)
As_INLA_WE_MASS_predict <- As_INLA_WE_MASS_predict[!duplicated(As_INLA_WE_MASS_predict),]

## Does this condition still hold true for comparison?
identical(As_TMB_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR, As_INLA_WE_MASS_predict$WATERBODY_CODE_SAMPLE_YEAR)

## Produces very similar resuWEs but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(exp(As_LMER_WE_MASS_boot_pred_res_sum$fit)  ~ As_INLA_WE_MASS_predict$`0.5quant`)
abline(0,1)

## 2.6) RSTAN - mixed model, Bayesian slow ---- 

## 2.x.x) Lake Trout ---- 

system.time(As_STAN_LT_MASS <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_LT, 
                                          cores = 4, chains = 4, iter = 8000, adapt_delta = 0.99))
saveRDS(As_STAN_LT_MASS, "./out_workspaces/As_STAN_LT_MASS.rds") ## 156s
As_STAN_LT_MASS <- readRDS("./out_workspaces/As_STAN_LT_MASS.rds")

# shinystan
#launch_shinystan(As_STAN_LT_MASS)

# Marginal r-squared (no random effects)
As_STAN_LT_MASS_rsq_marg <- bayes_R2(As_STAN_LT_MASS, re.form = NA)
median(As_STAN_LT_MASS_rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
As_STAN_LT_MASS_rsq_cond_lake <- bayes_R2(As_STAN_LT_MASS, re.form = ~ (WEIGHT_GRAM_LOG|WATERBODY_CODE) )
median(As_STAN_LT_MASS_rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
As_STAN_LT_MASS_fits <- data.frame(
  VALUE_LOG = As_LT$VALUE_LOG,
  As_STAN_FIT = fitted(As_STAN_LT_MASS),
  As_STAN_RESID = residuals(As_STAN_LT_MASS)
)

lw1_stan <- loess(As_STAN_RESID ~ As_STAN_FIT, data = As_STAN_LT_MASS_fits)
j <- order(As_STAN_LT_MASS_fits$As_STAN_FIT)
plot(As_STAN_RESID ~ As_STAN_FIT, data = As_STAN_LT_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ As_STAN_LT_MASS_fits$As_STAN_FIT[j], lwd = 2)

## Histogram of residuals
hist(As_STAN_LT_MASS_fits$As_STAN_RESID, breaks = 100)

## Alternative comparison to fitted and residuals
As_STAN_LT_MASS_posterior_full <- rstanarm::posterior_predict(As_STAN_LT_MASS)
As_STAN_LT_MASS_posterior_est_full <- apply(As_STAN_LT_MASS_posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

As_STAN_LT_MASS_fits$stan_posterior_q50 <- As_STAN_LT_MASS_posterior_est_full[2, ]
As_STAN_LT_MASS_fits$stan_posterior_q2p5 <- As_STAN_LT_MASS_posterior_est_full[1, ]
As_STAN_LT_MASS_fits$stan_posterior_q97p5 <- As_STAN_LT_MASS_posterior_est_full[3, ]

## Parameter estimates from posterior distribution
As_STAN_LT_MASS_sims <- as.data.frame(as.matrix(As_STAN_LT_MASS))
names(As_STAN_LT_MASS_sims)

a_quant <- apply(
  X = As_STAN_LT_MASS_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

As_STAN_LT_MASS_summary <- a_quant
As_STAN_LT_MASS_summary$param

As_STAN_LT_MASS_summary

## Variance parameters
As_STAN_LT_RE <- As_STAN_LT_MASS_sims[,grep("Sigma", names(As_STAN_LT_MASS_sims), ignore.case = TRUE)]

As_STAN_LT_RE <- apply(
  X = As_STAN_LT_RE,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.50, 0.025, 0.975)
)

sqrt(As_STAN_LT_RE)

## Predictions 

As_STAN_LT_MASS_predict <- unique(As_LT[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
As_STAN_LT_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

As_STAN_LT_MASS_predictions <- rstanarm::posterior_predict(As_STAN_LT_MASS,
                                                        newdata = As_STAN_LT_MASS_predict
)

As_STAN_LT_MASS_predictions <- apply(As_STAN_LT_MASS_predictions,
                                  MARGIN = 2,
                                  function(x) {
                                    quantile(x, probs = c(0.025, 0.5, 0.975))
                                  }
)

As_STAN_LT_MASS_predictions <- exp(t(As_STAN_LT_MASS_predictions))
colnames(As_STAN_LT_MASS_predictions) <- c("As_STAN_0.025quant", "As_STAN_0.5quant", "As_STAN_0.975quant")
As_STAN_LT_MASS_predictions <- cbind(As_STAN_LT_MASS_predict, As_STAN_LT_MASS_predictions)

par(mfrow = c(1,2))
plot(exp(As_LMER_LT_MASS_boot_pred_res_sum$fit), As_STAN_LT_MASS_predictions$As_STAN_0.5quant, xlab = "Maximum Likelihood Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))
plot((As_INLA_LT_MASS_prediction_posteriors$`0.5quant`), As_STAN_LT_MASS_predictions$As_STAN_0.5quant,  xlab = "Bayesian INLA Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))

## 2.x.x) Northern Pike ---- 

system.time(As_STAN_NP_MASS <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_NP,
                                          cores = 4, chains = 4, iter = 8000, adapt_delta = 0.99))
saveRDS(As_STAN_NP_MASS, "./out_workspaces/As_STAN_NP_MASS.rds") ## 936s
As_STAN_NP_MASS <- readRDS("./out_workspaces/As_STAN_NP_MASS.rds")

# shinystan
#launch_shinystan(STAN_NP_MASS)

# Marginal r-squared (no random effects)
As_STAN_NP_MASS_rsq_marg <- bayes_R2(As_STAN_NP_MASS, re.form = NA)
median(As_STAN_NP_MASS_rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
As_STAN_NP_MASS_rsq_cond_lake <- bayes_R2(As_STAN_NP_MASS, re.form = ~ (WEIGHT_GRAM_LOG|WATERBODY_CODE) )
median(As_STAN_NP_MASS_rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
As_STAN_NP_MASS_fits <- data.frame(
  VALUE_LOG = As_NP$VALUE_LOG,
  As_STAN_FIT = fitted(As_STAN_NP_MASS),
  As_STAN_RESID = residuals(As_STAN_NP_MASS)
)

lw1_stan <- loess(As_STAN_RESID ~ As_STAN_FIT, data = As_STAN_NP_MASS_fits)
j <- order(As_STAN_NP_MASS_fits$As_STAN_FIT)
plot(As_STAN_RESID ~ As_STAN_FIT, data = As_STAN_NP_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ As_STAN_NP_MASS_fits$As_STAN_FIT[j], lwd = 2)

## Histogram of residuals
hist(As_STAN_NP_MASS_fits$As_STAN_RESID, breaks = 100)

## Alternative comparison to fitted and residuals
As_STAN_NP_MASS_posterior_full <- rstanarm::posterior_predict(As_STAN_NP_MASS)
As_STAN_NP_MASS_posterior_est_full <- apply(As_STAN_NP_MASS_posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

As_STAN_NP_MASS_fits$stan_posterior_q50 <- As_STAN_NP_MASS_posterior_est_full[2, ]
As_STAN_NP_MASS_fits$stan_posterior_q2p5 <- As_STAN_NP_MASS_posterior_est_full[1, ]
As_STAN_NP_MASS_fits$stan_posterior_q97p5 <- As_STAN_NP_MASS_posterior_est_full[3, ]

## Parameter estimates from posterior distribution
As_STAN_NP_MASS_sims <- as.data.frame(as.matrix(As_STAN_NP_MASS))
names(As_STAN_NP_MASS_sims)

a_quant <- apply(
  X = As_STAN_NP_MASS_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

As_STAN_NP_MASS_summary <- a_quant
As_STAN_NP_MASS_summary$param

## Variance parameters
As_STAN_NP_RE <- As_STAN_NP_MASS_sims[,grep("Sigma", names(As_STAN_NP_MASS_sims), ignore.case = TRUE)]

As_STAN_NP_RE <- apply(
  X = As_STAN_NP_RE,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.50, 0.025, 0.975)
)

sqrt(As_STAN_NP_RE)

## Predictions 
As_STAN_NP_MASS_predict <- unique(As_NP[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
As_STAN_NP_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

As_STAN_NP_MASS_predictions <- rstanarm::posterior_predict(As_STAN_NP_MASS,
                                                        newdata = As_STAN_NP_MASS_predict
)

As_STAN_NP_MASS_predictions <- apply(As_STAN_NP_MASS_predictions,
                                  MARGIN = 2,
                                  function(x) {
                                    quantile(x, probs = c(0.025, 0.5, 0.975))
                                  }
)

As_STAN_NP_MASS_predictions <- exp(t(As_STAN_NP_MASS_predictions))
colnames(As_STAN_NP_MASS_predictions) <- c("As_STAN_0.025quant", "As_STAN_0.5quant", "As_STAN_0.975quant")
As_STAN_NP_MASS_predictions <- cbind(As_STAN_NP_MASS_predict, As_STAN_NP_MASS_predictions)

par(mfrow = c(1,2))
plot(exp(As_LMER_NP_MASS_boot_pred_res_sum$fit), As_STAN_NP_MASS_predictions$As_STAN_0.5quant, xlab = "Maximum Likelihood Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))
plot((As_INLA_NP_MASS_prediction_posteriors$`0.5quant`), As_STAN_NP_MASS_predictions$As_STAN_0.5quant,  xlab = "Bayesian INLA Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))

## 3.x.x) Walleye ----

system.time(As_STAN_WE_MASS <- stan_glmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_WE,
                                          cores = 4, chains = 4, iter = 8000, adapt_delta = 0.99))
saveRDS(As_STAN_WE_MASS, "./out_workspaces/As_STAN_WE_MASS.rds") ## 1081s
As_STAN_WE_MASS <- readRDS("./out_workspaces/As_STAN_WE_MASS.rds")

# shinystan
#launch_shinystan(As_STAN_WE_MASS)

# Marginal r-squared (no random effects)
As_STAN_WE_MASS_rsq_marg <- bayes_R2(As_STAN_WE_MASS, re.form = NA)
median(As_STAN_WE_MASS_rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
As_STAN_WE_MASS_rsq_cond_lake <- bayes_R2(As_STAN_WE_MASS, re.form = ~ (WEIGHT_GRAM_LOG|WATERBODY_CODE) )
median(As_STAN_WE_MASS_rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
As_STAN_WE_MASS_fits <- data.frame(
  VALUE_LOG = As_WE$VALUE_LOG,
  As_STAN_FIT = fitted(As_STAN_WE_MASS),
  As_STAN_RESID = residuals(As_STAN_WE_MASS)
)

lw1_stan <- loess(As_STAN_RESID ~ As_STAN_FIT, data = As_STAN_WE_MASS_fits)
j <- order(As_STAN_WE_MASS_fits$As_STAN_FIT)
plot(As_STAN_RESID ~ As_STAN_FIT, data = As_STAN_WE_MASS_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ As_STAN_WE_MASS_fits$As_STAN_FIT[j], lwd = 2)

## Histogram of residuals
hist(As_STAN_WE_MASS_fits$As_STAN_RESID, breaks = 100)

## Alternative comparison to fitted and residuals
As_STAN_WE_MASS_posterior_full <- rstanarm::posterior_predict(As_STAN_WE_MASS)
As_STAN_WE_MASS_posterior_est_full <- apply(As_STAN_WE_MASS_posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

As_STAN_WE_MASS_fits$stan_posterior_q50 <- As_STAN_WE_MASS_posterior_est_full[2, ]
As_STAN_WE_MASS_fits$stan_posterior_q2p5 <- As_STAN_WE_MASS_posterior_est_full[1, ]
As_STAN_WE_MASS_fits$stan_posterior_q97p5 <- As_STAN_WE_MASS_posterior_est_full[3, ]

## Parameter estimates from posterior distribution
As_STAN_WE_MASS_sims <- as.data.frame(as.matrix(As_STAN_WE_MASS))
names(As_STAN_WE_MASS_sims)

a_quant <- apply(
  X = As_STAN_WE_MASS_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

As_STAN_WE_MASS_summary <- a_quant
As_STAN_WE_MASS_summary$param

As_STAN_WE_MASS_summary

## Variance parameters
As_STAN_WE_RE <- As_STAN_WE_MASS_sims[,grep("Sigma", names(As_STAN_WE_MASS_sims), ignore.case = T)]

As_STAN_WE_RE <- apply(
  X = As_STAN_WE_RE,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.50, 0.025, 0.975)
)

sqrt(As_STAN_WE_RE)

## Predictions 
As_STAN_WE_MASS_predict <- unique(As_WE[,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
As_STAN_WE_MASS_predict$WEIGHT_GRAM_LOG <- log(1000)

As_STAN_WE_MASS_predictions <- rstanarm::posterior_predict(As_STAN_WE_MASS,
                                                        newdata = As_STAN_WE_MASS_predict
)

As_STAN_WE_MASS_predictions <- apply(As_STAN_WE_MASS_predictions,
                                  MARGIN = 2,
                                  function(x) {
                                    quantile(x, probs = c(0.025, 0.5, 0.975))
                                  }
)

As_STAN_WE_MASS_predictions <- exp(t(As_STAN_WE_MASS_predictions))
colnames(As_STAN_WE_MASS_predictions) <- c("As_STAN_0.025quant", "As_STAN_0.5quant", "As_STAN_0.975quant")
As_STAN_WE_MASS_predictions <- cbind(As_STAN_WE_MASS_predict, As_STAN_WE_MASS_predictions)

par(mfrow = c(1,2))
plot(exp(As_LMER_WE_MASS_boot_pred_res_sum$fit), As_STAN_WE_MASS_predictions$As_STAN_0.5quant, xlab = "Maximum Likelihood Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))
plot((As_INLA_WE_MASS_prediction_posteriors$`0.5quant`), As_STAN_WE_MASS_predictions$As_STAN_0.5quant,  xlab = "Bayesian INLA Predictions", ylab = "Bayesian MCMC Predictions", xlim = c(0,2), ylim = c(0,2))

## ***************
## MANUSCRIPT ----
## *************** 

## Weights for species of two datasets 

quantile(Hg_LT$WEIGHT_GRAM, probs = c(0.5, 0.025, 0.975))
quantile(Hg_NP$WEIGHT_GRAM, probs = c(0.5, 0.025, 0.975))
quantile(Hg_WE$WEIGHT_GRAM, probs = c(0.5, 0.025, 0.975))

quantile(As_LT$WEIGHT_GRAM, probs = c(0.5, 0.025, 0.975))
quantile(As_NP$WEIGHT_GRAM, probs = c(0.5, 0.025, 0.975))
quantile(As_WE$WEIGHT_GRAM, probs = c(0.5, 0.025, 0.975))

## Table XX Model Results ----

## Hg - SER ----

## LT
LT_hg_sub_preds_mass <- subset(hg_sub_preds_mass, SPECIES_CODE == "081")
median(LT_hg_sub_preds_mass$RMSE) 
min(LT_hg_sub_preds_mass$RMSE) 
max(LT_hg_sub_preds_mass$RMSE) 

median(LT_hg_sub_preds_mass$r2) 
min(LT_hg_sub_preds_mass$r2) 
max(LT_hg_sub_preds_mass$r2) 

sd(LT_hg_sub_preds_mass$int)
sd(LT_hg_sub_preds_mass$slp)

## NP
NP_hg_sub_preds_mass <- subset(hg_sub_preds_mass, SPECIES_CODE == "131")
median(NP_hg_sub_preds_mass$RMSE) 
min(NP_hg_sub_preds_mass$RMSE) 
max(NP_hg_sub_preds_mass$RMSE) 

median(NP_hg_sub_preds_mass$r2) 
min(NP_hg_sub_preds_mass$r2) 
max(NP_hg_sub_preds_mass$r2) 

sd(NP_hg_sub_preds_mass$int)
sd(NP_hg_sub_preds_mass$slp)

## WE
WE_hg_sub_preds_mass <- subset(hg_sub_preds_mass, SPECIES_CODE == "334")
median(WE_hg_sub_preds_mass$RMSE) 
min(WE_hg_sub_preds_mass$RMSE) 
max(WE_hg_sub_preds_mass$RMSE) 

median(WE_hg_sub_preds_mass$r2) 
min(WE_hg_sub_preds_mass$r2) 
max(WE_hg_sub_preds_mass$r2) 

sd(WE_hg_sub_preds_mass$int)
sd(WE_hg_sub_preds_mass$slp)

## Hg - ML ----

## LT
sqrt(mean(resid(Hg_LMER_LT_MASS)^2))
summary(lm(Hg_LT$VALUE_LOG ~ fitted(Hg_LMER_LT_MASS)))$adj.r.squared
LMER_boot_est_summary(Hg_LMER_LT_MASS_boot_est)


## NP
sqrt(mean(resid(Hg_LMER_NP_MASS)^2))
summary(lm(Hg_NP$VALUE_LOG ~ fitted(Hg_LMER_NP_MASS)))$adj.r.squared
LMER_boot_est_summary(Hg_LMER_NP_MASS_boot_est)

## WE
sqrt(mean(resid(Hg_LMER_WE_MASS)^2))
summary(lm(Hg_WE$VALUE_LOG ~ fitted(Hg_LMER_WE_MASS)))$adj.r.squared
LMER_boot_est_summary(Hg_LMER_WE_MASS_boot_est)

## Hg - AB ---- 

## LT
sqrt(mean(Hg_INLA_LT_MASS_fits$resid_inla^2))
MySqrt(Hg_INLA_LT_MASS$summary.hyperpar)
summary(lm(Hg_INLA_LT_MASS_fits$VALUE_LOG ~ Hg_INLA_LT_MASS_fits$inla_posterior_q50))

## NP
sqrt(mean(Hg_INLA_NP_MASS_fits$resid_inla^2))
MySqrt(Hg_INLA_NP_MASS$summary.hyperpar)
summary(lm(Hg_INLA_NP_MASS_fits$VALUE_LOG ~ Hg_INLA_NP_MASS_fits$inla_posterior_q50))

## WE
sqrt(mean(Hg_INLA_WE_MASS_fits$resid_inla^2))
MySqrt(Hg_INLA_WE_MASS$summary.hyperpar)
summary(lm(Hg_INLA_WE_MASS_fits$VALUE_LOG ~ Hg_INLA_WE_MASS_fits$inla_posterior_q50))

## Hg - MB ---- 

## LT
sqrt(mean(Hg_STAN_LT_MASS_fits$Hg_STAN_RESID^2))
summary(lm(Hg_STAN_LT_MASS_fits$VALUE_LOG ~ Hg_STAN_LT_MASS_fits$stan_posterior_q50))
sqrt(Hg_STAN_LT_RE)

sqrt(mean(Hg_STAN_NP_MASS_fits$Hg_STAN_RESID^2))
summary(lm(Hg_STAN_NP_MASS_fits$VALUE_LOG ~ Hg_STAN_NP_MASS_fits$stan_posterior_q50))
sqrt(Hg_STAN_NP_RE)

sqrt(mean(Hg_STAN_WE_MASS_fits$Hg_STAN_RESID^2))
sqrt(Hg_STAN_WE_RE)
summary(lm(Hg_STAN_WE_MASS_fits$VALUE_LOG ~ Hg_STAN_WE_MASS_fits$stan_posterior_q50))

## As - SER ----

## LT
LT_As_sub_preds_mass <- subset(As_sub_preds_mass, SPECIES_CODE == "LT")
median(LT_As_sub_preds_mass$RMSE) 
min(LT_As_sub_preds_mass$RMSE) 
max(LT_As_sub_preds_mass$RMSE) 

median(LT_As_sub_preds_mass$r2) 
min(LT_As_sub_preds_mass$r2) 
max(LT_As_sub_preds_mass$r2) 

sd(LT_As_sub_preds_mass$int)
sd(LT_As_sub_preds_mass$slp)

## NP
NP_As_sub_preds_mass <- subset(As_sub_preds_mass, SPECIES_CODE == "NP")
median(NP_As_sub_preds_mass$RMSE) 
min(NP_As_sub_preds_mass$RMSE) 
max(NP_As_sub_preds_mass$RMSE) 

median(NP_As_sub_preds_mass$r2) 
min(NP_As_sub_preds_mass$r2) 
max(NP_As_sub_preds_mass$r2) 

sd(NP_As_sub_preds_mass$int)
sd(NP_As_sub_preds_mass$slp)

## WE
WE_As_sub_preds_mass <- subset(As_sub_preds_mass, SPECIES_CODE == "WALL")
median(WE_As_sub_preds_mass$RMSE) 
min(WE_As_sub_preds_mass$RMSE) 
max(WE_As_sub_preds_mass$RMSE) 

median(WE_As_sub_preds_mass$r2) 
min(WE_As_sub_preds_mass$r2) 
max(WE_As_sub_preds_mass$r2) 

sd(WE_As_sub_preds_mass$int)
sd(WE_As_sub_preds_mass$slp)

## As - ML ----

## LT
sqrt(mean(resid(As_LMER_LT_MASS)^2))
summary(lm(As_LT$VALUE_LOG ~ fitted(As_LMER_LT_MASS)))$adj.r.squared
LMER_boot_est_summary(As_LMER_LT_MASS_boot_est)

## NP
sqrt(mean(resid(As_LMER_NP_MASS)^2))
summary(lm(As_NP$VALUE_LOG ~ fitted(As_LMER_NP_MASS)))$adj.r.squared
LMER_boot_est_summary(As_LMER_NP_MASS_boot_est)

## WE
sqrt(mean(resid(As_LMER_WE_MASS)^2))
summary(lm(As_WE$VALUE_LOG ~ fitted(As_LMER_WE_MASS)))$adj.r.squared
LMER_boot_est_summary(As_LMER_WE_MASS_boot_est)

## As - AB ---- 

## LT
sqrt(mean(As_INLA_LT_MASS_fits$resid_inla^2))
summary(lm(As_INLA_LT_MASS_fits$VALUE_LOG ~ As_INLA_LT_MASS_fits$inla_posterior_q50))$r.squared
MySqrt(As_INLA_LT_MASS$summary.hyperpar)

## NP
sqrt(mean(As_INLA_NP_MASS_fits$resid_inla^2))
summary(lm(As_INLA_NP_MASS_fits$VALUE_LOG ~ As_INLA_NP_MASS_fits$inla_posterior_q50))$r.squared
MySqrt(As_INLA_NP_MASS$summary.hyperpar)

## WE
sqrt(mean(As_INLA_WE_MASS_fits$resid_inla^2))
summary(lm(As_INLA_WE_MASS_fits$VALUE_LOG ~ As_INLA_WE_MASS_fits$inla_posterior_q50))$r.squared
MySqrt(As_INLA_WE_MASS$summary.hyperpar)

## As - MB ---- 

## LT
sqrt(mean(As_STAN_LT_MASS_fits$As_STAN_RESID^2))
summary(lm(As_STAN_LT_MASS_fits$VALUE_LOG ~ As_STAN_LT_MASS_fits$stan_posterior_q50))$r.squared
sqrt(As_STAN_LT_RE)


sqrt(mean(As_STAN_NP_MASS_fits$As_STAN_RESID^2))
summary(lm(As_STAN_NP_MASS_fits$VALUE_LOG ~ As_STAN_NP_MASS_fits$stan_posterior_q50))$r.squared
sqrt(As_STAN_NP_RE)

sqrt(mean(As_STAN_WE_MASS_fits$As_STAN_RESID^2))
summary(lm(As_STAN_WE_MASS_fits$VALUE_LOG ~ As_STAN_WE_MASS_fits$stan_posterior_q50))$r.squared
sqrt(As_STAN_WE_RE)



















LT_As_sub_preds_mass <- subset(As_sub_preds_mass, SPECIES_CODE == "081")
median(LT_As_sub_preds_mass$RMSE) 
min(LT_As_sub_preds_mass$RMSE) 
max(LT_As_sub_preds_mass$RMSE) 



As_sub_preds_mass$SPECIES_CODE

median(As_sub_preds_mass$RMSE)



As_sub_preds_mass

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

Hg_LT_pt_SER <- subset(hg_sub_preds_mass, (WATERBODY_CODE_SAMPLE_YEAR %in% Hg_LT_ind$WATERBODY_CODE_SAMPLE_YEAR) &
                         SPECIES_CODE == "081")

Hg_LT_SER_test <- merge(Hg_LT_ind, Hg_LT_pt_SER[, c("WATERBODY_CODE_SAMPLE_YEAR", "fit")])
nrow(Hg_LT_SER_test)

Hg_LT_SER_test$fit_LOG <- log(as.numeric(Hg_LT_SER_test$fit))
plot(Hg_LT_SER_test$fit_LOG, Hg_LT_SER_test$VALUE_LOG)

Hg_LT_SER_test_lm <- lm(Hg_LT_SER_test$VALUE_LOG ~ Hg_LT_SER_test$fit_LOG)
summary(Hg_LT_SER_test_lm)
sqrt(mean(Hg_LT_SER_test_lm$residuals^2))

Hg_LT_pred_list[["SER"]] <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = Hg_LT_SER_test$WATERBODY_CODE_SAMPLE_YEAR, 
                                       OBS = Hg_LT_SER_test$VALUE_LOG,
                                       PRED = Hg_LT_SER_test$fit_LOG)

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








## Sparse data test ----

## Find lake-year with most amount of samples for LT, NP, and WE for Hg and As 
## Run models randomly adding a new observation until exhausted
## Track prediction and confidence interval for all SER, ML, AB, but not B models (takes too long) 

## LT
Hg_DT <- Hg_LT

Hg_sparse_test <- function(Hg_DT){
  
  ## find sampling even with most  
  (DT_max_lake_year <- table(Hg_DT$WATERBODY_CODE_SAMPLE_YEAR))
  (DT_max_lake_year_ind <- which(DT_max_lake_year == max(DT_max_lake_year)))
  (DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR <- names(DT_max_lake_year)[DT_max_lake_year_ind[1]]) # take first
  
  ## assign random number to each observation 
  Hg_DT_sparse <- Hg_DT 
  
  Hg_DT_sparse$randint <- NA
  max_ind <- which(Hg_DT_sparse$WATERBODY_CODE_SAMPLE_YEAR == DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
  
  ## replicability
  set.seed(101)
  max_ind_ind <- sample(1:length(max_ind), length(max_ind), replace = FALSE)
  max_ind <- max_ind[max_ind_ind]
  (max_ind <- sample(max_ind, length(max_ind), replace = FALSE))
  
  ## 5 samples, 100 nsim for bootstrap and 4 cores = 282s
  ## 5 samples, 1000 nsim for bootstram and 4 cores = 
  system.time(sparse_seq <- foreach(i = 1:length(max_ind), .errorhandling = "pass") %do% {
    
    message("sample ", i, " @ ", Sys.time())
    
    loc <- which(max_ind[i] == max_ind)
    (loc <- 1:loc)
    
    (max_it <- max_ind[loc])
    
    ## SER 
    
    if(length(max_it) >= 4){
      
      Hg_SER_sub <- Hg_DT_sparse[max_it,]  
      Hg_SER_mod <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = Hg_SER_sub)
      Hg_SER_mod_pred <- predict(Hg_SER_mod, interval = "confidence", 
                                 newdata = data.frame(WEIGHT_GRAM_LOG = log(1000)))
      
      Hg_SER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    Hg_SER_mod_pred)
      
    } else {
      Hg_SER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
    }
    Hg_SER_mod_pred
    
    ## ML 
    (max_remove <- max_ind[!max_ind %in% max_it])
    
    if(length(max_remove)!=0){
    
      Hg_LMER_mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = Hg_DT_sparse[-max_remove,], 
                          control = lmerControl(optimizer = "Nelder_Mead")) 
    
    } else {
      
      Hg_LMER_mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = Hg_DT_sparse, 
                          control = lmerControl(optimizer = "Nelder_Mead")) 
    
    }
    
    ## Generate predictions for each waterbody based on the median size of DT in whole dataset 
    Hg_LMER_mod_predict <- unique(Hg_DT_sparse[max_it,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
    Hg_LMER_mod_predict$WEIGHT_GRAM_LOG <- log(1000)
    
    ## Bootstrap confidence intervals from bootMer
    Hg_LMER_boot_pred <- function(., newdata) {
      predict(., newdata = Hg_LMER_mod_predict)
    }
    
    clust <- LMER_boot_initiate(varlist = "Hg_LMER_mod_predict", nclust = 10)
    system.time(Hg_LMER_mod_boot_pred_res <- lme4::bootMer(Hg_LMER_mod, Hg_LMER_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 10)) # 150 
    
    (Hg_LMER_mod_boot_pred_res_sum <- LMER_boot_pred_summary(Hg_LMER_mod_boot_pred_res))
    
    (Hg_LMER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    Hg_LMER_mod_boot_pred_res_sum))
    
    ## INLA 
    
    ## Generate predictions for each waterbody based on the median size of DT in whole dataset 
    (Hg_INLA_mod_predict <- Hg_DT_sparse[max_it,])
    (Hg_INLA_mod_predict <- Hg_INLA_mod_predict[1,])
    (Hg_INLA_mod_predict$WEIGHT_GRAM_LOG <- log(1000))
    (Hg_INLA_mod_predict$VALUE_LOG <- NA)
    
    if(length(max_remove)!=0){
      
      Hg_INLA_mod_prediction <- rbind(Hg_DT_sparse[-max_remove,], 
                                      Hg_INLA_mod_predict)
    } else {
      
      Hg_INLA_mod_prediction <- rbind(Hg_DT_sparse, 
                                      Hg_INLA_mod_predict)
    }
    
    Hg_INLA_mod_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * Hg_DT$n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = Hg_INLA_mod_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
    )
    
    res <- Hg_INLA_mod_prediction_model$summary.fitted.values[nrow(Hg_INLA_mod_prediction_model$summary.fitted.values), c("0.5quant", "0.025quant", "0.975quant")]
    names(res) <- c("fit", "lwr", "upr")
    
    (Hg_INLA_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    res))
    
    res_list <- list(
      Hg_SER_mod_pred,
      Hg_LMER_mod_pred,
      Hg_INLA_mod_pred)
    
    return(res_list)
    
  })
  
  return(sparse_seq)
  
  
  
}

system.time(Hg_LT_sparse_test <- Hg_sparse_test(Hg_LT))
saveRDS(Hg_LT_sparse_test, "./out_workspaces/Hg_LT_sparse_test.RDS")
system.time(Hg_NP_sparse_test <- Hg_sparse_test(Hg_NP))
saveRDS(Hg_NP_sparse_test, "./out_workspaces/Hg_NP_sparse_test.RDS")
system.time(Hg_WE_sparse_test <- Hg_sparse_test(Hg_WE))
saveRDS(Hg_WE_sparse_test, "./out_workspaces/Hg_WE_sparse_test.RDS")

As_DT <- As_NP

As_sparse_test <- function(As_DT){
  
  ## find sampling even with most  
  (DT_max_lake_year <- table(As_DT$WATERBODY_CODE_SAMPLE_YEAR))
  (DT_max_lake_year_ind <- which(DT_max_lake_year == max(DT_max_lake_year)))
  (DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR <- names(DT_max_lake_year)[DT_max_lake_year_ind[1]]) # take first
  
  ## assign random number to each observation 
  As_DT_sparse <- As_DT 
  
  As_DT_sparse$randint <- NA
  max_ind <- which(As_DT_sparse$WATERBODY_CODE_SAMPLE_YEAR == DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
  
  ## replicability
  set.seed(105)
  max_ind_ind <- sample(1:length(max_ind), length(max_ind), replace = FALSE)
  max_ind <- max_ind[max_ind_ind]
  (max_ind <- sample(max_ind, length(max_ind), replace = FALSE))
  
  ## 5 samples, 100 nsim for bootstrap and 4 cores = 282s
  ## 5 samples, 1000 nsim for bootstram and 4 cores = 
 
  system.time(sparse_seq <- foreach(i = 1:length(max_ind), .errorhandling = "pass") %do% {
    
    message("sample ", i, " @ ", Sys.time())
    
    loc <- which(max_ind[i] == max_ind)
    (loc <- 1:loc)
    
    (max_it <- max_ind[loc])
    
    ## SER 
    
    if(length(max_it) >= 4){
      
      As_SER_sub <- As_DT_sparse[max_it,]  
      As_SER_mod <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = As_SER_sub)
      As_SER_mod_pred <- predict(As_SER_mod, interval = "confidence", 
                                 newdata = data.frame(WEIGHT_GRAM_LOG = log(1000)))
      
      As_SER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    As_SER_mod_pred)
      
    } else {
      As_SER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
    }
    
    ## ML 
    (max_remove <- max_ind[!max_ind %in% max_it])
    
    if(length(max_remove)!=0){
      
      As_LMER_mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_DT_sparse[-max_remove,], 
                          control = lmerControl(optimizer = "Nelder_Mead")) 
      
    } else {
      
      As_LMER_mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                            (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
                            (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_DT_sparse, 
                          control = lmerControl(optimizer = "Nelder_Mead")) 
      
    }
    
    ## Generate predictions for each waterbody based on the median size of DT in whole dataset 
    As_LMER_mod_predict <- unique(As_DT_sparse[max_it,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
    As_LMER_mod_predict$WEIGHT_GRAM_LOG <- log(1000)
    
    ## Bootstrap confidence intervals from bootMer
    As_LMER_boot_pred <- function(., newdata) {
      predict(., newdata = As_LMER_mod_predict)
    }
    
    clust <- LMER_boot_initiate(varlist = "As_LMER_mod_predict", nclust = 4)
    system.time(As_LMER_mod_boot_pred_res <- lme4::bootMer(As_LMER_mod, As_LMER_boot_pred, 
                                                           nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
                                                           parallel = "snow", 
                                                           cl = clust, 
                                                           ncpus = 4)) # 150 
    
    (As_LMER_mod_boot_pred_res_sum <- LMER_boot_pred_summary(As_LMER_mod_boot_pred_res))
    
    (As_LMER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    As_LMER_mod_boot_pred_res_sum))
    
    ## INLA 
    
    ## Generate unique waterbodies
    As_DT_n_waterbody <- dplyr::n_distinct(As_DT$WATERBODY_CODE) # for inla model
    
    ## Generate predictions for each waterbody based on the median size of DT in whole dataset 
    (As_INLA_mod_predict <- As_DT_sparse[max_it,])
    (As_INLA_mod_predict <- As_INLA_mod_predict[1,])
    (As_INLA_mod_predict$WEIGHT_GRAM_LOG <- log(1000))
    (As_INLA_mod_predict$VALUE_LOG <- NA)
    
    if(length(max_remove)!=0){
      
      As_INLA_mod_prediction <- rbind(As_DT_sparse[-max_remove,], 
                                      As_INLA_mod_predict)
    } else {
      
      As_INLA_mod_prediction <- rbind(As_DT_sparse, 
                                      As_INLA_mod_predict)
    }
    
    As_INLA_mod_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * As_DT_n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = As_INLA_mod_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
    )
    
    res <- As_INLA_mod_prediction_model$summary.fitted.values[nrow(As_INLA_mod_prediction_model$summary.fitted.values), c("0.5quant", "0.025quant", "0.975quant")]
    names(res) <- c("fit", "lwr", "upr")
    
    (As_INLA_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    res))
    
    res_list <- list(
      As_SER_mod_pred,
      As_LMER_mod_pred,
      As_INLA_mod_pred)
    
    return(res_list)
    
  })
  
  return(sparse_seq)
  
}

As_sparse_test_repair <- function(As_DT){
  
  ## find sampling even with most  
  (DT_max_lake_year <- table(As_DT$WATERBODY_CODE_SAMPLE_YEAR))
  (DT_max_lake_year_ind <- which(DT_max_lake_year == max(DT_max_lake_year)))
  (DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR <- names(DT_max_lake_year)[DT_max_lake_year_ind[1]]) # take first
  
  ## assign random number to each observation 
  As_DT_sparse <- As_DT 
  
  As_DT_sparse$randint <- NA
  max_ind <- which(As_DT_sparse$WATERBODY_CODE_SAMPLE_YEAR == DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
  
  ## replicability
  set.seed(101)
  max_ind_ind <- sample(1:length(max_ind), length(max_ind), replace = FALSE)
  max_ind <- max_ind[max_ind_ind]
  (max_ind <- sample(max_ind, length(max_ind), replace = FALSE))
  
  ## 5 samples, 100 nsim for bootstrap and 4 cores = 282s
  ## 5 samples, 1000 nsim for bootstram and 4 cores = 
  
  system.time(sparse_seq <- foreach(i = 1:length(max_ind), .errorhandling = "pass") %do% {
    
    message("sample ", i, " @ ", Sys.time())
    
    loc <- which(max_ind[i] == max_ind)
    (loc <- 1:loc)
    
    (max_it <- max_ind[loc])
    
    # ## SER 
    # 
    # if(length(max_it) >= 4){
    #   
    #   As_SER_sub <- As_DT_sparse[max_it,]  
    #   As_SER_mod <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = As_SER_sub)
    #   As_SER_mod_pred <- predict(As_SER_mod, interval = "confidence", 
    #                              newdata = data.frame(WEIGHT_GRAM_LOG = log(1000)))
    #   
    #   As_SER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
    #                                 As_SER_mod_pred)
    #   
    # } else {
    #   As_SER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
    # }
    # As_SER_mod_pred
    # 
    # ## ML 
     (max_remove <- max_ind[!max_ind %in% max_it])
    # 
    # if(length(max_remove)!=0){
    #   
    #   As_LMER_mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
    #                         (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
    #                         (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_DT_sparse[-max_remove,], 
    #                       control = lmerControl(optimizer = "Nelder_Mead")) 
    #   
    # } else {
    #   
    #   As_LMER_mod <- lmer(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
    #                         (WEIGHT_GRAM_LOG|WATERBODY_CODE) + 
    #                         (1|WATERBODY_CODE_SAMPLE_YEAR), data = As_DT_sparse, 
    #                       control = lmerControl(optimizer = "Nelder_Mead")) 
    #   
    # }
    # 
    # ## Generate predictions for each waterbody based on the median size of DT in whole dataset 
    # As_LMER_mod_predict <- unique(As_DT_sparse[max_it,c("WATERBODY_CODE", "SAMPLE_YEAR", "WATERBODY_CODE_SAMPLE_YEAR")])
    # As_LMER_mod_predict$WEIGHT_GRAM_LOG <- log(1000)
    # 
    # ## Bootstrap confidence intervals from bootMer
    # As_LMER_boot_pred <- function(., newdata) {
    #   predict(., newdata = As_LMER_mod_predict)
    # }
    # 
    # clust <- LMER_boot_initiate(varlist = "As_LMER_mod_predict", nclust = 4)
    # system.time(As_LMER_mod_boot_pred_res <- lme4::bootMer(As_LMER_mod, As_LMER_boot_pred, 
    #                                                        nsim=2000, use.u = TRUE, .progress = "txt",  ## 99 before
    #                                                        parallel = "snow", 
    #                                                        cl = clust, 
    #                                                        ncpus = 10)) # 150 
    # 
    # (As_LMER_mod_boot_pred_res_sum <- LMER_boot_pred_summary(As_LMER_mod_boot_pred_res))
    # 
    # (As_LMER_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
    #                                 As_LMER_mod_boot_pred_res_sum))
    # 
    ## INLA 
    
    ## Generate unique waterbodies
    As_DT_n_waterbody <- dplyr::n_distinct(As_DT$WATERBODY_CODE) # for inla model
    
    ## Generate predictions for each waterbody based on the median size of DT in whole dataset 
    (As_INLA_mod_predict <- As_DT_sparse[max_it,])
    (As_INLA_mod_predict <- As_INLA_mod_predict[1,])
    (As_INLA_mod_predict$WEIGHT_GRAM_LOG <- log(1000))
    (As_INLA_mod_predict$VALUE_LOG <- NA)
    
    if(length(max_remove)!=0){
      
      As_INLA_mod_prediction <- rbind(As_DT_sparse[-max_remove,], 
                                      As_INLA_mod_predict)
    } else {
      
      As_INLA_mod_prediction <- rbind(As_DT_sparse, 
                                      As_INLA_mod_predict)
    }
    
    As_INLA_mod_prediction_model <- inla(VALUE_LOG ~ WEIGHT_GRAM_LOG + 
                                           f(WATERBODY_CODE1, n = 2 * As_DT_n_waterbody, model = "iid2d") +   
                                           f(WATERBODY_CODE2, WEIGHT_GRAM_LOG, copy = "WATERBODY_CODE1") + 
                                           f(WATERBODY_CODE_SAMPLE_YEAR, model = "iid"),
                                         data = As_INLA_mod_prediction, 
                                         control.predictor = list(
                                           compute = TRUE, 
                                           quantiles = c(0.025, 0.5, 0.975)
                                         ),
                                         control.compute = list(
                                           cpo = TRUE
                                         )
    )
    
    res <- As_INLA_mod_prediction_model$summary.fitted.values[nrow(As_INLA_mod_prediction_model$summary.fitted.values), c("0.5quant", "0.025quant", "0.975quant")]
    names(res) <- c("fit", "lwr", "upr")
    
    (As_INLA_mod_pred <- data.frame(WATERBODY_CODE_SAMPLE_YEAR = DT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR, 
                                    res))
    
    As_SER_mod_pred <- NULL
    As_LMER_mod_pred <- NULL
    
    res_list <- list(
      As_SER_mod_pred,
      As_LMER_mod_pred,
      As_INLA_mod_pred)
    
    return(res_list)
    
  })
  
  return(sparse_seq)

}

As_LT_sparse_repair <- As_sparse_test_repair(As_LT)
As_NP_sparse_repair <- As_sparse_test_repair(As_NP)
As_WE_sparse_repair <- As_sparse_test_repair(As_WE)

for (i in 1:length(As_LT_sparse_test)){
  As_LT_sparse_test[[i]][[3]] <- As_LT_sparse_repair[[i]][[3]] 
  As_NP_sparse_test[[i]][[3]] <- As_NP_sparse_repair[[i]][[3]] 
  As_WE_sparse_test[[i]][[3]] <- As_WE_sparse_repair[[i]][[3]] 
}

As_LT_sparse_test <- As_sparse_test(As_LT)
saveRDS(As_LT_sparse_test, "./out_workspaces/As_LT_sparse_test.RDS")
As_NP_sparse_test <- As_sparse_test(As_NP)
saveRDS(As_NP_sparse_test, "./out_workspaces/As_NP_sparse_test.RDS")
As_WE_sparse_test <- As_sparse_test(As_WE)
saveRDS(As_WE_sparse_test, "./out_workspaces/As_WE_sparse_test.RDS")

As_NP_sparse_test
As_LT_sparse_test
As_WE_sparse_test

## 
Hg_LT_sparse_test <- readRDS("./out_workspaces/Hg_LT_sparse_test.RDS")
Hg_NP_sparse_test <- readRDS("./out_workspaces/Hg_NP_sparse_test.RDS")
Hg_WE_sparse_test <- readRDS("./out_workspaces/Hg_WE_sparse_test.RDS")

As_LT_sparse_test <- readRDS("./out_workspaces/As_LT_sparse_test.RDS")
As_NP_sparse_test <- readRDS("./out_workspaces/As_NP_sparse_test.RDS")
As_WE_sparse_test <- readRDS("./out_workspaces/As_WE_sparse_test.RDS")

plot_this <- function(){

par(mfrow = c(3,1))

## LT

Hg_LT_SER_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[1]})
Hg_LT_SER_fit <- dplyr::bind_rows(Hg_LT_SER_fit)
Hg_LT_SER_fit$n <- seq(1:nrow(Hg_LT_SER_fit))

Hg_LT_LMER_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[2]})
Hg_LT_LMER_fit <- dplyr::bind_rows(Hg_LT_LMER_fit)
Hg_LT_LMER_fit$n <- seq(1:nrow(Hg_LT_LMER_fit))

Hg_LT_INLA_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[3]})
Hg_LT_INLA_fit <- dplyr::bind_rows(Hg_LT_INLA_fit)
Hg_LT_INLA_fit$n <- seq(1:nrow(Hg_LT_INLA_fit))

with(Hg_LT_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "LT"))
with(Hg_LT_LMER_fit, arrows(n, fit, n, lwr, length = 0.025, angle = 90, col = "orange"))
with(Hg_LT_LMER_fit, arrows(n, fit, n, upr, length = 0.025, angle = 90, col = "orange"))
with(Hg_LT_LMER_fit, points(n, fit, pch = 16, col = "orange"))

with(Hg_LT_INLA_fit, arrows(n+0.25, fit, n+0.25, lwr, length = 0.025, angle = 90, col = "blue"))
with(Hg_LT_INLA_fit, arrows(n+0.25, fit, n+0.25, upr, length = 0.025, angle = 90, col = "blue"))
with(Hg_LT_INLA_fit, points(n+0.25, fit, pch = 16, col = "blue"))

with(Hg_LT_SER_fit, arrows(n-0.25, fit, n-0.25, lwr, length = 0.025, angle = 90, col = "red"))
with(Hg_LT_SER_fit, arrows(n-0.25, fit, n-0.25, upr, length = 0.025, angle = 90, col = "red"))
with(Hg_LT_SER_fit, points(n-0.25, fit, pch = 16, col = "red"))

axis(1, at = c(0:max(Hg_LT_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2)

## NP 

Hg_NP_SER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[1]})
Hg_NP_SER_fit <- dplyr::bind_rows(Hg_NP_SER_fit)
Hg_NP_SER_fit$n <- seq(1:nrow(Hg_NP_SER_fit))

Hg_NP_LMER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[2]})
Hg_NP_LMER_fit <- dplyr::bind_rows(Hg_NP_LMER_fit)
Hg_NP_LMER_fit$n <- seq(1:nrow(Hg_NP_LMER_fit))

Hg_NP_INLA_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[3]})
Hg_NP_INLA_fit <- dplyr::bind_rows(Hg_NP_INLA_fit)
Hg_NP_INLA_fit$n <- seq(1:nrow(Hg_NP_INLA_fit))

with(Hg_NP_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "NP"))
with(Hg_NP_LMER_fit, arrows(n, fit, n, lwr, length = 0.025, angle = 90, col = "orange"))
with(Hg_NP_LMER_fit, arrows(n, fit, n, upr, length = 0.025, angle = 90, col = "orange"))
with(Hg_NP_LMER_fit, points(n, fit, pch = 16, col = "orange"))

with(Hg_NP_INLA_fit, arrows(n+0.25, fit, n+0.25, lwr, length = 0.025, angle = 90, col = "blue"))
with(Hg_NP_INLA_fit, arrows(n+0.25, fit, n+0.25, upr, length = 0.025, angle = 90, col = "blue"))
with(Hg_NP_INLA_fit, points(n+0.25, fit, pch = 16, col = "blue"))

with(Hg_NP_SER_fit, arrows(n-0.25, fit, n-0.25, lwr, length = 0.025, angle = 90, col = "red"))
with(Hg_NP_SER_fit, arrows(n-0.25, fit, n-0.25, upr, length = 0.025, angle = 90, col = "red"))
with(Hg_NP_SER_fit, points(n-0.25, fit, pch = 16, col = "red"))

axis(1, at = c(0:max(Hg_NP_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2)

## WE

Hg_WE_SER_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[1]})
Hg_WE_SER_fit <- dplyr::bind_rows(Hg_WE_SER_fit)
Hg_WE_SER_fit$n <- seq(1:nrow(Hg_WE_SER_fit))

Hg_WE_LMER_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[2]})
Hg_WE_LMER_fit <- dplyr::bind_rows(Hg_WE_LMER_fit)
Hg_WE_LMER_fit$n <- seq(1:nrow(Hg_WE_LMER_fit))

Hg_WE_INLA_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[3]})
Hg_WE_INLA_fit <- dplyr::bind_rows(Hg_WE_INLA_fit)
Hg_WE_INLA_fit$n <- seq(1:nrow(Hg_WE_INLA_fit))

with(Hg_WE_LMER_fit, plot(fit ~ n, ylim = c(-6, -1), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "WE"))
with(Hg_WE_LMER_fit, arrows(n, fit, n, lwr, length = 0.025, angle = 90, col = "orange"))
with(Hg_WE_LMER_fit, arrows(n, fit, n, upr, length = 0.025, angle = 90, col = "orange"))
with(Hg_WE_LMER_fit, points(n, fit, pch = 16, col = "orange"))

with(Hg_WE_INLA_fit, arrows(n+0.25, fit, n+0.25, lwr, length = 0.025, angle = 90, col = "blue"))
with(Hg_WE_INLA_fit, arrows(n+0.25, fit, n+0.25, upr, length = 0.025, angle = 90, col = "blue"))
with(Hg_WE_INLA_fit, points(n+0.25, fit, pch = 16, col = "blue"))

with(Hg_WE_SER_fit, arrows(n-0.25, fit, n-0.25, lwr, length = 0.025, angle = 90, col = "red"))
with(Hg_WE_SER_fit, arrows(n-0.25, fit, n-0.25, upr, length = 0.025, angle = 90, col = "red"))
with(Hg_WE_SER_fit, points(n-0.25, fit, pch = 16, col = "red"))

axis(1, at = c(0:max(Hg_WE_LMER_fit$n)))
axis(2, at = seq(-6, -1, by = 0.25), las =2) 

}

## OR ---- 

Hg_LT_SER_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[1]})
Hg_LT_SER_fit <- dplyr::bind_rows(Hg_LT_SER_fit)
Hg_LT_SER_fit$n <- seq(1:nrow(Hg_LT_SER_fit))

Hg_LT_LMER_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[2]})
Hg_LT_LMER_fit <- dplyr::bind_rows(Hg_LT_LMER_fit)
Hg_LT_LMER_fit$n <- seq(1:nrow(Hg_LT_LMER_fit))

Hg_LT_INLA_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[3]})
Hg_LT_INLA_fit <- dplyr::bind_rows(Hg_LT_INLA_fit)
Hg_LT_INLA_fit$n <- seq(1:nrow(Hg_LT_INLA_fit))

Hg_NP_SER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[1]})
Hg_NP_SER_fit <- dplyr::bind_rows(Hg_NP_SER_fit)
Hg_NP_SER_fit$n <- seq(1:nrow(Hg_NP_SER_fit))

Hg_NP_LMER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[2]})
Hg_NP_LMER_fit <- dplyr::bind_rows(Hg_NP_LMER_fit)
Hg_NP_LMER_fit$n <- seq(1:nrow(Hg_NP_LMER_fit))

Hg_NP_INLA_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[3]})
Hg_NP_INLA_fit <- dplyr::bind_rows(Hg_NP_INLA_fit)
Hg_NP_INLA_fit$n <- seq(1:nrow(Hg_NP_INLA_fit))

Hg_WE_SER_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[1]})
Hg_WE_SER_fit <- dplyr::bind_rows(Hg_WE_SER_fit)
Hg_WE_SER_fit$n <- seq(1:nrow(Hg_WE_SER_fit))

Hg_WE_LMER_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[2]})
Hg_WE_LMER_fit <- dplyr::bind_rows(Hg_WE_LMER_fit)
Hg_WE_LMER_fit$n <- seq(1:nrow(Hg_WE_LMER_fit))

Hg_WE_INLA_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[3]})
Hg_WE_INLA_fit <- dplyr::bind_rows(Hg_WE_INLA_fit)
Hg_WE_INLA_fit$n <- seq(1:nrow(Hg_WE_INLA_fit))

As_LT_SER_fit <- lapply(As_LT_sparse_test, function(xx){xx[1]})
As_LT_SER_fit <- dplyr::bind_rows(As_LT_SER_fit)
As_LT_SER_fit$n <- seq(1:nrow(As_LT_SER_fit))

As_LT_LMER_fit <- lapply(As_LT_sparse_test, function(xx){xx[2]})
As_LT_LMER_fit <- dplyr::bind_rows(As_LT_LMER_fit)
As_LT_LMER_fit$n <- seq(1:nrow(As_LT_LMER_fit))

As_LT_INLA_fit <- lapply(As_LT_sparse_test, function(xx){xx[3]})
As_LT_INLA_fit <- dplyr::bind_rows(As_LT_INLA_fit)
As_LT_INLA_fit$n <- seq(1:nrow(As_LT_INLA_fit))

As_NP_SER_fit <- lapply(As_NP_sparse_test, function(xx){xx[1]})
As_NP_SER_fit <- dplyr::bind_rows(As_NP_SER_fit)
As_NP_SER_fit$n <- seq(1:nrow(As_NP_SER_fit))

As_NP_LMER_fit <- lapply(As_NP_sparse_test, function(xx){xx[2]})
As_NP_LMER_fit <- dplyr::bind_rows(As_NP_LMER_fit)
As_NP_LMER_fit$n <- seq(1:nrow(As_NP_LMER_fit))

As_NP_INLA_fit <- lapply(As_NP_sparse_test, function(xx){xx[3]})
As_NP_INLA_fit <- dplyr::bind_rows(As_NP_INLA_fit)
As_NP_INLA_fit$n <- seq(1:nrow(As_NP_INLA_fit))

As_WE_SER_fit <- lapply(As_WE_sparse_test, function(xx){xx[1]})
As_WE_SER_fit <- dplyr::bind_rows(As_WE_SER_fit)
As_WE_SER_fit$n <- seq(1:nrow(As_WE_SER_fit))

As_WE_LMER_fit <- lapply(As_WE_sparse_test, function(xx){xx[2]})
As_WE_LMER_fit <- dplyr::bind_rows(As_WE_LMER_fit)
As_WE_LMER_fit$n <- seq(1:nrow(As_WE_LMER_fit))

As_WE_INLA_fit <- lapply(As_WE_sparse_test, function(xx){xx[3]})
As_WE_INLA_fit <- dplyr::bind_rows(As_WE_INLA_fit)
As_WE_INLA_fit$n <- seq(1:nrow(As_WE_INLA_fit))

par(mfrow = c(3,2), mar = c(3,3,3,3))

## LT

## Plot 1
with(Hg_LT_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - LT"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_LT_SER_fit <- Hg_LT_SER_fit[complete.cases(Hg_LT_SER_fit),]

with(Hg_LT_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(Hg_LT_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(Hg_LT_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(Hg_LT_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(Hg_LT_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(Hg_LT_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(Hg_LT_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2, pos = 0)

## Plot 2
with(As_LT_LMER_fit, plot(fit ~ n, ylim = c(-5, -3), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - LT"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_LT_SER_fit <- As_LT_SER_fit[complete.cases(As_LT_SER_fit),]

with(As_LT_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(As_LT_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(As_LT_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(As_LT_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(As_LT_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(As_LT_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(As_LT_LMER_fit$n)))
axis(2, at = seq(-5, -3, by = 0.5), las =2, pos = 0)

## NP

## Plot 1
with(Hg_NP_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_NP_SER_fit <- Hg_NP_SER_fit[complete.cases(Hg_NP_SER_fit),]

with(Hg_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(Hg_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(Hg_NP_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(Hg_NP_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(Hg_NP_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2, pos = 0)

## Plot 2
with(As_NP_LMER_fit, plot(fit ~ n, ylim = c(-10, -1.5), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_NP_SER_fit <- As_NP_SER_fit[complete.cases(As_NP_SER_fit),]

with(As_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(As_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(As_NP_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(As_NP_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(As_NP_LMER_fit$n)))
axis(2, at = seq(-3, -1.5, by = 0.5), las =2)

summary(As_LMER_NP_MASS)
summary(As_INLA_NP_MASS)

As_INLA_NP_MASS$summary.random$WATERBODY_CODE1

## NP

with(Hg_NP_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - NP"))

Hg_NP_SER_fit <- Hg_NP_SER_fit[complete.cases(Hg_NP_SER_fit),]

with(Hg_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(Hg_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(Hg_NP_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(Hg_NP_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(Hg_NP_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2)


## WE

with(Hg_WE_LMER_fit, plot(fit ~ n, ylim = c(-6, 1), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - WE"))

Hg_WE_SER_fit <- Hg_WE_SER_fit[complete.cases(Hg_WE_SER_fit),]

with(Hg_WE_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(Hg_WE_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(Hg_WE_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(Hg_WE_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(Hg_WE_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(Hg_WE_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(Hg_WE_LMER_fit$n)))
axis(2, at = seq(-6, -1, by = 0.5), las =2)





sapply(Hg_LT_sparse_test, class)

for (i in 1:length(Hg_WE_sparse_test)){
  print(Hg_WE_sparse_test[[i]])
  Sys.sleep(2)
}

lapply(this, print)

str(this)

this[[15]]

500 * 24


## Run SER 
LT_SER_sub <- subset(Hg_LT_sparse, WATERBODY_CODE_SAMPLE_YEAR == LT_max_lake_year_WATERBODY_CODE_SAMPLE_YEAR)
LT_SER_sub$randint <- order(LT_SER_sub$randint)






## merge with 



## for i in 1:length,  

## 

## Run SER 

ser_mod <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = )

Hg_LT$WEIGHT_GRAM_LOG









Hg_LT_pt_SER





Hg_LT_pt_SER$fit

## PLOTS ---- 

## 1) Histograms ---- 

tiff("./out_figs/Fig1_HistrogramComparison.tif", 
     width = 8, height = 8, units = "in", res = 300)

par(mfrow = c(3, 2), mar = c(4,3,2,1))

h1 <- hist(as.numeric(subset(hg_sub_preds_mass, SPECIES_CODE == "081")$r2), 
           breaks = 50, ylim = c(0, 50), xlim = c(0,1), las = 2, 
           main = "Hg - LT", xlab = "", ylab = "",
           col = viridis::viridis(10)[1], border = viridis::viridis(10)[1], 
           axes = F)
axis(1, at = seq(0,1,by = 0.25), pos = 0)
axis(2, at = seq(0,50,by = 10), las = 2, pos = 0)
title(xlab = expression(R^2), ylab = "Frequency", line = 2)

legend("top", legend = c("SER", "ML", "AB", "MB"), bty = "n", lty = 1, 
       col = c("red", "orange", "blue", "black"), horiz = T)

arrows(0.54, 0, 0.54, 25, col = adjustcolor("red"), lwd = 2, length = 0)
arrows(0.83, 0, 0.83, 25, col = adjustcolor("orange"), lwd = 2, length = 0)
arrows(0.83, 0, 0.83, 25, col = adjustcolor("blue"), lwd = 2, length = 0, lty = 2)
arrows(0.83, 0, 0.83, 25, col = adjustcolor("black"), lwd = 2, length = 0, lty = 3)

h1 <- hist(as.numeric(subset(As_sub_preds_mass, SPECIES_CODE == "LT")$r2), 
           breaks = 50, ylim = c(0, 10), xlim = c(0,1), las = 2, 
           main = "As - LT", xlab = "", ylab = "",
           col = viridis::viridis(10)[5], border = viridis::viridis(10)[5], 
           axes = F)
axis(1, at = seq(0,1,by = 0.25), pos = 0)
axis(2, at = seq(0,10,by = 2), las = 2, pos = 0)
title(xlab = expression(R^2), ylab = "Frequency", line = 2)

legend("top", legend = c("SER", "ML", "AB", "MB"), bty = "n", lty = 1, 
       col = c("red", "orange", "blue", "black"), horiz = T)

arrows(0.22, 0, 0.22, 5, col = adjustcolor("red"), lwd = 2, length = 0)
arrows(0.79, 0, 0.79, 5, col = adjustcolor("orange"), lwd = 2, length = 0)
arrows(0.77, 0, 0.77, 5, col = adjustcolor("blue"), lwd = 2, length = 0, lty = 2)
arrows(0.77, 0, 0.77, 5, col = adjustcolor("black"), lwd = 2, length = 0, lty = 3)

h2 <- hist(as.numeric(subset(hg_sub_preds_mass, SPECIES_CODE == "131")$r2), 
           breaks = 50, ylim = c(0, 50), xlim = c(0,1), las = 2, 
           main = "Hg - NP", xlab = "", ylab = "",
           col = viridis::viridis(10)[1], border = viridis::viridis(10)[1], 
           axes = F)
axis(1, at = seq(0,1,by = 0.25), pos = 0)
axis(2, at = seq(0,50,by = 10), las = 2, pos = 0)
title(xlab = expression(R^2), ylab = "Frequency", line = 2)

arrows(0.60, 0, 0.60, 25, col = adjustcolor("red"), lwd = 2, length = 0)
arrows(0.81, 0, 0.81, 25, col = adjustcolor("orange"), lwd = 2, length = 0)
arrows(0.81, 0, 0.81, 25, col = adjustcolor("blue"), lwd = 2, length = 0, lty = 2)
arrows(0.81, 0, 0.81, 25, col = adjustcolor("black"), lwd = 2, length = 0, lty = 3)

h2 <- hist(as.numeric(subset(As_sub_preds_mass, SPECIES_CODE == "NP")$r2), 
           breaks = 50, ylim = c(0, 10), xlim = c(0,1), las = 2, 
           main = "As - NP", xlab = "", ylab = "",
           col = viridis::viridis(10)[5], border = viridis::viridis(10)[5], 
           axes = F)
axis(1, at = seq(0,1,by = 0.25), pos = 0)
axis(2, at = seq(0,10,by = 2), las = 2, pos = 0)
title(xlab = expression(R^2), ylab = "Frequency", line = 2)

arrows(0.34, 0, 0.34, 5, col = adjustcolor("red"), lwd = 2, length = 0)
arrows(0.86, 0, 0.86, 5, col = adjustcolor("orange"), lwd = 2, length = 0)
arrows(0.84, 0, 0.84, 5, col = adjustcolor("blue"), lwd = 2, length = 0, lty = 2)
arrows(0.85, 0, 0.85, 5, col = adjustcolor("black"), lwd = 2, length = 0, lty = 3)

h3 <- hist(as.numeric(subset(hg_sub_preds_mass, SPECIES_CODE == "334")$r2), 
           breaks = 50, ylim = c(0, 80), xlim = c(0,1), las = 2, 
           main = "Hg - WE", xlab = "", ylab = "",
           col = viridis::viridis(10)[1], border = viridis::viridis(10)[1], 
           axes = F)
axis(1, at = seq(0,1,by = 0.25), pos = 0)
axis(2, at = seq(0,80,by = 10), las = 2, pos = 0)
title(xlab = expression(R^2), ylab = "Frequency", line = 2)

arrows(0.70, 0, 0.70, 25, col = adjustcolor("red"), lwd = 2, length = 0)
arrows(0.84, 0, 0.84, 25, col = adjustcolor("orange"), lwd = 2, length = 0)
arrows(0.84, 0, 0.84, 25, col = adjustcolor("blue"), lwd = 2, length = 0, lty = 2)
arrows(0.84, 0, 0.84, 25, col = adjustcolor("black"), lwd = 2, length = 0, lty = 3)

h2 <- hist(as.numeric(subset(As_sub_preds_mass, SPECIES_CODE == "WALL")$r2), 
           breaks = 50, ylim = c(0, 10), xlim = c(0,1), las = 2, 
           main = "As - WE", xlab = "", ylab = "",
           col = viridis::viridis(10)[5], border = viridis::viridis(10)[5], 
           axes = F)
axis(1, at = seq(0,1,by = 0.25), pos = 0)
axis(2, at = seq(0,10,by = 2), las = 2, pos = 0)
title(xlab = expression(R^2), ylab = "Frequency", line = 2)

arrows(0.27, 0, 0.27, 5, col = adjustcolor("red"), lwd = 2, length = 0)
arrows(0.77, 0, 0.77, 5, col = adjustcolor("orange"), lwd = 2, length = 0)
arrows(0.76, 0, 0.76, 5, col = adjustcolor("blue"), lwd = 2, length = 0, lty = 2)
arrows(0.75, 0, 0.75, 5, col = adjustcolor("black"), lwd = 2, length = 0, lty = 3)

dev.off()

## 2) Precision Demo ----

tiff("./out_figs/GLFCSeminar_PrecisionDemo.tif", 
     width = 10, height = 5, units = "in", res = 200)

par(mfrow = c(3,2), mar = c(4, 4, 2, 2))
#par(mfrow = c(1,2), mar = c(4,5,3,0))


## Plot 1
with(Hg_LT_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - LT"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_LT_SER_fit <- Hg_LT_SER_fit[complete.cases(Hg_LT_SER_fit),]

with(Hg_LT_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(Hg_LT_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(Hg_LT_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(Hg_LT_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(Hg_LT_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(Hg_LT_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(Hg_LT_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 20, 25, 31), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2, pos = 0)
title(ylab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")


## Plot 2
with(As_LT_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - LT"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_LT_SER_fit <- As_LT_SER_fit[complete.cases(As_LT_SER_fit),]

with(As_LT_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(As_LT_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(As_LT_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(As_LT_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(As_LT_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(As_LT_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(As_LT_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 22), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2)
title(ylab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")


## Plot 1
with(Hg_NP_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_NP_SER_fit <- Hg_NP_SER_fit[complete.cases(Hg_NP_SER_fit),]

with(Hg_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(Hg_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(Hg_NP_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(Hg_NP_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(Hg_NP_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 20, 25, 31), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2, pos = 0)
title(ylab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")


## Plot 2
with(As_NP_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_NP_SER_fit <- As_NP_SER_fit[complete.cases(As_NP_SER_fit),]

with(As_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(As_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(As_NP_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(As_NP_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(As_NP_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 22), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2)
title(ylab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")

## Plot 1
with(Hg_WE_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - WE"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_WE_SER_fit <- Hg_WE_SER_fit[complete.cases(Hg_WE_SER_fit),]

with(Hg_WE_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(Hg_WE_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(Hg_WE_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(Hg_WE_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(Hg_WE_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(Hg_WE_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(Hg_WE_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 20, 25, 31), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2, pos = 0)
title(ylab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")


## Plot 2
with(As_WE_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - WE"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_WE_SER_fit <- As_WE_SER_fit[complete.cases(As_WE_SER_fit),]

with(As_WE_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(As_WE_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(As_WE_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(As_WE_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(As_WE_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(As_WE_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(As_WE_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 22), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2)
title(ylab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")
















dev.off()













h3 <- hist(as.numeric(subset(hg_sub_preds_mass, SPECIES_CODE == "334")$adj_r), 
           breaks = 100, ylim = c(0, 50), xlim = c(-0.4, 1.0), 
           main = "WE", xlab = expression(R^2), col = "lightgrey", border = "lightgrey")

abline(v = 0.54, col = adjustcolor("red"), lwd = 2)
abline(v = 0.83, col = adjustcolor("orange"), lwd = 2)
abline(v = 0.83, col = adjustcolor("blue"), lwd = 2, lty = 2)
abline(v = 0.83, col = adjustcolor("black"), lwd = 2, lty = 3)

h3 <- hist(as.numeric(subset(As_sub_preds_mass, SPECIES_CODE == "WALL")$adj_r), 
           breaks = 100, ylim = c(0, 10), xlim = c(-0.4, 1.0), 
           main = "WE", xlab = expression(R^2), col = "lightgrey", border = "lightgrey")

abline(v = 0.54, col = adjustcolor("red"), lwd = 2)
abline(v = 0.83, col = adjustcolor("orange"), lwd = 2)
abline(v = 0.83, col = adjustcolor("blue"), lwd = 2, lty = 2)
abline(v = 0.83, col = adjustcolor("black"), lwd = 2, lty = 3)
















rn <- sample(1:638, size = 15, replace = F)

this1 <- PI2[rn,]
this2 <- boot_res_Hg_LMER_LT_MASS[rn,]
this3 <- Hg_INLA_LT_MASS_predict[rn, c("0.5quant", "0.975quant", "0.025quant")]
names(this3) <- names(this1)
this4 <- Hg_STAN_LT_MASS_predictions[rn, c("Hg_STAN_0.5quant", "Hg_STAN_0.975quant", "Hg_STAN_0.025quant")]
names(this4) <- names(this1)

names(this1)
names(this2)
names(this3)
names(this4)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(this1))-0.2, this1),
                   data.frame(Predict.Method="bootMER", x=(1:nrow(this2))-0.1, this2),
                   data.frame(Predict.Method="INLA", x=(1:nrow(this3)), this3),
                   data.frame(Predict.Method="STAN", x=(1:nrow(this4))+0.1, this4))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)


## 2 -




## GLFC Seminar Plots ---- 

## 1) Hg_NP and As_NP relationships ----

tiff("./out_figs/GLFCSeminar_HgAsConcentrationWeightRelationships.tif", 
     width = 10, height = 5, units = "in", res = 200)

par(mfrow = c(1,2), mar = c(4, 4, 2, 2))

with(Hg_NP, plot(VALUE_LOG ~ WEIGHT_GRAM_LOG, ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)), 
                 pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
title(ylab = expression("[Hg]"~"("*mu*"g/g wet weight)"), xlab = "Mass (g)")
title(main = "Northern Pike - Hg ~ Mass")

mod_Hg_NP <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = Hg_NP)
mod_Hg_NP_effects <- allEffects(mod_Hg_NP)

lines(mod_Hg_NP_effects$WEIGHT_GRAM_LOG$fit ~ mod_Hg_NP_effects$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)

with(As_NP, plot(VALUE_LOG ~ WEIGHT_GRAM_LOG, ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)), 
                 pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
title(ylab = expression("[As]"~"("*mu*"g/g wet weight)"), xlab = "Mass (g)")
title(main = "Northern Pike - As ~ Mass")

mod_As_NP <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = As_NP)
mod_As_NP_effects <- allEffects(mod_As_NP)

lines(mod_As_NP_effects$WEIGHT_GRAM_LOG$fit ~ mod_As_NP_effects$WEIGHT_GRAM_LOG$x$WEIGHT_GRAM_LOG, lwd = 2)

dev.off()

## 2) SER Demo ----

plot_ser_lm <- function(x, lwd1 = 1, lwd2 = 1, plot_pred = TRUE){
  
  mod1 <- lm(VALUE_LOG ~ WEIGHT_GRAM_LOG, data = x)
  mod1_x <- seq(min(x$WEIGHT_GRAM_LOG), 
                max(x$WEIGHT_GRAM_LOG), 
                by = 0.001)
  mod1_pred1 <- predict(mod1, newdata = data.frame(WEIGHT_GRAM_LOG = mod1_x)) 
  mod1_pred1000 <- predict(mod1, newdata = data.frame(WEIGHT_GRAM_LOG = log(1000))) 
  lines(mod1_pred1 ~ mod1_x, lwd = lwd1)
  
  if(plot_pred == TRUE){
    arrows(x0 = log(1000), y0 = log(0.002), x1 = log(1000), y1 = mod1_pred1000, length = 0, lty = 2, lwd = lwd2)
    arrows(x0 = log(20), y0 = mod1_pred1000, x1 = log(1000), y1 = mod1_pred1000, length = 0, lty = 2, lwd = lwd2)
  }
  
  
}

tiff("./out_figs/GLFCSeminar_HgAsConcentrationWeightRelationships_LakeByLake.tif", 
     width = 12, height = 6, units = "in", res = 200)

set.seed(110)
Hg_lakes <- sample(Hg_NP$WATERBODY_CODE_SAMPLE_YEAR, size = 3, replace = FALSE)

Hg_lakes <- lapply(Hg_lakes, function(xx){
  subset(Hg_NP, WATERBODY_CODE_SAMPLE_YEAR == xx)
})

par(mfrow = c(2, 4))

with(Hg_lakes[[1]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
title(ylab = expression("[Hg]"~"("*mu*"g/g wet weight)"))
title(ylab = expression("[Hg]"~"("*mu*"g/g wet weight)"))
plot_ser_lm(Hg_lakes[[1]])

with(Hg_lakes[[2]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(Hg_lakes[[2]])

with(Hg_lakes[[3]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(Hg_lakes[[3]])

h2 <- hist(as.numeric(subset(hg_sub_preds_mass, SPECIES_CODE == "131")$adj_r), 
           breaks = 50, ylim = c(0, 50), xlim = c(-0.4, 1.0), main = "", 
           xlab = "", col = viridis::viridis(10)[1], border = viridis::viridis(10)[1])


set.seed(108)
As_lakes <- sample(As_NP$WATERBODY_CODE_SAMPLE_YEAR, size = 3, replace = FALSE)

As_lakes <- lapply(As_lakes, function(xx){
  subset(As_NP, WATERBODY_CODE_SAMPLE_YEAR == xx)
})

with(As_lakes[[1]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
title(ylab = expression("[As]"~"("*mu*"g/g wet weight)"))
title(ylab = expression("[As]"~"("*mu*"g/g wet weight)"))
plot_ser_lm(As_lakes[[1]])

with(As_lakes[[2]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(As_lakes[[2]])

with(As_lakes[[3]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(As_lakes[[3]])

h2 <- hist(as.numeric(subset(As_sub_preds_mass, SPECIES_CODE == "NP")$adj_r), 
           breaks = 50, ylim = c(0, 10), xlim = c(-0.4, 1.0), main = "", 
           xlab = "", col = viridis::viridis(10)[5], border = viridis::viridis(10)[5])

dev.off()

##3) Mixed model demo ----

tiff("./out_figs/GLFCSeminar_HgAsConcentrationWeightRelationships_LakeByLake_MixedModel.tif", 
     width = 12, height = 6, units = "in", res = 200)

set.seed(110)
Hg_lakes <- sample(Hg_NP$WATERBODY_CODE_SAMPLE_YEAR, size = 3, replace = FALSE)

Hg_lakes <- lapply(Hg_lakes, function(xx){
  subset(Hg_NP, WATERBODY_CODE_SAMPLE_YEAR == xx)
})

par(mfrow = c(2, 4))

set.seed(110)
Hg_lakes <- sample(Hg_NP$WATERBODY_CODE_SAMPLE_YEAR, size = 3, replace = FALSE)

Hg_lakes <- lapply(Hg_lakes, function(xx){
  subset(Hg_NP, WATERBODY_CODE_SAMPLE_YEAR == xx)
})

par(mfrow = c(2, 4))

with(Hg_lakes[[1]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
title(ylab = expression("[Hg]"~"("*mu*"g/g wet weight)"))
title(ylab = expression("[Hg]"~"("*mu*"g/g wet weight)"))
plot_ser_lm(Hg_lakes[[1]])

with(Hg_lakes[[2]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(Hg_lakes[[2]])

with(Hg_lakes[[3]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(Hg_lakes[[3]])

with(Hg_lakes[[1]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))

with(Hg_lakes[[2]], points(VALUE_LOG ~ WEIGHT_GRAM_LOG, pch = 21, col = viridis::viridis(10)[1]))
with(Hg_lakes[[3]], points(VALUE_LOG ~ WEIGHT_GRAM_LOG, pch = 21, col = viridis::viridis(10)[1]))

plot_ser_lm(Hg_lakes[[1]])
plot_ser_lm(Hg_lakes[[2]])
plot_ser_lm(Hg_lakes[[3]])

Hg_lakes_full <- dplyr::bind_rows(Hg_lakes)
plot_ser_lm(Hg_lakes_full, lwd1 = 4, lwd2 = 1, plot_pred = FALSE)

set.seed(108)
As_lakes <- sample(As_NP$WATERBODY_CODE_SAMPLE_YEAR, size = 3, replace = FALSE)

As_lakes <- lapply(As_lakes, function(xx){
  subset(As_NP, WATERBODY_CODE_SAMPLE_YEAR == xx)
})

with(As_lakes[[1]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
title(ylab = expression("[As]"~"("*mu*"g/g wet weight)"))
plot_ser_lm(As_lakes[[1]])

with(As_lakes[[2]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(As_lakes[[2]])

with(As_lakes[[3]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))
plot_ser_lm(As_lakes[[3]])

with(As_lakes[[1]], plot(VALUE_LOG ~ WEIGHT_GRAM_LOG,
                         ylim = log(c(0.002, 8)), xlim = log(c(20, 15000)),
                         pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, line = 0)
axis(1, at = log(c(0.2, 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000)), 
     labels = c("", 25, 100, 250, 500, 1000, 2500, 5000, 10000, 15000), las = 2, pos = log(0.002))

with(As_lakes[[2]], points(VALUE_LOG ~ WEIGHT_GRAM_LOG, pch = 21, col = viridis::viridis(10)[5]))
with(As_lakes[[3]], points(VALUE_LOG ~ WEIGHT_GRAM_LOG, pch = 21, col = viridis::viridis(10)[5]))

plot_ser_lm(As_lakes[[1]])
plot_ser_lm(As_lakes[[2]])
plot_ser_lm(As_lakes[[3]])

As_lakes_full <- dplyr::bind_rows(As_lakes)

plot_ser_lm(As_lakes_full, lwd1 = 4, plot_pred = FALSE)

dev.off()

## 4) Prediction Demo ----

plot_prediction_test_lm <- function(x, lwd1 = 1, lwd2 = 1, plot_pred = TRUE){
  
  mod1 <- lm(OBS ~ PRED, data = x)
  mod1_x <- seq(min(x$PRED), 
                max(x$PRED), 
                by = 0.001)
  mod1_pred1 <- predict(mod1, newdata = data.frame(PRED = mod1_x)) 
  lines(mod1_pred1 ~ mod1_x, lwd = lwd1)
  print(summary(mod1)$r.squared)
}


tiff("./out_figs/GLFCSeminar_HgAs_1000gPredictionTest.tif", 
     width = 12, height = 6, units = "in", res = 200)

windows(10,10)
par(mfrow = c(2,4))

with(Hg_NP_pred_list[["SER"]], plot(OBS ~ PRED,
                         ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                         pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(ylab = expression("Observed [Hg]"~"("*mu*"g/g wet weight)"))
title(xlab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(Hg_NP_pred_list[["SER"]])

with(Hg_NP_pred_list[["ML"]], plot(OBS ~ PRED,
                                    ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                    pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(xlab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(Hg_NP_pred_list[["ML"]])

with(Hg_NP_pred_list[["AB"]], plot(OBS ~ PRED,
                                   ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                   pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(xlab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(Hg_NP_pred_list[["AB"]])

with(Hg_NP_pred_list[["MB"]], plot(OBS ~ PRED,
                                   ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                   pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(xlab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(Hg_NP_pred_list[["MB"]])

## As 

with(As_NP_pred_list[["SER"]], plot(OBS ~ PRED,
                                    ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                    pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(2, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(ylab = expression("Observed [As]"~"("*mu*"g/g wet weight)"))
title(xlab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(As_NP_pred_list[["SER"]])

with(As_NP_pred_list[["ML"]], plot(OBS ~ PRED,
                                   ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                   pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(xlab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(As_NP_pred_list[["ML"]])

with(As_NP_pred_list[["AB"]], plot(OBS ~ PRED,
                                   ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                   pch = 21, col = viridis::viridis(10)[5], axes = F, ylab = "", xlab = ""))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(xlab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(As_NP_pred_list[["AB"]])

with(As_NP_pred_list[["MB"]], plot(OBS ~ PRED,
                                   ylim = log(c(0.002, 8)), xlim = log(c(0.002, 8)),
                                   pch = 21, col = viridis::viridis(10)[1], axes = F, ylab = "", xlab = ""))
axis(1, at = log(c(0.002, 0.02, 0.2, 1, 2, 4, 8)), labels = c(0.002, 0.02, 0.2, 1, 2, 4, 8), las = 2, pos = log(0.002))
title(xlab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
plot_prediction_test_lm(As_NP_pred_list[["MB"]])

dev.off()

## 5) Precision Demo ---- 


## NP

Hg_NP_SER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[1]})
Hg_NP_SER_fit <- dplyr::bind_rows(Hg_NP_SER_fit)
Hg_NP_SER_fit$n <- seq(1:nrow(Hg_NP_SER_fit))

Hg_NP_LMER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[2]})
Hg_NP_LMER_fit <- dplyr::bind_rows(Hg_NP_LMER_fit)
Hg_NP_LMER_fit$n <- seq(1:nrow(Hg_NP_LMER_fit))

Hg_NP_INLA_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[3]})
Hg_NP_INLA_fit <- dplyr::bind_rows(Hg_NP_INLA_fit)
Hg_NP_INLA_fit$n <- seq(1:nrow(Hg_NP_INLA_fit))

As_NP_SER_fit <- lapply(As_NP_sparse_test, function(xx){xx[1]})
As_NP_SER_fit <- dplyr::bind_rows(As_NP_SER_fit)
As_NP_SER_fit$n <- seq(1:nrow(As_NP_SER_fit))

As_NP_LMER_fit <- lapply(As_NP_sparse_test, function(xx){xx[2]})
As_NP_LMER_fit <- dplyr::bind_rows(As_NP_LMER_fit)
As_NP_LMER_fit$n <- seq(1:nrow(As_NP_LMER_fit))

As_NP_INLA_fit <- lapply(As_NP_sparse_test, function(xx){xx[3]})
As_NP_INLA_fit <- dplyr::bind_rows(As_NP_INLA_fit)
As_NP_INLA_fit$n <- seq(1:nrow(As_NP_INLA_fit))







tiff("./out_figs/GLFCSeminar_PrecisionDemo.tif", 
     width = 10, height = 5, units = "in", res = 200)

par(mfrow = c(1,2), mar = c(4, 4, 2, 2))
#par(mfrow = c(1,2), mar = c(4,5,3,0))

## Plot 1
with(Hg_NP_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_NP_SER_fit <- Hg_NP_SER_fit[complete.cases(Hg_NP_SER_fit),]

with(Hg_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(Hg_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(Hg_NP_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(Hg_NP_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(Hg_NP_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 20, 25, 31), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2, pos = 0)
title(ylab = expression("Predicted [Hg]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")


## Plot 2
with(As_NP_LMER_fit, plot(fit ~ n, ylim = log(c(0.01, 3)), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_NP_SER_fit <- As_NP_SER_fit[complete.cases(As_NP_SER_fit),]

with(As_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor(viridis::viridis(10)[1], alpha.f = 0.25), 
                            border = NA))

with(As_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[5], alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor(viridis::viridis(10)[10], alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[10], type = "b"))
with(As_NP_INLA_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[5], type = "b"))
with(As_NP_SER_fit, points(n, fit, pch = 16, col = viridis::viridis(10)[1], type = "b"))

axis(1, at = 0:max(As_NP_LMER_fit$n), pos = log(0.01), labels = FALSE)
axis(1, at = c(5, 10, 15, 22), tick = FALSE, pos = log(0.01))
axis(2, at = log(c(0.01, 0.2, 1, 2, 3)), labels = c(0.01, 0.2, 1, 2, 3), las = 2)
title(ylab = expression("Predicted [As]"~"("*mu*"g/g wet weight)"))
title(xlab = "Number of fish")

dev.off()

MuMIn::r.squaredGLMM(Hg_LMER_LT_MASS)
MuMIn::r.squaredGLMM(Hg_LMER_NP_MASS)
MuMIn::r.squaredGLMM(Hg_LMER_WE_MASS)
MuMIn::r.squaredGLMM(As_LMER_LT_MASS)
MuMIn::r.squaredGLMM(As_LMER_NP_MASS)
MuMIn::r.squaredGLMM(As_LMER_WE_MASS)


?MuMIn::r.squaredGLMM


exp(c(-4, 0, 1))

log( c(0.01, 1, 3) )

Hg_LT_SER_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[1]})
Hg_LT_SER_fit <- dplyr::bind_rows(Hg_LT_SER_fit)
Hg_LT_SER_fit$n <- seq(1:nrow(Hg_LT_SER_fit))

Hg_LT_LMER_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[2]})
Hg_LT_LMER_fit <- dplyr::bind_rows(Hg_LT_LMER_fit)
Hg_LT_LMER_fit$n <- seq(1:nrow(Hg_LT_LMER_fit))

Hg_LT_INLA_fit <- lapply(Hg_LT_sparse_test, function(xx){xx[3]})
Hg_LT_INLA_fit <- dplyr::bind_rows(Hg_LT_INLA_fit)
Hg_LT_INLA_fit$n <- seq(1:nrow(Hg_LT_INLA_fit))

Hg_NP_SER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[1]})
Hg_NP_SER_fit <- dplyr::bind_rows(Hg_NP_SER_fit)
Hg_NP_SER_fit$n <- seq(1:nrow(Hg_NP_SER_fit))

Hg_NP_LMER_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[2]})
Hg_NP_LMER_fit <- dplyr::bind_rows(Hg_NP_LMER_fit)
Hg_NP_LMER_fit$n <- seq(1:nrow(Hg_NP_LMER_fit))

Hg_NP_INLA_fit <- lapply(Hg_NP_sparse_test, function(xx){xx[3]})
Hg_NP_INLA_fit <- dplyr::bind_rows(Hg_NP_INLA_fit)
Hg_NP_INLA_fit$n <- seq(1:nrow(Hg_NP_INLA_fit))

Hg_WE_SER_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[1]})
Hg_WE_SER_fit <- dplyr::bind_rows(Hg_WE_SER_fit)
Hg_WE_SER_fit$n <- seq(1:nrow(Hg_WE_SER_fit))

Hg_WE_LMER_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[2]})
Hg_WE_LMER_fit <- dplyr::bind_rows(Hg_WE_LMER_fit)
Hg_WE_LMER_fit$n <- seq(1:nrow(Hg_WE_LMER_fit))

Hg_WE_INLA_fit <- lapply(Hg_WE_sparse_test, function(xx){xx[3]})
Hg_WE_INLA_fit <- dplyr::bind_rows(Hg_WE_INLA_fit)
Hg_WE_INLA_fit$n <- seq(1:nrow(Hg_WE_INLA_fit))

As_LT_SER_fit <- lapply(As_LT_sparse_test, function(xx){xx[1]})
As_LT_SER_fit <- dplyr::bind_rows(As_LT_SER_fit)
As_LT_SER_fit$n <- seq(1:nrow(As_LT_SER_fit))

As_LT_LMER_fit <- lapply(As_LT_sparse_test, function(xx){xx[2]})
As_LT_LMER_fit <- dplyr::bind_rows(As_LT_LMER_fit)
As_LT_LMER_fit$n <- seq(1:nrow(As_LT_LMER_fit))

As_LT_INLA_fit <- lapply(As_LT_sparse_test, function(xx){xx[3]})
As_LT_INLA_fit <- dplyr::bind_rows(As_LT_INLA_fit)
As_LT_INLA_fit$n <- seq(1:nrow(As_LT_INLA_fit))

As_NP_SER_fit <- lapply(As_NP_sparse_test, function(xx){xx[1]})
As_NP_SER_fit <- dplyr::bind_rows(As_NP_SER_fit)
As_NP_SER_fit$n <- seq(1:nrow(As_NP_SER_fit))

As_NP_LMER_fit <- lapply(As_NP_sparse_test, function(xx){xx[2]})
As_NP_LMER_fit <- dplyr::bind_rows(As_NP_LMER_fit)
As_NP_LMER_fit$n <- seq(1:nrow(As_NP_LMER_fit))

As_NP_INLA_fit <- lapply(As_NP_sparse_test, function(xx){xx[3]})
As_NP_INLA_fit <- dplyr::bind_rows(As_NP_INLA_fit)
As_NP_INLA_fit$n <- seq(1:nrow(As_NP_INLA_fit))

As_WE_SER_fit <- lapply(As_WE_sparse_test, function(xx){xx[1]})
As_WE_SER_fit <- dplyr::bind_rows(As_WE_SER_fit)
As_WE_SER_fit$n <- seq(1:nrow(As_WE_SER_fit))

As_WE_LMER_fit <- lapply(As_WE_sparse_test, function(xx){xx[2]})
As_WE_LMER_fit <- dplyr::bind_rows(As_WE_LMER_fit)
As_WE_LMER_fit$n <- seq(1:nrow(As_WE_LMER_fit))

As_WE_INLA_fit <- lapply(As_WE_sparse_test, function(xx){xx[3]})
As_WE_INLA_fit <- dplyr::bind_rows(As_WE_INLA_fit)
As_WE_INLA_fit$n <- seq(1:nrow(As_WE_INLA_fit))

par(mfrow = c(3,2), mar = c(3,3,3,3))

## LT

## Plot 1
with(Hg_LT_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - LT"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_LT_SER_fit <- Hg_LT_SER_fit[complete.cases(Hg_LT_SER_fit),]

with(Hg_LT_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(Hg_LT_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(Hg_LT_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(Hg_LT_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(Hg_LT_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(Hg_LT_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(Hg_LT_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2, pos = 0)

## Plot 2
with(As_LT_LMER_fit, plot(fit ~ n, ylim = c(-5, -3), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - LT"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_LT_SER_fit <- As_LT_SER_fit[complete.cases(As_LT_SER_fit),]

with(As_LT_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(As_LT_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(As_LT_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(As_LT_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(As_LT_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(As_LT_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(As_LT_LMER_fit$n)))
axis(2, at = seq(-5, -3, by = 0.5), las =2, pos = 0)

## NP

## Plot 1
with(Hg_NP_LMER_fit, plot(fit ~ n, ylim = c(-1.25, 0.25), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "Hg - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

Hg_NP_SER_fit <- Hg_NP_SER_fit[complete.cases(Hg_NP_SER_fit),]

with(Hg_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(Hg_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(Hg_NP_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(Hg_NP_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(Hg_NP_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(Hg_NP_LMER_fit$n)))
axis(2, at = seq(-1.5, 0.25, by = 0.25), las =2, pos = 0)

## Plot 2
with(As_NP_LMER_fit, plot(fit ~ n, ylim = c(-10, -1.5), pch = 16, axes = F, 
                          ylab = "", xlab = "", type = "n", main = "As - NP"))

legend("top", legend = c("SER", "ML", "AB"), 
       fill = c(adjustcolor("red", alpha.f = 0.25), 
                adjustcolor("orange", alpha.f = 0.25), 
                adjustcolor("blue", alpha.f = 0.25)),
       border = NA,
       bty = "n", horiz = T)

As_NP_SER_fit <- As_NP_SER_fit[complete.cases(As_NP_SER_fit),]

with(As_NP_SER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                            col = adjustcolor("red", alpha.f = 0.25), 
                            border = NA))

with(As_NP_INLA_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("blue", alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, polygon(x = c(n, rev(n)), y = c(upr, rev(lwr)), 
                             col = adjustcolor("orange", alpha.f = 0.25), 
                             border = NA))

with(As_NP_LMER_fit, points(n, fit, pch = 16, col = "orange", type = "b"))
with(As_NP_INLA_fit, points(n, fit, pch = 16, col = "blue", type = "b"))
with(As_NP_SER_fit, points(n, fit, pch = 16, col = "red", type = "b"))

axis(1, at = c(0:max(As_NP_LMER_fit$n)))
axis(2, at = seq(-3, -1.5, by = 0.5), las =2)

