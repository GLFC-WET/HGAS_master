## ARS_FishAs_mixedmodels_inla.R
## Author(s): Brian Kielstra
## Originated: 2021-0
##
##
## Runs regression to standardize As concentration per lake
##
##
## Script information:
## Major comments with two hashtags (procedural)
## Minor comments with one hashtag (functional)
## Major code chunks delineated by:
##    ##***********
##    ## NAME"----" (delete quotes to be navigable in R Studio Document outline)
##    ## **********
## Code lines should not exceed 80 characters where possible

## Mixed model code taken from 
## https://rstudio-pubs-static.s3.amazonaws.com/93118_c8d057fc0f754e7dafbd61c17ac5708c.html

## **********
## SETUP ----
## **********

## Load libraries
library(INLA)
library(dplyr)
library(readxl)
library(glmmTMB)
library(effects)
library(sjPlot)

#install.packages("INLA", repos = c(getOption("repos"),
#                                   INLA = "https://inla.r-inla-download.org/R/testing"
#), dep = TRUE)

## *******************
## ANALYSIS SETUP ----
## *******************

## Read the saved file
fish <- read.csv("./data/CK_Total_As_23Jun21.csv")
fish_lake <- subset(fish, System_Type == "Lake")

## Some data cleanup and variable creation 
fish_lake <- dplyr::rename(fish_lake, As_ugg = Arsenic..Âµg.L.)
fish_lake$Sampling_Year <- as.factor(as.numeric(sapply(strsplit(fish_lake$Sampling_Date, split = "-"), function(x){x[[1]]})))
fish_lake$As_ugg_log <- log(fish_lake$As_ugg)
fish_lake$RWT_log <- log(fish_lake$RWT)
fish_lake$TLEN_log <- log(fish_lake$TLEN)
fish_lake$Waterbody_Sampling_Year <- with(fish_lake, paste(Waterbody, Sampling_Year, sep = "_"))

## A little bit of data exploration 
par(mfrow = c(1,2))
hist(fish_lake$As_ugg)
dotchart(fish_lake$As_ugg) # looks like there is a sequence of outliers

hist(log(fish_lake$As_ugg))
dotchart(log(fish_lake$As_ugg)) # a few outliers but doesn't seem crazy 

hist(fish_lake$TLEN)
dotchart(fish_lake$TLEN) #>1500 m fish_lake make sense?

hist(log(fish_lake$TLEN))
dotchart(log(fish_lake$TLEN)) #>1500 m fish_lake make sense?

hist(log(fish_lake$RWT))
dotchart(log(fish_lake$RWT)) #>1500 m fish_lake make sense?

hist(fish_lake$Sampling_Year)

## How often are lakes sampled?
WB_Sampling_Year <- unique(data.frame(fish_lake$Waterbody, fish_lake$Sampling_Year))
hist(table(WB_Sampling_Year$fish_lake.Waterbody)) 
## Most lakes sampled once and a handful sampled more than that 
## You could probably safely ignore the year effects but it might be worth taking the most recent sample or something like that? 
## Here, I include Waterbody_Sampling_Year as a random effect, for better or worse

## Generate species median lengths and mass per taxon across whole dataset

median_length <- sapply(unique(fish_lake$Taxon), function(x) {
  fish_lake_sub <- subset(fish_lake, Taxon == x)
  res <- median(fish_lake_sub$TLEN, na.rm = T)
  return(res)
}, USE.NAMES = TRUE)

median_mass <- sapply(unique(fish_lake$Taxon), function(x) {
  fish_lake_sub <- subset(fish_lake, Taxon == x)
  res <- median(fish_lake_sub$RWT, na.rm = TRUE)
  return(res)
}, USE.NAMES = TRUE)

median_frame <- data.frame(
  Taxon = unique(fish_lake$Taxon),
  TLEN = median_length,
  RWT = median_mass,
  stringsAsFactors = FALSE
)

## Backbone model 
names(fish)

with(fish_lake, plot(log(As_ugg) ~ log(TLEN))) ## geeewwww.....but doesn't mean there aren't lake-specific effects
with(fish_lake, plot(log(As_ugg) ~ log(RWT))) ## geeewww....but doesn't mean there aren't lake specific effects 

## Walleye
WE <- subset(fish_lake, Taxon == "WALL") 
WE$WATERBODY_CODE <- factor(WE$Waterbody) # for inla model
WE$WATERBODY_CODE1 <- as.integer(WE$WATERBODY_CODE) # for inla model
WE$WATERBODY_CODE2 <- WE$WATERBODY_CODE1 + max(WE$WATERBODY_CODE1) # for inla model
WE_n_waterbody <- dplyr::n_distinct(WE$WATERBODY_CODE) # for inla model

## *****************************************************
## ANALYSIS 1 - SAMPLING EVENT BASED REGRESSION, LM ----
## ***************************************************** 

## Unique Waterbody, Year, and Taxon combinations
(u_wbsp <- unique(fish_lake[, c("Waterbody", "Sampling_Year", "Taxon")]))

## Combine with median_frame 
u_wbsp <- merge(u_wbsp, median_frame, by = "Taxon")
u_wbsp$As_ugg <- NA
u_wbsp$Waterbody_Sampling_Year <- paste(u_wbsp$Waterbody, u_wbsp$Sampling_Year, sep = "_")
u_wbsp$Taxon <- as.character(u_wbsp$Taxon)

## Standardize by length
## For each waterbody, for each year, for each species, generate a As~TLEN regression and prediction to the median size fish across the whole dataset for those waterbodies having more than 5 individuals
no_samples <- 5

SER_TLEN <- lapply(unique(fish_lake$Waterbody), function(x) {
  
  ## Subsets the waterbody of interest
  (sub_waterbody <- subset(fish_lake, Waterbody == x))
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Sampling_Year), function(y) {
    
    ## Subsets a year within a given waterbody
    (sub_year <- subset(sub_waterbody, Sampling_Year == y))
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Taxon), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Taxon == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      if (nrow(sub_species) > no_samples) {
        
        ## log mercury concentrations and lengths (converted to mm)
        logm <- log(sub_species$As_ugg)
        logl <- log(sub_species$TLEN)
        
        percentile <- ecdf(logl)
        
        ## regression relationship
        rel <- lm(logm ~ logl)
        
        # plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)
        
        ## Generating context and regression summary statistics
        (Waterbody <- x) # waterbody
        (Sampling_Year <- y) # year
        (Taxon <- z) # species
        (n <- length(logl)) # number of individuals
        (int <- formatC(rel$coefficients[1], digits = 4, format = "f")) # estimated intercept
        (slp <- formatC(rel$coefficients[2], digits = 4, format = "f")) # estimated slope
        (int_confint <- paste(formatC(confint(rel)[1, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated intercept confidence interval
        (slp_confint <- paste(formatC(confint(rel)[2, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated slope confidence interval
        (adj_r2 <- formatC(summary(rel)$adj.r.squared, digits = 4, format = "f")) # adjusted R2
        
        ## Modify the predicted standardized length based on species
        pred_modify <- median_frame[median_frame$Taxon == z, "TLEN"]
        
        (target_size_percentile <- percentile(pred_modify))
        
        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
                             newdata = data.frame(logl = log(pred_modify)),
                             interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        (pred <- as.data.frame(pred))
        (pred <- dplyr::rename(pred, As_ugg_pred = fit, As_ugg_pred_lwr = lwr, As_ugg_pred_upr = upr))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Waterbody, Sampling_Year, Taxon, n,
                             int, int_confint, slp, slp_confint, adj_r2,
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
SER_TLEN <- SER_length[!sapply(SER_TLEN, is.null)]
SER_TLEN <- do.call(plyr::rbind.fill, SER_TLEN)

nrow(SER_TLEN) # 90 lakes for no_samples = 5; 24 lakes for no_samples = 10

## Standardize by RWT
## For each waterbody, for each year, for each species, generate a As~TLEN regression and prediction to the median size fish across the whole dataset for those waterbodies having more than 5 individuals
SER_RWT <- lapply(unique(fish_lake$Waterbody), function(x) {
  
  ## Subsets the waterbody of interest
  (sub_waterbody <- subset(fish_lake, Waterbody == x))
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Sampling_Year), function(y) {
    
    ## Subsets a year within a given waterbody
    (sub_year <- subset(sub_waterbody, Sampling_Year == y))
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Taxon), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Taxon == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      if (nrow(sub_species) > no_samples) {
        
        ## log mercury concentrations and lengths (converted to mm)
        logm <- log(sub_species$As_ugg)
        logl <- log(sub_species$RWT)
        
        percentile <- ecdf(logl)
        
        ## regression relationship
        rel <- lm(logm ~ logl)
        
        # plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)
        
        ## Generating context and regression summary statistics
        (Waterbody <- x) # waterbody
        (Sampling_Year <- y) # year
        (Taxon <- z) # species
        (n <- length(logl)) # number of individuals
        (int <- formatC(rel$coefficients[1], digits = 4, format = "f")) # estimated intercept
        (slp <- formatC(rel$coefficients[2], digits = 4, format = "f")) # estimated slope
        (int_confint <- paste(formatC(confint(rel)[1, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated intercept confidence interval
        (slp_confint <- paste(formatC(confint(rel)[2, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated slope confidence interval
        (adj_r2 <- formatC(summary(rel)$adj.r.squared, digits = 4, format = "f")) # adjusted R2
        
        ## Modify the predicted standardized length based on species
        pred_modify <- median_frame[median_frame$Taxon == z, "RWT"]
        
        (target_size_percentile <- percentile(pred_modify))
        
        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
                             newdata = data.frame(logl = log(pred_modify)),
                             interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        (pred <- as.data.frame(pred))
        (pred <- dplyr::rename(pred, As_ugg_pred = fit, As_ugg_pred_lwr = lwr, As_ugg_pred_upr = upr))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Waterbody, Sampling_Year, Taxon, n,
                             int, int_confint, slp, slp_confint, adj_r2,
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
SER_RWT <- SER_RWT[!sapply(SER_RWT, is.null)]
SER_RWT <- do.call(plyr::rbind.fill, SER_RWT)

nrow(SER_RWT) # 90 lakes for no_samples = 5; 24 lakes for no_samples = 10

## *****************************************************************
## ANALYSIS 2 - MIXED MODEL REGRESSION, GLMMTMB (~= LME4::LMER) ---- 
## *****************************************************************

## glmmTMB, frequentist version (i.e., no confidence on random effects, predictions do not take into account error in random effects) 
## https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions
## "Note that none of the following approaches takes the uncertainty of the random effects parameters into account … if you want to take RE parameter uncertainty into account, a Bayesian approach is probably the easiest way to do it."

TMB_WE_TLEN <- glmmTMB(As_ugg_log ~ TLEN_log + (TLEN_log|Waterbody) + (1|Waterbody_Sampling_Year), data = WE) 
summary(TMB_WE_TLEN)

TMB_WE_RWT <- glmmTMB(As_ugg_log ~ RWT_log + (RWT_log|Waterbody) + (1|Waterbody_Sampling_Year), data = WE) 
summary(TMB_WE_RWT)

hist(resid(TMB_WE_TLEN)) # seems ok
hist(resid(TMB_WE_RWT)) # seems ok

plot(resid(TMB_WE_TLEN) ~ fitted(TMB_WE_TLEN)) # seems reasonable; a few very high residuals 
plot(resid(TMB_WE_RWT) ~ fitted(TMB_WE_RWT)) # seems reasonable; a few very high residuals 

## Generate predictions for each waterbody based on the median size of WE in whole dataset 

head(u_wbsp)

WE_predict <- unique(WE[,c("Waterbody", "Sampling_Year", "Waterbody_Sampling_Year")])
WE_predict$TLEN_log <- median(WE_complete$TLEN_log)
WE_predict$RWT_log <- median(WE_complete$RWT_log)
WE_predict <- WE_predict[complete.cases(WE_predict),]

TMB_WE_predictions_TLEN <- predict(TMB_WE_TLEN, newdata = WE_predict, se.fit = TRUE)
TMB_WE_predictions_RWT <- predict(TMB_WE_RWT, newdata = WE_predict, se.fit = TRUE)

WE_predict$TMB_WE_predictions_TLEN <- exp(TMB_WE_predictions_TLEN$fit)
WE_predict$TMB_WE_predictions_RWT <- exp(TMB_WE_predictions_RWT$fit)

## Both producing similar predictions? (not really expected but we'll see)
plot(WE_predictions_TLEN$fit, WE_predictions_RWT$fit) ## nearly identical 

## Useful diagnostic plots
## https://cran.r-project.org/web/packages/glmmTMB/vignettes/model_evaluation.pdf

## Easy to generate effects plots with glmmTMB 

ae_TLEN <- allEffects(TMB_WE_TLEN)
ae_RWT <- allEffects(TMB_WE_RWT)

par(mfrow = c(1,2))
plot(ae_TLEN$TLEN_log$data$As_ugg_log ~ ae_TLEN$TLEN_log$data$TLEN_log, main = "TLEN")  
lines(ae_TLEN$TLEN_log$fit ~ ae_TLEN$TLEN_log$x$TLEN_log, lwd = 2)
lines(ae_TLEN$TLEN_log$upper ~ ae_TLEN$TLEN_log$x$TLEN_log, lwd = 2, lty = 2)
lines(ae_TLEN$TLEN_log$lower ~ ae_TLEN$TLEN_log$x$TLEN_log, lwd = 2, lty = 2)

plot(ae_RWT$RWT_log$data$As_ugg_log ~ ae_RWT$RWT_log$data$RWT_log, main = "RWT")  
lines(ae_RWT$RWT_log$fit ~ ae_RWT$RWT_log$x$RWT_log, lwd = 2)
lines(ae_RWT$RWT_log$upper ~ ae_RWT$RWT_log$x$RWT_log, lwd = 2, lty = 2)
lines(ae_RWT$RWT_log$lower ~ ae_RWT$RWT_log$x$RWT_log, lwd = 2, lty = 2)

(pp <- plot_model(TMB_WE_RWT,type=c("pred"),
                  terms=c("RWT_log","Waterbody"), pred.type="re", show.legend = FALSE, show.values = TRUE, ci.lvl = NA))
## Most of the lines follow that line and some are flat. A few are steeper and a few might be reversed. 
## Overall, the characteristics of a weak relationship. But at least we are accounting for between-lake differences

## **********************************************
## ANALYSIS 3 - MIXED MODEL REGRESSION, INLA ----
## ********************************************** 

## Slightly more complicated to set up but very similar results and better estimates on the variance components 
## which I feel is an important diagnostic 

## Set up precision --> standard deviation formula; Bayesian models use precision (tau) where sd = 1/sqrt(tau) 
MySqrt <- function(x) {
  1 / sqrt(x)
}

## Only select those variables in u_wbsp

## Run model, same model as GLMMTMB as above
INLA_WE_TLEN <- inla(As_ugg_log ~ TLEN_log + 
                  f(WATERBODY_CODE1, n = 2 * WE_n_waterbody, model = "iid2d") +   
                  f(WATERBODY_CODE2, TLEN_log, copy = "WATERBODY_CODE1") + 
                  f(Waterbody_Sampling_Year, model = "iid"),
                data = WE, 
                control.predictor = list(
                  compute = TRUE, 
                  quantiles = c(0.025, 0.5, 0.975)
                ),
                control.compute = list(
                  cpo = TRUE
                )
)

summary(INLA_WE_TLEN)

## There is a lot of information in the INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
INLA_WE_TLEN_fits <- data.frame(WATERBODY_CODE = WE$WATERBODY_CODE, 
                      As_ugg_log = WE$As_ugg_log)

## Alternative comparison to fitted and residuals
INLA_WE_TLEN_fits$inla_posterior_q50 <- INLA_WE_TLEN$summary.fitted.values[, "0.5quant"]
INLA_WE_TLEN_fits$inla_posterior_q2p5 <- INLA_WE_TLEN$summary.fitted.values[, "0.025quant"]
INLA_WE_TLEN_fits$inla_posterior_q97p5 <- INLA_WE_TLEN$summary.fitted.values[, "0.975quant"]
INLA_WE_TLEN_fits$resid_inla <- INLA_WE_TLEN_fits$As_ugg_log - INLA_WE_TLEN_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = INLA_WE_TLEN_fits)
j <- order(INLA_WE_TLEN_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = INLA_WE_TLEN_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ INLA_WE_TLEN_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(INLA_WE_TLEN_fits$As_ugg_log ~ INLA_WE_TLEN_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(INLA_WE_TLEN_fits$As_ugg_log ~ INLA_WE_TLEN_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(INLA_WE_TLEN_fits$resid_inla) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
INLA_WE_TLEN_fixed <- data.frame(
  ID = rownames(INLA_WE_TLEN$summary.fixed),
  INLA_WE_TLEN$summary.fixed, stringsAsFactors = FALSE
)
names(INLA_WE_TLEN_fixed) <- c("ID", names(INLA_WE_TLEN$summary.fixed))
INLA_WE_TLEN_fixed$Type <- "Fixed"
head(INLA_WE_TLEN_fixed)

INLA_WE_TLEN_random_intercept <- INLA_WE_TLEN$summary.random$WATERBODY_CODE1
INLA_WE_TLEN_random_intercept$ID <- as.character(INLA_WE_TLEN_random_intercept$ID)
INLA_WE_TLEN_random_intercept$WATERBODY_CODE1 <- INLA_WE_TLEN_random_intercept$ID
INLA_WE_TLEN_random_intercept <- merge(INLA_WE_TLEN_random_intercept, WE[,c("WATERBODY_CODE1", "Waterbody")], no.dups = TRUE)
INLA_WE_TLEN_random_intercept <- INLA_WE_TLEN_random_intercept[!duplicated(INLA_WE_TLEN_random_intercept),]
INLA_WE_TLEN_random_intercept$Type <- "Random Intercept - Waterbody"
head(INLA_WE_TLEN_random_intercept)

INLA_WE_TLEN_random_slope <- INLA_WE_TLEN$summary.random$WATERBODY_CODE2
INLA_WE_TLEN_random_slope$ID <- as.character(INLA_WE_TLEN_random_slope$ID)
INLA_WE_TLEN_random_slope$WATERBODY_CODE2 <- INLA_WE_TLEN_random_slope$ID
INLA_WE_TLEN_random_slope <- merge(INLA_WE_TLEN_random_slope, WE[,c("WATERBODY_CODE2", "Waterbody")])
INLA_WE_TLEN_random_slope <- INLA_WE_TLEN_random_slope[!duplicated(INLA_WE_TLEN_random_slope),]
INLA_WE_TLEN_random_slope$Type <- "Random Slope - Waterbody"
head(INLA_WE_TLEN_random_slope)

INLA_WE_TLEN_random_Waterbody_Sampling_Year <- INLA_WE_TLEN$summary.random$Waterbody_Sampling_Year
INLA_WE_TLEN_random_Waterbody_Sampling_Year$Waterbody_Sampling_Year <- INLA_WE_TLEN_random_Waterbody_Sampling_Year$ID
INLA_WE_TLEN_random_Waterbody_Sampling_Year <- INLA_WE_TLEN_random_Waterbody_Sampling_Year[!duplicated(INLA_WE_TLEN_random_Waterbody_Sampling_Year),]
INLA_WE_TLEN_random_Waterbody_Sampling_Year$Type <- "Random Intercept - Waterbody Sampling Year "
head(INLA_WE_TLEN_random_Waterbody_Sampling_Year)

INLA_WE_TLEN_summary <- list(
  INLA_WE_TLEN_fixed,
  INLA_WE_TLEN_random_intercept,
  INLA_WE_TLEN_random_slope,
  INLA_WE_TLEN_random_Waterbody_Sampling_Year
)

INLA_WE_TLEN_summary <- dplyr::bind_rows(INLA_WE_TLEN_summary)
head(INLA_WE_TLEN_summary)
tail(INLA_WE_TLEN_summary)

## The estimated coefficient for the grand intercept is 
## -5.70 (-6.81, -4.68) using the 50th percentile from the posterior distribution and the 0.025 and 0.975 quantiles for the error 

## The estimated slope coefficient is 
## 0.39 (0.22, 0.57) using the same values 

## Using Type, you can hone in on individual estimate parameters for each lake 

subset(INLA_WE_TLEN_summary, Waterbody == "Agnew Lake")
## The random intercept coefficient is 0.021 (-5.286, 5.299); using the same values
## The random slope coefficent is 0.158 (-0.71, 1.04); using the same values 

## These suggest that Agnew Lake does not significantly differ by intercept or slope from the grand relationship since their values overlap 0.

## Variance parameters, difficult to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
INLA_WE_TLEN$marginals.hyperpar$`Precision for the Gaussian observations`
INLA_WE_TLEN$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
INLA_WE_TLEN$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
INLA_WE_TLEN$marginals.hyperpar$`Precision for Sampling_Year`

tau_WATERBODY_CODE1 <- INLA_WE_TLEN$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- INLA_WE_TLEN$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_Waterbody_Sampling_Year <- INLA_WE_TLEN$marginals.hyperpar$`Precision for Waterbody_Sampling_Year`
tau_residual <- INLA_WE_TLEN$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_Waterbody_Sampling_Year <- inla.emarginal(MySqrt, tau_Waterbody_Sampling_Year))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(WE_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_Waterbody_Sampling_Year)^2)

## TMB vs. INLA, comapring sigmas (i.e., Std.Dev. from TMB; below)
c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE, sigma_Waterbody_Sampling_Year)
summary(TMB_WE_TLEN)

## Random intercept: 3.4148 vs. 2.6115 
## Random slope: 0.5636 vs. 0.4459 
## Random intercept waterbody sampling year: 0.1934 vs. 0.1290 
## Residual: 0.3943 vs. 0.3810
## 

## This is admittedly more different than my Hg relationships but it still preserves the variability 
## Waterbody absolute variation >> Waterbody As~Length relation > Waterbody Sampling Year > residual variation 

## Now predictions 

## Predictions in INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the As_ugg values are set to NA
## Then you extract the posterior information 

## Take TMB predictions but fit the Waterbody Codes to the dataframe 
head(WE_predict) 

WE_predict_inla <- merge(WE_predict, WE[,c("Waterbody_Sampling_Year", "WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
WE_predict_inla <- WE_predict_inla[!duplicated(WE_predict_inla),]
WE_predict_inla$As_ugg_log <- NA
nrow(WE_predict)
nrow(WE_predict_inla)

## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(WE_predict$Waterbody_Sampling_Year, WE_predict_inla$Waterbody_Sampling_Year)

WE_inla_prediction <- WE[names(WE_predict_inla)]
WE_inla_prediction <- rbind(WE_inla_prediction, WE_predict_inla)

head(WE_inla_prediction_model)
tail(WE_inla_prediction_model)

INLA_WE_TLEN_prediction_model <- inla(As_ugg_log ~ TLEN_log + 
                       f(WATERBODY_CODE1, n = 2 * WE_n_waterbody, model = "iid2d") +   
                       f(WATERBODY_CODE2, TLEN_log, copy = "WATERBODY_CODE1") + 
                       f(Waterbody_Sampling_Year, model = "iid"),
                     data = WE_inla_prediction, 
                     control.predictor = list(
                       compute = TRUE, 
                       quantiles = c(0.025, 0.5, 0.975)
                     ),
                     control.compute = list(
                       cpo = TRUE
                     )
)

INLA_WE_TLEN_predict_posteriors <- INLA_WE_TLEN_prediction_model$summary.fitted.values[
  (nrow(WE) + 1):nrow(WE_inla_prediction),
  c("0.025quant", "0.5quant", "0.975quant")
]

nrow(INLA_WE_TLEN_predict_posteriors)

INLA_WE_TLEN_predict_posteriors <- exp(INLA_WE_TLEN_predict_posteriors)

WE_predict_inla <- cbind(WE_predict_inla, INLA_WE_TLEN_predict_posteriors)
WE_predict_inla <- WE_predict_inla[!duplicated(WE_predict_inla),]

## Does this condition still hold true for comparison?
identical(WE_predict_inla$Waterbody_Sampling_Year, WE_predict$Waterbody_Sampling_Year )

## Produces very similar results but less so at the high end - expected since both models don't seem to fit well. Still, these are very close!!
plot(WE_predict$TMB_WE_predictions_TLEN ~ WE_predict_inla$`0.5quant`)

## WE_predict_inla would feed into any waterbody-level regression 

## Overall: 

## INLA produces comparable results to the GLMMTMB model if you are going to go the mixed model route 
## INLA takes into account uncertainty in the random effects using a Bayesian approximating approach 
## INLA allows for an information gain for lakes that haven't been sampled as well ... but ... the relationship isn't as strong between As ~ TLEN which makes inference more difficult 

## ******************************************************
## COMPARING SAMPLING EVENT BASED REGRESSION TO INLA ----
## ****************************************************** 

SER_TLEN_WE <- subset(SER_TLEN, Taxon == "WALL")
nrow(SER_TLEN_WE) # 19 walleye lakes met the 5 individual criteria

SER_TLEN_WE$Waterbody_Sampling_Year <- paste(SER_TLEN_WE$Waterbody, SER_TLEN_WE$Sampling_Year, sep = "_")

SER_TLEN_WE <- merge(SER_TLEN_WE, WE_predict_inla[,c("Waterbody_Sampling_Year", "0.025quant", "0.5quant", "0.975quant")])

with(SER_TLEN_WE, plot(As_ugg_pred ~ `0.5quant`, pch = 16, col = "black"))

## I guess you would want to decide which one you "trust" more.