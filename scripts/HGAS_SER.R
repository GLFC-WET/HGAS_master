## HGAS_SER
## Author(s): Brian Kielstra
## Originated: 2022-06-04
##
##
## Runs sampling event regressions
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

## *************
## ANALYSIS ----
## *************

## *******
## HG ----
## *******

Hg_dat <- readRDS("./out_workspaces/HGAS_Hg_LTNPWEData.rds")
Hg_dat_full <- dplyr::bind_rows(Hg_dat)

Hg_dat_full_SER <- lapply(unique(Hg_dat_full$WATERBODY_CODE), function(x) {
  
  ## Subsets the waterbody of interest
  sub_waterbody <- subset(Hg_dat_full, WATERBODY_CODE == x)
  
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
Hg_dat_full_SER <- Hg_dat_full_SER[!sapply(Hg_dat_full_SER, is.null)]
Hg_dat_full_SER_results <- dplyr::bind_rows(Hg_dat_full_SER)

saveRDS(Hg_dat_full_SER_results, "./out_workspaces/HGAS_Hg_dat_full_SER_results.rds")

## *******
## As ----
## *******

As_dat <- readRDS("./out_workspaces/HGAS_As_LTNPWEData.rds")
As_dat_full <- dplyr::bind_rows(As_dat)

As_dat_full_SER <- lapply(unique(As_dat_full$WATERBODY_CODE), function(x) {
  
  ## Subsets the waterbody of interest
  sub_waterbody <- subset(As_dat_full, WATERBODY_CODE == x)
  
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
As_dat_full_SER <- As_dat_full_SER[!sapply(As_dat_full_SER, is.null)]
As_dat_full_SER_results <- dplyr::bind_rows(As_dat_full_SER)

saveRDS(As_dat_full_SER_results, "./out_workspaces/HGAS_As_dat_full_SER_results.rds")