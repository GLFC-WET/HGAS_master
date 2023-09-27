## Merc_FishMercStandardizeTesting
## Author(s): Brian Kielstra
## Originated: 2019-11-19
##
##
## Runs regression to standardize mercury concentration per lake
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
## Setup ----
## **********

## Generate directories
proj_dir <- "F:/Project_CEON/MERC"
sp_dir <- "D:/CEON"

## Load libraries
library(sf)
library(gdata)
library(sp)
library(rstanarm)
library(INLA)
library(viridis)
library(classInt)
library(dplyr)
library(quantregForest)

install.packages("INLA", repos = c(getOption("repos"),
  INLA = "https://inla.r-inla-download.org/R/testing"
), dep = TRUE)

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

## cleaned and matched data from Merc_FishHgNetLinkage.R
st_layers("./out_spatial/Merc.gpkg")

m1 <- st_read(dsn = "./out_spatial/Merc.gpkg", layer = "MercBSM_Merged")

unique(m1$Hg_species_name)

## For all species
m1$Hg_length_mm <- m1$Hg_length_cm * 10
m1$Hg_value_log <- log(m1$Hg_value)
m1$Hg_length_mm_log <- log(m1$Hg_length_mm)
m1$Hg_weight_gram_log <- log(m1$Hg_weight_gram)
m1$Hg_waterbody_code <- as.character(m1$Hg_waterbody_code)

## All species
names(m1)
m2 <- m1[, c(
  "Hg_value_log", "Hg_length_mm_log", "Hg_waterbody_code",
  "Hg_species_name", "Hg_species_code", "Hg_fishcode_lake", "Hg_sample_year"
)]
m2 <- st_drop_geometry(m2)
m2 <- m2[complete.cases(m2), ]

## Outliers in X and Y
table(m1$Hg_value) # in 0.01 unit intervals
hist(m1$Hg_value, ylim = c(0, 15000), xlim = c(0, 7)) # strongly right skewed with majority of data less than 1 ug/g wet
dotchart(m2$Hg_value_log)
(check <- m1[m1$Hg_value_log < -3.5, ]) # Hg < 4 appear to still be measureable with no analytical flags

table(m1$Hg_length_mm) # in 1 unit intervals
hist(m1$Hg_length_mm, ylim = c(0, 10000), xlim = c(0, 2000)) # quite normally distrubted with a mean of about 500 mm = 50 cm
dotchart(m2$Hg_length_mm) # H
(check <- m1[m1$Hg_length_mm_log > 7, ])
## Choosing not to eliminate any outliers as they're certainly within the ranges of possibility
## Taking the log does squeeze some of this variation to avoid affecting the results too much
## Will evaluate after modelling to determine any outlier effects

## Expected X and Y relationship
plot(m2$Hg_value_log ~ m2$Hg_length_mm_log) # modelling this relationship
cor.test(m2$Hg_value_log, m2$Hg_length_mm_log) # 0.48

mod.lm <- lm(m2$Hg_value_log ~ m2$Hg_length_mm_log) # 23% of the variation in simple lm
abline(mod.lm)
hist(resid(mod.lm))
plot(fitted(mod.lm), resid(mod.lm))
## There are some "striping patterns but difficult to know if this is cause for concern until
## we have fit the final model.

## 1) Target species w/ standardized lengths ----

## From V. Danco's thesis
## For each lake, for each species, run a power regression
## Predict the value for a given species at a standard length
## Year might be different

## dataframe of target species and standardized lengths (i.e. sport fish)
species <- c(334, 131, 81, 80, 316)
std_lengths <- c(500, 650, 600, 300, 400)
std_mass <- 1000

(spec_standards <- data.frame(
  species_code = species,
  std_lengths = std_lengths,
  std_mass = std_mass
))

(spec_standards_med <- data.frame(
  species_code = species,
  std_lengths = sapply(species, function(xx){
    
    sel <- subset(m2, Hg_species_code == xx)
    sel_median <- median(sel$Hg_length_mm, na.rm = TRUE)  
    
  }),
  std_mass = sapply(species, function(xx){
    
    sel <- subset(m2, Hg_species_code == xx)
    sel_median <- median(sel$Hg_weight_gram, na.rm = TRUE)  
    
  })
))

## For condition, against the median condition of the fish across the
## whole dataset
spec_standards$condition <- sapply(species, function(x) {
  spec_subset <- subset(m2, Hg_species_code == x)
  k <- 100 * (spec_subset$Hg_weight_gram / spec_subset$Hg_length_cm^3)
  median(k, na.rm = T)
})

spec_standards_med$condition <- sapply(species, function(x) {
  spec_subset <- subset(m2, Hg_species_code == x)
  k <- 100 * (spec_subset$Hg_weight_gram / spec_subset$Hg_length_cm^3)
  median(k, na.rm = T)
})

## Only select those organsisms in the spec_standards table
m2 <- subset(m1, Hg_species_code %in% spec_standards$species_code)

## Generate 'condition'
m2$condition <- with(m2, 100 * (Hg_weight_gram / Hg_length_cm^3))

## 2) All species w/ dataset-based median length, mass, and condition ----

median_length <- sapply(unique(m1$Hg_species_name), function(x) {
  m1_sub <- subset(m1, Hg_species_name == x)
  res <- median(m1_sub$Hg_length_mm_log, na.rm = T)
  return(res)
}, USE.NAMES = TRUE)

median_mass <- sapply(unique(m1$Hg_species_name), function(x) {
  m1_sub <- subset(m1, Hg_species_name == x)
  res <- median(m1_sub$Hg_weight_gram_log, na.rm = TRUE)
  return(res)
}, USE.NAMES = TRUE)

median_condition <- sapply(unique(m1$Hg_species_name), function(x) {
  m1_sub <- subset(m1, Hg_species_name == x)
  m1_sub
  k <- 100 * (median(m1_sub$Hg_weight_gram, na.rm = T) / median(m1_sub$Hg_length_cm^3, na.rm = T))
  return(k)
}, USE.NAMES = TRUE)

## Get median size for all species and construct prediction dataframe
median_frame <- data.frame(
  Hg_species_name = unique(m1$Hg_species_name),
  Hg_length_mm_log = median_length,
  Hg_weight_gram_log = median_mass,
  Hg_condition = median_condition,
  stringsAsFactors = FALSE
)
median_frame

## Target dataframe
m1_nogeom <- st_drop_geometry(m1)

(u_wbsp <- unique(m1_nogeom[, c("Hg_waterbody_code", "Hg_sample_year", "Hg_species_name", "Hg_species_code")]))
u_wbsp <- merge(u_wbsp, median_frame, by = "Hg_species_name")
u_wbsp$Hg_value_log <- NA
u_wbsp$YRWB <- paste(u_wbsp$Hg_waterbody_code, u_wbsp$Hg_sample_year, sep = "_")
u_wbsp$Hg_species_name <- as.character(u_wbsp$Hg_species_name)

## Troubleshooter for lapply statements below
x <- 50559335 # waterbody
y <- 2015 # year
z <- 131 # species
# rm(x,y,z)

## *************************************
## SAMPLING EVENT-BASED REGRESSIONS ----
## *************************************

## Greater than 4 ----

## Standardize by length, >4 individuals needed and using Rob Mackereth's standardized values
m2_preds_length <- lapply(unique(m2$Hg_waterbody_code), function(x) {

  ## Subsets the waterbody of interest
  (sub_waterbody <- subset(m2, Hg_waterbody_code == x))

  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Hg_sample_year), function(y) {

    ## Subsets a year within a given waterbody
    (sub_year <- subset(sub_waterbody, Hg_sample_year == y))
    # print(y)

    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Hg_species_code), function(z) {

      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Hg_species_code == z))
      print(paste(x, y, z, sep = " - "))

      ## If number of individuals is greater than 4, compute regression
      if (nrow(sub_species) > 4) {

        ## log mercury concentrations and lengths (converted to mm)
        logm <- log(sub_species$Hg_value)
        logl <- log(sub_species$Hg_length_cm * 10)

        percentile <- ecdf(sub_species$Hg_length_cm * 10)

        ## regression relationship
        rel <- lm(logm ~ logl)

        # plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)

        ## Generating context and regression summary statistics
        (Hg_waterbody_code <- x) # waterbody
        (Hg_sample_year <- y) # year
        (Hg_species_code <- z) # species
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
        pred_modify <- spec_standards[spec_standards$species_code == z, "std_lengths"]

        (target_size_percentile <- percentile(pred_modify))

        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
          newdata = data.frame(logl = log(pred_modify)),
          interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))

        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Hg_waterbody_code, Hg_sample_year, Hg_species_code, n,
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
m2_preds_length <- m2_preds_length[!sapply(m2_preds_length, is.null)]
m2_preds_length <- do.call(plyr::rbind.fill, m2_preds_length)

## Standardize by mass, , >4 individuals needed and using Rob Mackereth's standardized values
m2_preds_mass <- lapply(unique(m2$Hg_waterbody_code), function(x) {

  ## Subsets the waterbody of interest
  sub_waterbody <- subset(m2, Hg_waterbody_code == x)

  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Hg_sample_year), function(y) {

    ## Subsets a year within a given waterbody
    sub_year <- subset(sub_waterbody, Hg_sample_year == y)
    # print(y)

    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Hg_species_code), function(z) {

      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Hg_species_code == z))
      print(paste(x, y, z, sep = " - "))

      ## If number of individuals is greater than 4, compute regression
      ## the sum evaluation statement is meant to make sure that all
      ## masses are not NAs which was a case in these data
      if (nrow(sub_species) > 4 & !sum(sub_species$Hg_weight_gram, na.rm = T) == 0) {

        ## log mercury concentrations and masss (converted to mm)
        (logm <- log(sub_species$Hg_value))
        (logl <- log(sub_species$Hg_weight_gram))

        percentile <- ecdf(sub_species$Hg_weight_gram)

        ## regression relationship
        rel <- lm(logm ~ logl)

        # plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)

        ## Generating context and regression summary statistics
        (Hg_waterbody_code <- x) # waterbody
        (Hg_sample_year <- y) # year
        (Hg_species_code <- z) # species
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

        ## Modify the predicted standardized mass based on species
        pred_modify <- spec_standards[spec_standards$species_code == z, "std_mass"]

        target_size_percentile <- percentile(pred_modify)

        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
          newdata = data.frame(logl = log(pred_modify)),
          interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))

        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Hg_waterbody_code, Hg_sample_year, Hg_species_code, n,
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
m2_preds_mass <- m2_preds_mass[!sapply(m2_preds_mass, is.null)]
m2_preds_mass <- do.call(plyr::rbind.fill, m2_preds_mass)

## Standardize by condition factor, , >4 individuals needed and using Rob Mackereth's standardized values
m2_preds_cond <- lapply(unique(m2$Hg_waterbody_code), function(x) {

  ## Subsets the waterbody of interest
  sub_waterbody <- subset(m2, Hg_waterbody_code == x)

  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Hg_sample_year), function(y) {

    ## Subsets a year within a given waterbody
    sub_year <- subset(sub_waterbody, Hg_sample_year == y)
    # print(y)

    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Hg_species_code), function(z) {

      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Hg_species_code == z))
      print(paste(x, y, z, sep = " - "))

      ## If number of individuals is greater than 4, compute regression
      ## the sum evaluation statement is meant to make sure that all
      ## condes are not NAs which was a case in these data
      if (nrow(sub_species) > 4 & !sum(sub_species$Hg_weight_gram, na.rm = T) == 0) {

        ## mercury concentrations as a function of condition
        (k <- 100 * (sub_species$Hg_weight_gram / sub_species$Hg_length_cm^3))

        percentile <- ecdf(k)

        ## regression relationship
        rel <- lm(sub_species$Hg_value ~ k)

        # summary(rel)
        # plot(sub_species$value ~ k, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)

        ## Generating context and regression summary statistics
        (Hg_waterbody_code <- x) # waterbody
        (Hg_sample_year <- y) # year
        (Hg_species_code <- z) # species
        (n <- length(k)) # number of individuals
        (int <- formatC(rel$coefficients[1], digits = 4, format = "f")) # estimated intercept
        (slp <- formatC(rel$coefficients[2], digits = 4, format = "f")) # estimated slope
        (int_confint <- paste(formatC(confint(rel)[1, ], digits = 4, format = "f"),
          collapse = " - "
        )) # estimated intercept confidence interval
        (slp_confint <- paste(formatC(confint(rel)[2, ], digits = 4, format = "f"),
          collapse = " - "
        )) # estimated slope confidence interval
        (adj_r2 <- formatC(summary(rel)$adj.r.squared, digits = 4, format = "f")) # adjusted R2

        ## Modify the predicted standardized cond based on species
        (pred_modify <- spec_standards[spec_standards$species_code == z, "condition"])
        (target_size_percentile <- percentile(pred_modify))

        ## Generated exponentiated prediction interval
        (pred <- (predict(rel,
          newdata = data.frame(k = (pred_modify)),
          interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))

        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Hg_waterbody_code, Hg_sample_year, Hg_species_code, n,
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
m2_preds_cond <- m2_preds_cond[!sapply(m2_preds_cond, is.null)]
m2_preds_cond <- do.call(plyr::rbind.fill, m2_preds_cond)

par(mfrow = c(2, 3))

## Boxplot of R2 by fish species, standardizing by length or mass
boxplot(as.numeric(m2_preds_length$adj_r2) ~ m2_preds_length$Hg_species_code,
  ylab = "Adjusted R2", xlab = "Species", main = "Length"
)
boxplot(as.numeric(m2_preds_mass$adj_r2) ~ m2_preds_mass$Hg_species_code,
  ylab = "Adjusted R2", xlab = "Species", main = "Mass"
)
boxplot(as.numeric(m2_preds_cond$adj_r2) ~ m2_preds_cond$Hg_species_code,
  ylab = "Adjusted R2", xlab = "Species", main = "Condition"
)

## Boxplot of percentile by fish species, standardizing by length
boxplot(as.numeric(m2_preds_length$target_size_percentile) ~ m2_preds_length$Hg_species_code,
  ylab = "Target size percentile of fish caught", xlab = "Species"
)
boxplot(as.numeric(m2_preds_mass$target_size_percentile) ~ m2_preds_mass$Hg_species_code,
  ylab = "Target size percentile of fish caught", xlab = "Species"
)
boxplot(as.numeric(m2_preds_cond$target_size_percentile) ~ m2_preds_cond$Hg_species_code,
  ylab = "Target size percentile of fish caught", xlab = "Species"
)

## Greater than 10 ---- 

## Standardize by length, >4 individuals needed
m2_preds_length_gte10 <- lapply(unique(m2$Hg_waterbody_code), function(x) {
  
  ## Subsets the waterbody of interest
  (sub_waterbody <- subset(m2, Hg_waterbody_code == x))
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Hg_sample_year), function(y) {
    
    ## Subsets a year within a given waterbody
    (sub_year <- subset(sub_waterbody, Hg_sample_year == y))
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Hg_species_code), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Hg_species_code == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      if (nrow(sub_species) > 9) {
        
        ## log mercury concentrations and lengths (converted to mm)
        logm <- log(sub_species$Hg_value)
        logl <- log(sub_species$Hg_length_cm * 10)
        
        percentile <- ecdf(sub_species$Hg_length_cm * 10)
        
        ## regression relationship
        rel <- lm(logm ~ logl)
        
        # plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)
        
        ## Generating context and regression summary statistics
        (Hg_waterbody_code <- x) # waterbody
        (Hg_sample_year <- y) # year
        (Hg_species_code <- z) # species
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
        pred_modify <- spec_standards_med[spec_standards_med$species_code == z, "std_lengths"]
        
        (target_size_percentile <- percentile(pred_modify))
        
        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
                             newdata = data.frame(logl = log(pred_modify)),
                             interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Hg_waterbody_code, Hg_sample_year, Hg_species_code, n,
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
m2_preds_length_gte10 <- m2_preds_length_gte10[!sapply(m2_preds_length_gte10, is.null)]
m2_preds_length_gte10 <- do.call(plyr::rbind.fill, m2_preds_length_gte10)

## Standardize by mass, >4 individuals needed
m2_preds_mass_gte10 <- lapply(unique(m2$Hg_waterbody_code), function(x) {
  
  ## Subsets the waterbody of interest
  sub_waterbody <- subset(m2, Hg_waterbody_code == x)
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Hg_sample_year), function(y) {
    
    ## Subsets a year within a given waterbody
    sub_year <- subset(sub_waterbody, Hg_sample_year == y)
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Hg_species_code), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Hg_species_code == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      ## the sum evaluation statement is meant to make sure that all
      ## masses are not NAs which was a case in these data
      if (nrow(sub_species) > 9 & !sum(sub_species$Hg_weight_gram, na.rm = T) == 0) {
        
        ## log mercury concentrations and masss (converted to mm)
        (logm <- log(sub_species$Hg_value))
        (logl <- log(sub_species$Hg_weight_gram))
        
        percentile <- ecdf(sub_species$Hg_weight_gram)
        
        ## regression relationship
        rel <- lm(logm ~ logl)
        
        # plot(logm ~ logl, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)
        
        ## Generating context and regression summary statistics
        (Hg_waterbody_code <- x) # waterbody
        (Hg_sample_year <- y) # year
        (Hg_species_code <- z) # species
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
        
        ## Modify the predicted standardized mass based on species
        pred_modify <- spec_standards_med[spec_standards_med$species_code == z, "std_mass"]
        
        target_size_percentile <- percentile(pred_modify)
        
        ## Generated exponentiated prediction interval
        (pred <- exp(predict(rel,
                             newdata = data.frame(logl = log(pred_modify)),
                             interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Hg_waterbody_code, Hg_sample_year, Hg_species_code, n,
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
m2_preds_mass_gte10 <- m2_preds_mass_gte10[!sapply(m2_preds_mass_gte10, is.null)]
m2_preds_mass_gte10 <- do.call(plyr::rbind.fill, m2_preds_mass_gte10)

## Standardize by condition factor, >4 needed
m2_preds_cond_gte10 <- lapply(unique(m2$Hg_waterbody_code), function(x) {
  
  ## Subsets the waterbody of interest
  sub_waterbody <- subset(m2, Hg_waterbody_code == x)
  
  ## Subsets a given year within the waterbody and then outputs a table
  ## for each species where a regression was fit
  sub_year_output <- lapply(unique(sub_waterbody$Hg_sample_year), function(y) {
    
    ## Subsets a year within a given waterbody
    sub_year <- subset(sub_waterbody, Hg_sample_year == y)
    # print(y)
    
    ## Subsets a species for a given year and waterbody and calculates
    ## a log-log regression when the number of individuals exceeds 4
    sub_species_output <- lapply(unique(sub_year$Hg_species_code), function(z) {
      
      ## Subsets the species-level dataframe
      (sub_species <- subset(sub_year, Hg_species_code == z))
      print(paste(x, y, z, sep = " - "))
      
      ## If number of individuals is greater than 4, compute regression
      ## the sum evaluation statement is meant to make sure that all
      ## condes are not NAs which was a case in these data
      if (nrow(sub_species) > 9 & !sum(sub_species$Hg_weight_gram, na.rm = T) == 0) {
        
        ## mercury concentrations as a function of condition
        (k <- 100 * (sub_species$Hg_weight_gram / sub_species$Hg_length_cm^3))
        
        percentile <- ecdf(k)
        
        ## regression relationship
        rel <- lm(sub_species$Hg_value ~ k)
        
        # summary(rel)
        # plot(sub_species$value ~ k, pch=16, col="black", main=paste(x,y,z, sep=" - "))
        # abline(rel)
        
        ## Generating context and regression summary statistics
        (Hg_waterbody_code <- x) # waterbody
        (Hg_sample_year <- y) # year
        (Hg_species_code <- z) # species
        (n <- length(k)) # number of individuals
        (int <- formatC(rel$coefficients[1], digits = 4, format = "f")) # estimated intercept
        (slp <- formatC(rel$coefficients[2], digits = 4, format = "f")) # estimated slope
        (int_confint <- paste(formatC(confint(rel)[1, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated intercept confidence interval
        (slp_confint <- paste(formatC(confint(rel)[2, ], digits = 4, format = "f"),
                              collapse = " - "
        )) # estimated slope confidence interval
        (adj_r2 <- formatC(summary(rel)$adj.r.squared, digits = 4, format = "f")) # adjusted R2
        
        ## Modify the predicted standardized cond based on species
        (pred_modify <- spec_standards_med[spec_standards_med$species_code == z, "condition"])
        (target_size_percentile <- percentile(pred_modify))
        
        ## Generated exponentiated prediction interval
        (pred <- (predict(rel,
                          newdata = data.frame(k = (pred_modify)),
                          interval = "prediction"
        )))
        (pred <- formatC(pred, digits = 4, format = "f"))
        
        ## Bring context and summary statistics together into dataframe
        (frame <- data.frame(Hg_waterbody_code, Hg_sample_year, Hg_species_code, n,
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
m2_preds_cond_gte10 <- m2_preds_cond_gte10[!sapply(m2_preds_cond_gte10, is.null)]
m2_preds_cond_gte10 <- do.call(plyr::rbind.fill, m2_preds_cond_gte10)

par(mfrow = c(2, 3))

## Boxplot of R2 by fish species, standardizing by length or mass
boxplot(as.numeric(m2_preds_length_gte10$adj_r2) ~ m2_preds_length_gte10$Hg_species_code,
        ylab = "Adjusted R2", xlab = "Species", main = "Length"
)
boxplot(as.numeric(m2_preds_mass_gte10$adj_r2) ~ m2_preds_mass_gte10$Hg_species_code,
        ylab = "Adjusted R2", xlab = "Species", main = "Mass"
)
boxplot(as.numeric(m2_preds_cond_gte10$adj_r2) ~ m2_preds_cond_gte10$Hg_species_code,
        ylab = "Adjusted R2", xlab = "Species", main = "Condition"
)

## Boxplot of percentile by fish species, standardizing by length
boxplot(as.numeric(m2_preds_length_gte10$target_size_percentile) ~ m2_preds_length_gte10$Hg_species_code,
        ylab = "Target size percentile of fish caught", xlab = "Species"
)
boxplot(as.numeric(m2_preds_mass_gte10$target_size_percentile) ~ m2_preds_mass_gte10$Hg_species_code,
        ylab = "Target size percentile of fish caught", xlab = "Species"
)
boxplot(as.numeric(m2_preds_cond_gte10$target_size_percentile) ~ m2_preds_cond_gte10$Hg_species_code,
        ylab = "Target size percentile of fish caught", xlab = "Species"
)

head(m2_preds_mass_gte10)

sd(m2_preds_mass_gte10$slp) # more variability in intercepts than slopes 
sd(m2_preds_mass_gte10$int) # 
hist(m2_preds_mass_gte10$n) # number of fish per lake
hist(as.numeric(m2_preds_mass_gte10$adj_r2))

spc_80 <- subset(m2_preds_mass_gte10, Hg_species_code == 80)
nrow(spc_80)
median(spc_80$n); min(spc_80$n); max(spc_80$n) 
sd(spc_80$int)
sd(spc_80$slp)
median(spc_80$adj_r2); min(spc_80$adj_r2); max(spc_80$adj_r2) 

spc_81 <- subset(m2_preds_mass_gte10, Hg_species_code == 81)
nrow(spc_81)
median(spc_81$n); min(spc_81$n); max(spc_81$n) 
sd(spc_81$int)
sd(spc_81$slp)
median(spc_81$adj_r2); min(spc_81$adj_r2); max(spc_81$adj_r2) 

spc_131 <- subset(m2_preds_mass_gte10, Hg_species_code == 131)
nrow(spc_131)
median(spc_131$n, na.rm = T); min(spc_131$n, na.rm = T); max(spc_131$n, na.rm = T) 
sd(spc_131$int, na.rm = T)
sd(spc_131$slp, na.rm = T)
median(as.numeric(spc_131$adj_r2), na.rm = T); min(as.numeric(spc_131$adj_r2), na.rm = T); max(as.numeric(spc_131$adj_r2), na.rm = T) 

spc_316 <- subset(m2_preds_mass_gte10, Hg_species_code == 316)
nrow(spc_316)
median(spc_316$n, na.rm = T); min(spc_316$n, na.rm = T); max(spc_316$n, na.rm = T) 
sd(spc_316$int, na.rm = T)
sd(spc_316$slp, na.rm = T)
median(as.numeric(spc_316$adj_r2), na.rm = T); min(as.numeric(spc_316$adj_r2), na.rm = T); max(as.numeric(spc_316$adj_r2), na.rm = T) 

spc_334 <- subset(m2_preds_mass_gte10, Hg_species_code == 334)
nrow(spc_334)
median(spc_334$n, na.rm = T); min(spc_334$n, na.rm = T); max(spc_334$n, na.rm = T) 
sd(spc_334$int, na.rm = T)
sd(spc_334$slp, na.rm = T)
median(as.numeric(spc_334$adj_r2), na.rm = T); min(as.numeric(spc_334$adj_r2), na.rm = T); max(as.numeric(spc_334$adj_r2), na.rm = T) 

spc_all <- m2_preds_mass_gte10
nrow(spc_all)
median(spc_all$n, na.rm = T); min(spc_all$n, na.rm = T); max(spc_all$n, na.rm = T) 
sd(spc_all$int, na.rm = T)
sd(spc_all$slp, na.rm = T)
median(as.numeric(spc_all$adj_r2), na.rm = T); min(as.numeric(spc_all$adj_r2), na.rm = T); max(as.numeric(spc_all$adj_r2), na.rm = T) 

length(unique(paste0(spc_all$Hg_waterbody_code, spc_all$Hg_sample_year)))





## ***************
## STAN MODEL ----
## ***************

names(m1)
m3 <- m1[, c(
  "Hg_value_log", "Hg_length_mm_log", "Hg_waterbody_code",
  "Hg_species_name", "Hg_species_code", "Hg_fishcode_lake", "Hg_sample_year"
)]
m3 <- st_drop_geometry(m3)
m3 <- m3[complete.cases(m3), ]

## Run full model, using data m2 with all species and waterbodies
## Takes about 16 hours to run
system.time(mod_full_length_stan <- stan_glmer(Hg_value_log ~ Hg_length_mm_log * Hg_species_name + (1 | Hg_waterbody_code / Hg_sample_year),
  data = m3, cores = 8, chains = 8
))
# saveRDS(mod_full_length_stan, "./out_workspaces/mod_full_length_stan.rds") #17 hours to run
mod_full_length_stan <- readRDS("./out_workspaces/mod_full_length_stan.rds")

## Model diagnostics

# shinystan
launch_shinystan(mod_full_length_stan)

# Marginal r-squared (no random effects)
rsq_marg <- bayes_R2(mod_full_length_stan, re.form = NA)
median(rsq_marg)

# Marginal r-squared (random effects) and so largely driven by fish and not their lake/sample year?
rsq_cond_lake <- bayes_R2(mod_full_length_stan, re.form = ~ (1 | Hg_waterbody_code / Hg_sample_year))
median(rsq_cond_lake)

## Residuals vs_ Fitted Values -- fitted values are medians
fits <- data.frame(
  Hg_value_log = m3$Hg_value_log,
  fitted_stan = fitted(mod_full_length_stan),
  resid_stan = residuals(mod_full_length_stan)
)

lw1_stan <- loess(resid_stan ~ fitted_stan, data = fits)
j <- order(fits$fitted_stan)
plot(resid_stan ~ fitted_stan, data = fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_stan$fitted[j] ~ fits$fitted_stan[j], lwd = 2)

## Histogram of residuals
hist(fits$resid_stan)

## Alternative comparison to fitted and residuals
posterior_full <- rstanarm::posterior_predict(mod_full_length_stan)
posterior_est_full <- apply(posterior_full, MARGIN = 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

fits$stan_posterior_q50 <- posterior_est_full[2, ]
fits$stan_posterior_q2p5 <- posterior_est_full[1, ]
fits$stan_posterior_q97p5 <- posterior_est_full[3, ]

## Parameter estimates from posterior distribution
mod_full_length_sims <- as.data.frame(as.matrix(mod_full_length_stan))
names(mod_full_length_sims)

a_quant <- apply(
  X = mod_full_length_sims,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2_5", "Q50", "Q97_5")

a_quant$param <- rownames(a_quant)

mod_full_length_stan_summary <- a_quant
mod_full_length_stan_summary$param

## Variance parameters
mod_full_length_stan

## Waterbody/Year = 0.15^2
## Waterbody = 0.45^2
## Residual = 0.36^2

## ************

## ************
## R-INLA  ----
## ************

## Code nested random effect
m3$YRWB <- paste(m3$Hg_waterbody_code, m3$Hg_sample_year, sep = "_")

mod_full_length_inla <- inla(Hg_value_log ~ Hg_length_mm_log * Hg_species_name + f(Hg_waterbody_code, model = "iid") +
  f(YRWB, model = "iid"),
data = m3,
control.predictor = list(
  compute = TRUE,
  quantiles = c(0.025, 0.5, 0.975)
),
control.compute = list(
  cpo = TRUE
)
)
mod_full_length_inla # 31 seconds

## Residuals vs_ Fitted Values -- fitted values are medians
fits$fitted_inla <- mod_full_length_inla$summary.fitted.values[, "0.5quant"]
fits$resid_inla <- fits$Hg_value_log - fits$fitted_inla

lw1_inla <- loess(resid_inla ~ fitted_inla, data = fits)
j <- order(fits$fitted_inla)
plot(resid_inla ~ fitted_inla, data = fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ fits$fitted_inla[j], lwd = 2)

## Alternative comparison to fitted and residuals
fits$inla_posterior_q50 <- mod_full_length_inla$summary.fitted.values[, "0.5quant"]
fits$inla_posterior_q2p5 <- mod_full_length_inla$summary.fitted.values[, "0.025quant"]
fits$inla_posterior_q97p5 <- mod_full_length_inla$summary.fitted.values[, "0.975quant"]

## R2
summary(lm(m3$Hg_value_log ~ fits$inla_posterior_q50)) # 0.80

## Histogram of residuals
hist(fits$resid_inla)

## Parameter estimates from posterior distribution
length_inla_fixed <- data.frame(
  ID = rownames(mod_full_length_inla$summary.fixed),
  mod_full_length_inla$summary.fixed, stringsAsFactors = FALSE
)
names(length_inla_fixed) <- names(mod_full_length_inla$summary.random$YRWB)

mod_full_length_inla_summary <- list(
  length_inla_fixed,
  mod_full_length_inla$summary.random$YRWB,
  mod_full_length_inla$summary.random$Hg_waterbody_code
)

(mod_full_length_inla_summary <- do.call("rbind", mod_full_length_inla_summary))

# High CPO implies good model fit
dotchart(mod_full_length_inla$cpo$cpo)
summary(mod_full_length_inla$cpo$cpo)

quantile(mod_full_length_inla$cpo$cpo, probs = c(0.0275, 0.5, 0.975))

## Another way is to look at posterior predictive check
pval <- rep(NA, nrow(m3))
for (i in 1:nrow(m3)) {
  print(i)
  pval[i] <- inla.pmarginal(
    q = fits$Hg_value_log[i],
    marginal = mod_full_length_inla$marginals.fitted.values[[i]]
  )
}

hist(pval) # apparently this indicates a poor model fit?

# Variance parameters

mod_full_length_inla$marginals.hyperpar$`Precision for the Gaussian observations`

tau_waterbody <- mod_full_length_inla$marginals.hyperpar$`Precision for Hg_waterbody_code`
tau_YRWB <- mod_full_length_inla$marginals.hyperpar$`Precision for YRWB`
tau_residual <- mod_full_length_inla$marginals.hyperpar$`Precision for the Gaussian observations`

MySqrt <- function(x) {
  1 / sqrt(x)
}

(sigma_waterbody <- inla.emarginal(MySqrt, tau_waterbody))^2
(sigma_YRWB <- inla.emarginal(MySqrt, tau_YRWB))^2
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))^2

summary

## ***************************
## RANDOM INTERCEPT/SLOPE ----
## *************************** 

## Taken from KRIG code

## Code nested random effect
m3$YRWB <- paste(m3$Hg_waterbody_code, m3$Hg_sample_year, sep = "_")
m3$WATERBODY_CODE <- factor(m3$Hg_waterbody_code)
m3$WATERBODY_CODE1 <- as.integer(m3$WATERBODY_CODE)
m3$WATERBODY_CODE2 <- m3$WATERBODY_CODE1 + max(m3$WATERBODY_CODE1)
m3_n_waterbody <- n_distinct(m3$WATERBODY_CODE)

## Null model ---- 

mod_rirs_null <- INLA::inla(Hg_value_log ~ 1 +
                              f(WATERBODY_CODE1, model = "iid") + 
                              f(YRWB, model = "iid"),
                            data = m3,
                            control.predictor = list(
                              compute = TRUE,
                              quantiles = c(0.025, 0.5, 0.975)
                            ),
                            control.compute = list(
                              cpo = TRUE,
                              config = TRUE
                            )
)

# Set up precision --> standard deviation formula 
MySqrt <- function(x) {
  1 / sqrt(x)
}

rirs_summary <- summary(mod_rirs_null)

fits_rirs <- data.frame(WATERBODY_CODE = m3$WATERBODY_CODE, 
                        Hg_value_log = m3$Hg_value_log)


## Alternative comparison to fitted and residuals
fits_rirs$inla_posterior_q50 <- mod_rirs_null$summary.fitted.values[, "0.5quant"]
fits_rirs$inla_posterior_q2p5 <- mod_rirs_null$summary.fitted.values[, "0.025quant"]
fits_rirs$inla_posterior_q97p5 <- mod_rirs_null$summary.fitted.values[, "0.975quant"]
fits_rirs$resid_inla <- fits_rirs$Hg_value_log - fits_rirs$inla_posterior_q50

## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = fits_rirs)
j <- order(fits_rirs$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = fits_rirs, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ fits_rirs$inla_posterior_q50[j], lwd = 2)
## Looks like trouble at the low concentration values (i.e., overpredicting?)

## R2
plot(fits_rirs$Hg_value_log ~ fits_rirs$inla_posterior_q50, pch = 16) 
fitted_test <- lm(fits_rirs$Hg_value_log ~ fits_rirs$inla_posterior_q50)
summary(fitted_test) 
abline(fitted_test)

# Variance parameters
mod_rirs_null$marginals.hyperpar$`Precision for the Gaussian observations`
mod_rirs_null$marginals.hyperpar$`Precision for WATERBODY_CODE1`
mod_rirs_null$marginals.hyperpar$`Precision for YRWB`

tau_WATERBODY_CODE1 <- mod_rirs_null$marginals.hyperpar$`Precision for WATERBODY_CODE1`
tau_YRWB <- mod_rirs_null$marginals.hyperpar$`Precision for YRWB`
tau_residual <- mod_rirs_null$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_YRWB <- inla.emarginal(MySqrt, tau_YRWB))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate posteriors 
(posterior_WATERBODY_CODE1 <- inla.tmarginal(MySqrt, tau_WATERBODY_CODE1))
(posterior_YRWB <- inla.tmarginal(MySqrt, tau_YRWB))
(posterior_residual <- inla.tmarginal(MySqrt, tau_residual))

par(mfrow = c(1,3))
plot(posterior_WATERBODY_CODE1, type = "l", main = "Lake")
plot(posterior_YRWB, type = "l", main = "YWRB")
plot(posterior_residual, type = "l", main = "residual")

## Generate variance components 
rirs_variance_components <- c(sigma_residual, sigma_YRWB, sigma_WATERBODY_CODE1)^2
rirs_variance_components

## Random slope/random intercept model ----

## Code nested random effect
m3$YRWB <- paste(m3$Hg_waterbody_code, m3$Hg_sample_year, sep = "_")
m3$WATERBODY_CODE <- factor(m3$Hg_waterbody_code)
m3$WATERBODY_CODE1 <- as.integer(m3$WATERBODY_CODE)
m3$WATERBODY_CODE2 <- m3$WATERBODY_CODE1 + max(m3$WATERBODY_CODE1)
m3_n_waterbody <- n_distinct(m3$WATERBODY_CODE)

mod_rirs <- INLA::inla(Hg_value_log ~ Hg_length_mm_log * Hg_species_name +
                         f(WATERBODY_CODE1, n = 2 * m3_n_waterbody, model = "iid2d") + 
                         f(WATERBODY_CODE2, Hg_length_mm_log, copy = "WATERBODY_CODE1") + 
                         f(YRWB, model = "iid"),
                       data = m3,
                       control.predictor = list(
                         compute = TRUE,
                         quantiles = c(0.025, 0.5, 0.975)
                       ),
                       control.compute = list(
                         cpo = TRUE,
                         config = TRUE
                       )
)

summary(mod_rirs)

length(unique(m3$YRWB))

this <- lapply(unique(m3$YRWB), function(xx){
  
  sel_yrwb <- subset(m3, YRWB == xx)
  sel_yrwb_spc_count <- table(sel_yrwb$Hg_species_code)  
})
median(unlist(this))
min(unlist(this))
max(unlist(this))

## Set up precision --> standard deviation formula 
MySqrt <- function(x) {
  1 / sqrt(x)
}

rirs_summary <- summary(mod_rirs)

fits_rirs <- data.frame(WATERBODY_CODE = m3$WATERBODY_CODE, 
                        Hg_value_log = m3$Hg_value_log)

## Alternative comparison to fitted and residuals
fits_rirs$inla_posterior_q50 <- mod_rirs$summary.fitted.values[, "0.5quant"]
fits_rirs$inla_posterior_q2p5 <- mod_rirs$summary.fitted.values[, "0.025quant"]
fits_rirs$inla_posterior_q97p5 <- mod_rirs$summary.fitted.values[, "0.975quant"]
fits_rirs$resid_inla <- fits_rirs$Hg_value_log - fits_rirs$inla_posterior_q50

## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = fits_rirs)
j <- order(fits_rirs$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = fits_rirs, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ fits_rirs$inla_posterior_q50[j], lwd = 2)
## Looks like trouble at the low concentration values (i.e., overpredicting?)

## R2
plot(fits_rirs$Hg_value_log ~ fits_rirs$inla_posterior_q50, pch = 16) 
fitted_test <- lm(fits_rirs$Hg_value_log ~ fits_rirs$inla_posterior_q50)
summary(fitted_test) 
abline(fitted_test)

## Histogram of residuals
hist(fits_rirs$resid_inla) # approximately normally distributed

## Histogram of random effect estimates 
hist(mod_rirs$summary.random$WATERBODY_CODE1$`0.5quant`) # approximately normally distributed 
hist(mod_rirs$summary.random$WATERBODY_CODE2$`0.5quant`) # approximately normally distributed 
hist(mod_rirs$summary.random$YRWB$`0.5quant`) # approximately normally distributed 

## Parameter estimates from posterior distribution
length_inla_fixed <- data.frame(
  ID = rownames(mod_rirs$summary.fixed),
  mod_rirs$summary.fixed, stringsAsFactors = FALSE
)

names(length_inla_fixed) <- c("ID", names(mod_rirs$summary.fixed))
length_inla_fixed

mod_rirs_summary <- list(
  length_inla_fixed,
  mod_rirs$summary.random$WATERBODY_CODE1,
  mod_rirs$summary.random$WATERBODY_CODE2
)

(mod_rirs_summary <- do.call("rbind", mod_rirs_summary))

# High CPO implies good model fit
dotchart(mod_rirs$cpo$cpo)
summary(mod_rirs$cpo$cpo)

quantile(mod_rirs$cpo$cpo, probs = c(0.0275, 0.5, 0.975))

## Another way is to look at posterior predictive check
pval <- rep(NA, nrow(m3))
for (i in 1:nrow(m3)) {
  print(i)
  pval[i] <- inla.pmarginal(
    q = fits_rirs$Hg_value_log[i],
    marginal = mod_rirs$marginals.fitted.values[[i]]
  )
}

hist(pval) # apparently this indicates a poor model fit?

# Variance parameters
mod_rirs$marginals.hyperpar$`Precision for the Gaussian observations`
mod_rirs$marginals.hyperpar$`Precision for YRWB`
mod_rirs$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
mod_rirs$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
mod_rirs$marginals.hyperpar$`Rho1:2 for WATERBODY_CODE1`

tau_WATERBODY_CODE1 <- mod_rirs$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- mod_rirs$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_YRWB <- mod_rirs$marginals.hyperpar$`Precision for YRWB`
tau_residual <- mod_rirs$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_YRWB <- inla.emarginal(MySqrt, tau_YRWB))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate posteriors 
(posterior_WATERBODY_CODE1 <- inla.tmarginal(MySqrt, tau_WATERBODY_CODE1))
(posterior_WATERBODY_CODE1_RSLOPE <- inla.tmarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(posterior_YRWB <- inla.tmarginal(MySqrt, tau_YRWB))
(posterior_residual <- inla.tmarginal(MySqrt, tau_residual))

par(mfrow = c(1,4))
plot(posterior_WATERBODY_CODE1, type = "l", main = "Lake")
plot(posterior_WATERBODY_CODE1_RSLOPE, type = "l", main = "Lake Slope")
plot(posterior_YRWB, type = "l", main = "Lake Slope")
plot(posterior_residual, type = "l", main = "residual")

## Generate variance components 
rirs_variance_components <- c(sigma_residual, sigma_YRWB, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE)^2
rirs_variance_components

rirs_variance_components[1] / sum(rirs_variance_components)
rirs_variance_components[2] / sum(rirs_variance_components)
rirs_variance_components[3] / sum(rirs_variance_components) 

#rirs_prediction_values <- c(454) ## TBD

rirs_predict <- aggregate(Hg_value_log ~ Hg_species_name + Hg_species_code + WATERBODY_CODE1 + WATERBODY_CODE2 + YRWB, data = m3, FUN = mean )
rirs_predict <- merge(rirs_predict, median_frame, by = "Hg_species_name")
rirs_predict$Hg_value_log_orig <- rirs_predict$Hg_value_log
rirs_predict$Hg_value_log <- NA

m4 <- bind_rows(m3, rirs_predict)

mod_rirs_predict <- INLA::inla(Hg_value_log ~ Hg_length_mm_log * Hg_species_name +
                                 f(WATERBODY_CODE1, n = 2 * m3_n_waterbody, model = "iid2d") + 
                                 f(WATERBODY_CODE2, Hg_length_mm_log, copy = "WATERBODY_CODE1") + 
                                 f(YRWB, model = "iid"),
                               data = m4,
                               control.predictor = list(
                                 compute = TRUE,
                                 quantiles = c(0.025, 0.5, 0.975)
                               ),
                               control.compute = list(
                                 cpo = TRUE,
                                 config = TRUE
                               )
)

mod_rirs_predict_posteriors <- mod_rirs_predict$summary.fitted.values[
  (nrow(m3) + 1):nrow(m4),
  c("0.025quant", "0.5quant", "0.975quant")
  ]
mod_rirs_predict_posteriors <- exp(mod_rirs_predict_posteriors)

(rirs_predictions <- cbind(rirs_predict, mod_rirs_predict_posteriors))

# 
# 
# 
# predictions <- merge(predictions, m3[,c("Hg_waterbody_code", "WATERBODY_CODE1")], 
#                      by = "WATERBODY_CODE1")
# predictions <- predictions[!duplicated(predictions),]
# 
# 
# plot(exp(predictions$Hg_value_log_orig) ~ predictions$`0.5quant`, 
#      pch = 21, col = "lightgrey", bg = adjustcolor("black", alpha.f = 0.5), 
#      xlab = "[Hg], INLA 0.5 quantile", ylab = "[Hg], observed")
# 
# names(predictions)
# 
# names(predictions)[10:12] <- c("Inla_0_025", "Inla_0_5", "Inla_0_975")
# 
# write.csv(predictions,
#           "./out_tables/Merc_FishHgStandardize_Predictions_RIRSModel.csv",
#           row.names = FALSE
# )

## Random slope/random intercept model - sampling event comparison ----

unique(m2$Hg_species_code)

m2_comp <- m2

## Code nested random effect
m2_comp$YRWB <- paste(m2_comp$Hg_waterbody_code, m2_comp$Hg_sample_year, sep = "_")
m2_comp$WATERBODY_CODE <- factor(m2_comp$Hg_waterbody_code)
m2_comp$WATERBODY_CODE1 <- as.integer(m2_comp$WATERBODY_CODE)
m2_comp$WATERBODY_CODE2 <- m2_comp$WATERBODY_CODE1 + max(m2_comp$WATERBODY_CODE1)
m2_comp_n_waterbody <- n_distinct(m2_comp$WATERBODY_CODE)

mod_rirs_comp <- INLA::inla(Hg_value_log ~ Hg_length_mm_log * Hg_species_name +
                              #f(WATERBODY_CODE1, model = "iid") +
                              f(WATERBODY_CODE1, n = 2 * m2_comp_n_waterbody, model = "iid2d") + 
                              f(WATERBODY_CODE2, Hg_length_mm_log, copy = "WATERBODY_CODE1") + 
                              f(YRWB, model = "iid"),
                            data = m2_comp,
                            control.predictor = list(
                              compute = TRUE,
                              quantiles = c(0.025, 0.5, 0.975)
                            ),
                            control.compute = list(
                              cpo = TRUE,
                              config = TRUE
                            )
)

summary(mod_rirs_comp)

length(unique(m2_comp$YRWB))

## For each year waterbody, for each species, get count

xx <- m2_comp$YRWB[1]

this <- lapply(unique(m2_comp$YRWB), function(xx){
  
  sel_yrwb <- subset(m2_comp, YRWB == xx)
  sel_yrwb_spc_count <- table(sel_yrwb$Hg_species_code)  
})
median(unlist(this))
min(unlist(this))
max(unlist(this))


## Set up precision --> standard deviation formula 
MySqrt <- function(x) {
  1 / sqrt(x)
}

rirs_comp_summary <- summary(mod_rirs_comp)

fits_rirs_comp <- data.frame(WATERBODY_CODE = m2_comp$WATERBODY_CODE, 
                             Hg_value_log = m2_comp$Hg_value_log, 
                             Hg_species_code = m2_comp$Hg_species_code, 
                             Hg_species_name = m2_comp$Hg_species_name)

## Alternative comparison to fitted and residuals
fits_rirs_comp$inla_posterior_q50 <- mod_rirs_comp$summary.fitted.values[, "0.5quant"]
fits_rirs_comp$inla_posterior_q2p5 <- mod_rirs_comp$summary.fitted.values[, "0.025quant"]
fits_rirs_comp$inla_posterior_q97p5 <- mod_rirs_comp$summary.fitted.values[, "0.975quant"]
fits_rirs_comp$resid_inla <- fits_rirs_comp$Hg_value_log - fits_rirs_comp$inla_posterior_q50

## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = fits_rirs_comp)
j <- order(fits_rirs_comp$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = fits_rirs_comp, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ fits_rirs_comp$inla_posterior_q50[j], lwd = 2)
## Looks like trouble at the low concentration values (i.e., overpredicting?)

## R2
plot(fits_rirs_comp$Hg_value_log ~ fits_rirs_comp$inla_posterior_q50, pch = 16) 
fitted_test <- lm(fits_rirs_comp$Hg_value_log ~ fits_rirs_comp$inla_posterior_q50)
summary(fitted_test) 
abline(fitted_test)

## Histogram of residuals
hist(fits_rirs_comp$resid_inla) # approximately normally distributed

## Histogram of random effect estimates 
hist(mod_rirs_comp$summary.random$WATERBODY_CODE1$`0.5quant`) # approximately normally distributed 
hist(mod_rirs_comp$summary.random$WATERBODY_CODE2$`0.5quant`) # approximately normally distributed 
hist(mod_rirs_comp$summary.random$YRWB$`0.5quant`) # approximately normally distributed 

## Parameter estimates from posterior distribution
length_inla_fixed <- data.frame(
  ID = rownames(mod_rirs_comp$summary.fixed),
  mod_rirs_comp$summary.fixed, stringsAsFactors = FALSE
)

names(length_inla_fixed) <- c("ID", names(mod_rirs_comp$summary.fixed))
length_inla_fixed

mod_rirs_comp_summary <- list(
  length_inla_fixed,
  mod_rirs_comp$summary.random$WATERBODY_CODE1,
  mod_rirs_comp$summary.random$WATERBODY_CODE2
)

(mod_rirs_comp_summary <- do.call("rbind", mod_rirs_comp_summary))

# High CPO implies good model fit
dotchart(mod_rirs_comp$cpo$cpo)
summary(mod_rirs_comp$cpo$cpo)

quantile(mod_rirs_comp$cpo$cpo, probs = c(0.0275, 0.5, 0.975))

## Another way is to look at posterior predictive check
pval <- rep(NA, nrow(m2_comp))
for (i in 1:nrow(m2_comp)) {
  print(i)
  pval[i] <- inla.pmarginal(
    q = fits_rirs_comp$Hg_value_log[i],
    marginal = mod_rirs_comp$marginals.fitted.values[[i]]
  )
}

hist(pval) # apparently this indicates a poor model fit?

# Variance parameters
mod_rirs_comp$marginals.hyperpar$`Precision for the Gaussian observations`
mod_rirs_comp$marginals.hyperpar$`Precision for YRWB`
mod_rirs_comp$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
mod_rirs_comp$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
mod_rirs_comp$marginals.hyperpar$`Rho1:2 for WATERBODY_CODE1`

tau_WATERBODY_CODE1 <- mod_rirs_comp$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- mod_rirs_comp$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_YRWB <- mod_rirs_comp$marginals.hyperpar$`Precision for YRWB`
tau_residual <- mod_rirs_comp$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1))
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_YRWB <- inla.emarginal(MySqrt, tau_YRWB))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate posteriors 
(posterior_WATERBODY_CODE1 <- inla.tmarginal(MySqrt, tau_WATERBODY_CODE1))
(posterior_WATERBODY_CODE1_RSLOPE <- inla.tmarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(posterior_YRWB <- inla.tmarginal(MySqrt, tau_YRWB))
(posterior_residual <- inla.tmarginal(MySqrt, tau_residual))

par(mfrow = c(1,4))
plot(posterior_WATERBODY_CODE1, type = "l", main = "Lake")
plot(posterior_WATERBODY_CODE1_RSLOPE, type = "l", main = "Lake Slope")
plot(posterior_YRWB, type = "l", main = "Lake Slope")
plot(posterior_residual, type = "l", main = "residual")

## Generate variance components 
rirs_comp_variance_components <- c(sigma_residual, sigma_YRWB, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE)^2
rirs_comp_variance_components

rirs_comp_variance_components[1] / sum(rirs_comp_variance_components)
rirs_comp_variance_components[2] / sum(rirs_comp_variance_components)
rirs_comp_variance_components[3] / sum(rirs_comp_variance_components) 





## *******************
## R-INLA VS STAN ----
## *******************

##
head(u_wbsp)

names(mod_full_length_stan$data) %in% names(u_wbsp)

stan_pred <- u_wbsp[, names(mod_full_length_stan$data)[-6]]

## Predictions, rstanarm
posterior_wbsp <- rstanarm::posterior_predict(mod_full_length_stan,
  newdata = stan_pred[, -1]
)
posterior_est_wbsp <- apply(posterior_wbsp,
  MARGIN = 2,
  function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }
)
posterior_est_wbsp_stan <- (t(posterior_est_wbsp))

## Predictions, f

inla_pred <- u_wbsp[, names(m3)[-6]]
inla_pred <- rbind(m3[, -6], inla_pred) # Hg_value_log has NA

mod_full_length_inla_preds <- inla(Hg_value_log ~ Hg_length_mm_log *
  Hg_species_name + f(Hg_waterbody_code, model = "iid") +
  f(YRWB, model = "iid"),
data = inla_pred,
control.predictor = list(
  compute = TRUE,
  quantiles = c(0.025, 0.5, 0.975)
),
control.compute = list(
  cpo = TRUE
), verbose = FALSE
)
mod_full_length_inla_preds # 31 seconds

nrow(mod_full_length_inla_preds$summary.fitted.values) - nrow(m2)

posterior_est_wbsp_inla <- mod_full_length_inla_preds$summary.fitted.values[
  (nrow(m3) + 1):nrow(inla_pred),
  c("0.025quant", "0.5quant", "0.975quant")
]

head(posterior_est_wbsp_inla)
head(posterior_est_wbsp_stan)

plot(posterior_est_wbsp_inla$`0.025quant`)
plot(posterior_est_wbsp_stan[, 2])

predictions <- cbind(u_wbsp, posterior_est_wbsp_inla, posterior_est_wbsp_stan)

## INLA
boxplot(predictions$`0.5quant` ~ predictions$Hg_sample_year)
boxplot(predictions$`0.5quant` ~ predictions$Hg_species_name, las = 2)

## Stan
boxplot(predictions$`50%` ~ predictions$Hg_sample_year)
boxplot(predictions$`50%` ~ predictions$Hg_species_name)

names(predictions)

c("Inla_0_025", "Inla_0_5", "Inla_0_975")

predictions_out <- predictions[, c(1:7, 10:12)]
names(predictions_out)[8:10] <- c("Inla_0_025", "Inla_0_5", "Inla_0_975")

write.csv(predictions_out,
  "./out_tables/Merc_FishHgStandardize_Predictions.csv",
  row.names = FALSE
)

## *************************************
## SAMPLING EVENT VS R-INLA VS STAN ----
## *************************************

## Prepare a comparison dataframe based on m2_preds_length
head(m2_preds_length)

com <- m2_preds_length[, c(
  "Hg_waterbody_code", "Hg_sample_year",
  "Hg_species_code"
)]

com <- merge(com, m2[,c("Hg_species_code", "Hg_species_name")], by = "Hg_species_code")
com

com <- merge(com, spec_standards_med, by.x = "Hg_species_code", by.y = "species_code")
com <- merge(com, unique(predictions[, c("Hg_species_code", "Hg_species_name")]))
com$Hg_length_mm_log <- log(com$std_lengths)

## Stan predictions
names(stan_pred) %in% names(com)
com_stan <- com[, names(stan_pred)[-1]] # remove Hg_value_log

posterior_com <- rstanarm::posterior_predict(mod_full_length_stan,
  newdata = com_stan
)
posterior_est_com <- apply(posterior_com,
  MARGIN = 2,
  function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }
)
posterior_est_com_stan <- (t(posterior_est_com))

## R-INLA predictions
names(inla_pred) %in% names(com_stan)

com_inla <- com
com_inla$YRWB <- paste(com_inla$Hg_waterbody_code, com_inla$Hg_sample_year, sep = "_")
com_inla$Hg_value_log <- NA
com_inla <- com_inla[, names(m3)[-6]] # remove fishcode
com_inla <- rbind(m3[, -6], com_inla) # Hg_value_log has NA

mod_full_length_inla_com <- inla(Hg_value_log ~ Hg_length_mm_log *
  Hg_species_name + f(Hg_waterbody_code, model = "iid") +
  f(YRWB, model = "iid"),
data = com_inla,
control.predictor = list(
  compute = TRUE,
  quantiles = c(0.025, 0.5, 0.975)
),
control.compute = list(
  cpo = TRUE
), verbose = FALSE
)
mod_full_length_inla_com # 31 seconds

posterior_est_com_inla <- mod_full_length_inla_com$summary.fitted.values[
  (nrow(m3) + 1):nrow(com_inla),
  c("0.025quant", "0.5quant", "0.975quant")
]

## Bringing together com, stan, and inla predictions
predictions_com <- cbind(com, posterior_est_com_inla, posterior_est_com_stan)
predictions_com$WBYRSPC <- paste(predictions_com$Hg_waterbody_code,
  predictions_com$Hg_sample_year,
  predictions_com$Hg_species_code,
  sep = "-"
)

m2_preds_length$WBYRSPC <- paste(m2_preds_length$Hg_waterbody_code,
  m2_preds_length$Hg_sample_year,
  m2_preds_length$Hg_species_code,
  sep = "-"
)

identical(sort(predictions_com$WBYRSPC), sort(m2_preds_length$WBYRSPC))

predictions_com_sub <- predictions_com[, c(
  "WBYRSPC", "0.025quant", "0.5quant",
  "0.975quant", "2.5%", "50%", "97.5%"
)]


names(predictions_com_sub)[2:7] <- c(
  "Inla_0_025", "Inla_0_5", "Inla_0_975",
  "Stan_0_025", "Stan_0_5", "Stan_0_975"
)

## Exponentiate Stan and INLA predictions
predictions_com_sub[2:7] <- lapply(predictions_com_sub[2:7], function(x) {
  exp(x)
})

head(predictions_com_sub)

final_com <- merge(m2_preds_length, predictions_com_sub)
head(final_com)

names(final_com)[6:14] <- paste0("SER_", names(final_com)[6:14])

plot(final_com$SER_fit ~ final_com$Inla_0_5)


## Comparison - Random slope random intercept ----

## Prepare a comparison dataframe based on m2_preds_length
head(m2_preds_length_gte10)

m2_preds_length_gte10$WYS <- paste(m2_preds_length_gte10$Hg_waterbody_code, 
                                   m2_preds_length_gte10$Hg_sample_year, 
                                   m2_preds_length_gte10$Hg_species_code, sep = "_")

head(rirs_predictions)

rirs_predictions$WYS <- paste(rirs_predictions$YRWB, 
                              rirs_predictions$Hg_species_code, sep = "_")

pred_compare <- merge(m2_preds_length_gte10, rirs_predictions[,c("WYS", "0.5quant")])

par(mfrow = c(2,3))

lapply(1:length(unique(pred_compare$Hg_species_code)), function(xx){
  
  sel <- subset(pred_compare, Hg_species_code == unique(pred_compare$Hg_species_code)[xx])
  
  plot(pred_compare$fit ~ pred_compare$`0.5quant`, pch = 16, col = "black", 
       main = unique(pred_compare$Hg_species_code)[xx], 
       xlab = "Bayesian", ylab = "Lake-level")
  pred_compare_lm <- lm(pred_compare$fit ~ pred_compare$`0.5quant`)
  abline(pred_compare_lm, lty = 2)
  
  points(sel$fit ~ sel$`0.5quant`, pch = 16, col = xx+1)
  sel_lm <- lm(sel$fit ~ sel$`0.5quant`)
  abline(sel_lm, lty = 2, col = xx+1)

})

saveRDS(pred_compare, "./out_workspaces/pred_compare.rds")

getwd()
plot(this$fit, this$`0.5quant`)

cor.test(as.numeric(this$fit), this$`0.5quant`) #0.92



this2 <- lm(this$fi)

head(this)


com_inla[,-c(4:6)]

names(inla_pred) %in% names(com_stan)

com_inla <- com
com_inla$YRWB <- paste(com_inla$Hg_waterbody_code, com_inla$Hg_sample_year, sep = "_")
com_inla$Hg_value_log <- NA

com_inla <- left_join(com_inla, m3[,c("YRWB", "WATERBODY_CODE", "WATERBODY_CODE1", "WATERBODY_CODE2")], 
                  by = "YRWB")

com_inla <- com_inla[!duplicated(com_inla),]

com_inla_check <- rbind(m3[,names(com_inla)[-c(4:6)]], com_inla[,-c(4:6)])
com_inla_check


head(com_inla_check)

com_inla <- com_inla[, names(m3)[-6]] # remove fishcode
com_inla <- rbind(m3[, -6], com_inla) # Hg_value_log has NA


## ************************************
## POOR PREDICTOR MODEL COMPARISON ----
## ************************************

head(pred_compare)

pred_compare$slp_confint

slp_confint_splt <- strsplit(pred_compare$slp_confint, split = " ")

(confint_lower <- as.numeric(sapply(slp_confint_splt, function(x){x[1]})))
(confint_upper <- as.numeric(sapply(slp_confint_splt, function(x){x[3]})))

zero_slp <- which((confint_upper - confint_lower)>confint_upper)
zero_slp <- pred_compare[zero_slp,]

zero_slp$Hg_value_log_pop_mean <-   sapply(1:nrow(zero_slp), function(x){
  
  sel <- subset(m2, Hg_waterbody_code == zero_slp$Hg_waterbody_code[x] & 
             Hg_sample_year == zero_slp$Hg_sample_year[x] & 
             Hg_species_code == zero_slp$Hg_species_code[x])
  
  (mean(sel$Hg_value))  
})

par(mfrow=c(1,3))
cor.test(as.numeric(zero_slp$`0.5quant`), as.numeric(zero_slp$fit)) #0.93
plot(zero_slp$fit, zero_slp$`0.5quant`, ylim = c(0,3), xlim = c(0,3), 
     xlab = "Sampling event regressions", 
     ylab = "Bayesian model"); abline(0,1, lty = 2); text(2.5, 0.5, "0.93")

cor.test(as.numeric(zero_slp$`0.5quant`), as.numeric(zero_slp$Hg_value_log_pop_mean)) #0.86
plot(zero_slp$Hg_value_log_pop_mean~zero_slp$`0.5quant`, ylim = c(0,3), xlim = c(0,3), 
     xlab = "Bayesian model", 
     ylab = "Population means"); abline(0,1, lty = 2); text(2.5, 0.5, "0.86")

cor.test(as.numeric(zero_slp$fit), as.numeric(zero_slp$Hg_value_log_pop_mean)) #0.95
plot(zero_slp$Hg_value_log_pop_mean~zero_slp$fit, ylim = c(0,3), xlim = c(0,3), 
     xlab = "Sampling event regressions", 
     ylab = "Population means"); abline(0,1, lty = 2); text(2.5, 0.5, "0.95")


nzero_slp <- which(!((confint_upper - confint_lower)>confint_upper))
nzero_slp <- pred_compare[nzero_slp,]

nzero_slp$Hg_value_log_pop_mean <-   sapply(1:nrow(nzero_slp), function(x){
  
  sel <- subset(m2, Hg_waterbody_code == nzero_slp$Hg_waterbody_code[x] & 
                  Hg_sample_year == nzero_slp$Hg_sample_year[x] & 
                  Hg_species_code == nzero_slp$Hg_species_code[x])
  
  (mean(sel$Hg_value))  
})

cor.test(as.numeric(nzero_slp$`0.5quant`), as.numeric(nzero_slp$fit)) #0.93
plot(nzero_slp$fit, nzero_slp$`0.5quant`, ylim = c(0,3), xlim = c(0,3), 
     xlab = "Sampling event regressions", 
     ylab = "Bayesian model"); abline(0,1, lty = 2); text(2.5, 0.5, "0.93")

cor.test(as.numeric(nzero_slp$`0.5quant`), as.numeric(nzero_slp$Hg_value_log_pop_mean)) #0.83
plot(nzero_slp$Hg_value_log_pop_mean~nzero_slp$`0.5quant`, ylim = c(0,3), xlim = c(0,3), 
     xlab = "Bayesian model", 
     ylab = "Population means"); abline(0,1, lty = 2); text(2.5, 0.5, "0.83")

cor.test(as.numeric(nzero_slp$fit), as.numeric(nzero_slp$Hg_value_log_pop_mean)) #0.80
plot(nzero_slp$Hg_value_log_pop_mean~nzero_slp$fit, ylim = c(0,3), xlim = c(0,3), 
     xlab = "Sampling event regressions", 
     ylab = "Population means"); abline(0,1, lty = 2); text(2.5, 0.5, "0.80")











plot(zero_slp$fit, zero_slp$Hg_value_log_pop_mean, ylim = c(0,3), xlim = c(0,3)); abline(0,1, lty = 2)
plot(zero_slp$`0.5quant`, zero_slp$Hg_value_log_pop_mean); abline(0,1, lty = 2)

cor.test(as.numeric(zero_slp$fit), zero_slp$Hg_value_log_pop_mean)
cor.test(as.numeric(zero_slp$`0.5quant`), zero_slp$Hg_value_log_pop_mean)


zero_slp$Hg_value_log_pop_mean


hist(as.numeric(pred_compare[zero_slp,]$adj_r2))

## ******************************
## SPATIAL DISTANCE ANALYSIS ----
## ******************************

head(predictions_com)

## Subset by species
predictions_spatial <- merge(predictions, m1[, c("Hg_waterbody_code")], by = "Hg_waterbody_code")
predictions_spatial <- st_as_sf(predictions_spatial)

head(predictions_spatial)
names(predictions_spatial)

names(predictions_spatial)[c(10:12, 13:15)] <- c(
  "Inla_0_025", "Inla_0_5", "Inla_0_975",
  "Stan_0_025", "Stan_0_5", "Stan_0_975"
)

names(predictions_spatial)

predictions_spatial <- predictions_spatial[!duplicated(data_frame(predictions_spatial)), ]
predictions_spatial$Inla_0_5_exp <- exp(predictions_spatial$Inla_0_5)
jenks_breaks <- classIntervals(predictions_spatial$Inla_0_5_exp, style = "jenks", n = 5)

predictions_spatial$jenks_col <- viridis(5)[cut(predictions_spatial$Inla_0_5_exp,
  breaks = jenks_breaks$brks
)]

predictions_spatial_np <- subset(predictions_spatial, Hg_species_name == "Northern Pike") # 639
predictions_spatial_wa <- subset(predictions_spatial, Hg_species_name == "Walleye") # 618
predictions_spatial_sb <- subset(predictions_spatial, Hg_species_name == "Smallmouth Bass") # 327

## Northern Pike

De <- as.data.frame(st_coordinates(predictions_spatial_np))
De <- dist(De, method = "euclidean", diag = T, upper = T)
mems <- (adespatial::dbmem(De, MEM.autocor = "positive"))

predictions_spatial_np <- cbind(data.frame(predictions_spatial_np), mems)
predictions_spatial_np <- st_as_sf(predictions_spatial_np)

nrow(predictions_spatial_np)

## Quantile random forest regression

names(predictions_spatial_np)
p_voi <- seq(from = 18, to = 153)

rf_frame <- predictions_spatial_np
rf_frame <- st_drop_geometry(rf_frame)

head(rf_frame)

rf <- quantregForest(
  y = rf_frame[, "Inla_0_5_exp"],
  x = rf_frame[, p_voi],
  mtry = 6, importance = TRUE, keep.inbag = T, ntree = 2000,
  keep.forest = TRUE, nodesize = 5
)

rf

rf2 <- rf
class(rf2) <- "randomForest"

rf2 # 44% of the variation in Hg can be related to spatial distance properties

par(mfrow = c(1, 2))
plot(rf2$importance[, 1], type = "l", pch = 16, col = "blue", lwd = 2)
varImpPlot(rf2, type = 1)

partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM115")
partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM26")
partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM2")

## Walleye
De <- as.data.frame(st_coordinates(predictions_spatial_wa))
De <- dist(De, method = "euclidean", diag = T, upper = T)
mems <- (adespatial::dbmem(De, MEM.autocor = "positive"))

predictions_spatial_wa <- cbind(data.frame(predictions_spatial_wa), mems)
predictions_spatial_wa <- st_as_sf(predictions_spatial_wa)

nrow(predictions_spatial_wa)

## Quantile random forest regression

names(predictions_spatial_wa)
p_voi <- seq(from = 18, to = 149)

rf_frame <- predictions_spatial_wa
rf_frame <- st_drop_geometry(rf_frame)

head(rf_frame)

rf <- quantregForest(
  y = rf_frame[, "Inla_0_5_exp"],
  x = rf_frame[, p_voi],
  mtry = 6, importance = TRUE, keep.inbag = T, ntree = 2000,
  keep.forest = TRUE, nodesize = 5
)

rf

rf2 <- rf
class(rf2) <- "randomForest"

rf2 # 56% of the variation in Hg can be related to spatial distance properties

par(mfrow = c(1, 2))
plot(rf2$importance[, 1], type = "l", pch = 16, col = "blue", lwd = 2)
varImpPlot(rf2, type = 1)

partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM115")
partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM26")
partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM2")

## Walleye
De <- as.data.frame(st_coordinates(predictions_spatial_sb))
De <- dist(De, method = "euclidean", diag = T, upper = T)
mems <- (adespatial::dbmem(De, MEM.autocor = "positive"))

predictions_spatial_sb <- cbind(data.frame(predictions_spatial_sb), mems)
predictions_spatial_sb <- st_as_sf(predictions_spatial_sb)

## Quantile random forest regression

names(predictions_spatial_sb)
p_voi <- seq(from = 18, to = 70)

rf_frame <- predictions_spatial_sb
rf_frame <- st_drop_geometry(rf_frame)

head(rf_frame)

rf <- quantregForest(
  y = rf_frame[, "Inla_0_5_exp"],
  x = rf_frame[, p_voi],
  mtry = 6, importance = TRUE, keep.inbag = T, ntree = 2000,
  keep.forest = TRUE, nodesize = 5
)

rf

rf2 <- rf
class(rf2) <- "randomForest"

rf2 # 35% of the variation in Hg can be related to spatial distance properties

par(mfrow = c(1, 2))
plot(rf2$importance[, 1], type = "l", pch = 16, col = "blue", lwd = 2)
varImpPlot(rf2, type = 1)

partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM53")
partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM26")
partialPlot(rf2, pred.data = rf_frame[, p_voi], x.var = "MEM2")

## **********
## PLOTS ----
## **********

## 1) R2 Histogram - Sampling event-based regression vs Bayesian regresion R2 ----
tiff(
  filename = paste0("./out_figs/", "Merc_SampleEventR2_vs_Stan", ".tif"),
  width = 10, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

hist(as.numeric(final_com$SER_adj_r2),
  ylim = c(0, 400),
  xlab = expression(R^2 ~ "(sampling event regressions)"),
  main = ""
)
abline(v = 0.80, col = "red", lwd = 2)
legend("topleft", legend = expression("Bayesian" ~ R^2), lwd = 2, bty = "n", col = "red")

t.test(as.numeric(final_com$SER_adj_r2), mu = 0.80)

dev.off()

sum(rowSums(table(m3$Hg_waterbody_code, m3$Hg_species_name) > 4)) # 1384 waterbody-species combinations
sum(rowSums(table(m3$Hg_waterbody_code, m3$Hg_species_name) > 1)) # 1730 waterbody-species combinations

## 2) Prediction intervals - Sampling event-based regresion vs. INLA ----
pred_width <- with(final_com, as.numeric(SER_upr) - as.numeric(SER_lwr))
bayes_width <- with(final_com, Inla_0_975 - Inla_0_025)

min(pred_width)
max(pred_width)
min(bayes_width)

hist(bayes_width, xlim = c(0, 10), ylim = c(0, 2000), breaks = seq(0, 10, by = 0.5), col = adjustcolor("lightgrey", alpha.f = 0.5), main = "", xlab = "Prediction interval width (Hg concentration)")
hist(pred_width, xlim = c(0, 10), ylim = c(0, 2000), breaks = seq(0, 35, by = 0.5), col = adjustcolor("lightblue", alpha.f = 0.5), add = T)
arrows(
  x0 = median(bayes_width), y0 = 0, x1 = median(bayes_width), y1 = 2000,
  angle = 0, length = 0, lwd = 4, col = "lightgrey"
)
arrows(
  x0 = median(pred_width), y0 = 0, x1 = median(pred_width), y1 = 2000,
  angle = 0, length = 0, lwd = 4, col = "lightblue"
)

t.test(bayes_width, pred_width) # bayes width significantly smaller

## 3) INLA species Hg~Length relationships ----

m_len_reg <- lapply(unique(m3$Hg_species_name), function(x) {
  print(x)
  (min_len <- min(subset(m3, Hg_species_name == x)$Hg_length_mm_log))
  (max_len <- max(subset(m3, Hg_species_name == x)$Hg_length_mm_log))

  (len_seq <- seq(round(min_len), round(max_len), by = 0.01))

  head(m3)

  pred_data <- data_frame(
    Hg_value_log = NA,
    Hg_length_mm_log = len_seq,
    Hg_waterbody_code = 44507658,
    Hg_species_name = x,
    Hg_species_code = NA,
    Hg_fishcode_lake = NA,
    Hg_sample_year = 2015
  )
  pred_data$YRWB <- paste(pred_data$Hg_waterbody_code, pred_data$Hg_sample_year)

  inla_sppred <- rbind(m3, pred_data)

  mod_full_length_inla_sppred <- inla(Hg_value_log ~ Hg_length_mm_log *
    Hg_species_name + f(Hg_waterbody_code, model = "iid") +
    f(YRWB, model = "iid"),
  data = inla_sppred,
  control.predictor = list(
    compute = TRUE,
    quantiles = c(0.025, 0.5, 0.975)
  ),
  control.compute = list(
    cpo = TRUE
  ), verbose = FALSE
  )
  mod_full_length_inla_sppred # 31 seconds

  posterior_est_sppred_inla <- mod_full_length_inla_sppred$summary.fitted.values[
    (nrow(m3) + 1):nrow(inla_sppred),
    c("0.025quant", "0.5quant", "0.975quant")
  ]

  nrow(posterior_est_sppred_inla)
  nrow(pred_data)

  pred_data$post_med <- posterior_est_sppred_inla[, 2]
  pred_data$post_2p5 <- posterior_est_sppred_inla[, 1]
  pred_data$post_97p5 <- posterior_est_sppred_inla[, 3]

  plot(exp(pred_data$post_med) ~ exp(pred_data$Hg_length_mm_log), type = "l", main = x)
  lines(exp(pred_data$post_2p5) ~ exp(pred_data$Hg_length_mm_log), type = "l", lty = 2)
  lines(exp(pred_data$post_97p5) ~ exp(pred_data$Hg_length_mm_log), type = "l", lty = 2)

  return(pred_data)
})

tiff(
  filename = paste0("./out_figs/", "Merc_SpeciesHg-LengthRelationships_INLA", ".tif"),
  width = 50, height = 25, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

windows(30, 20)
par(mfrow = c(5, 6))

lapply(m_len_reg, function(x) {
  plot(exp(x$post_med) ~ exp(x$Hg_length_mm_log),
    type = "b", main = x[1, "Hg_species_name"],
    lwd = 2, ylim = c(0, 5), ylab = "[Hg]", xlab = "Length (mm)", xlim = c(0, 1000), cex = 0.5
  )
  lines(exp(x$post_2p5) ~ exp(x$Hg_length_mm_log), type = "l", lty = 2)
  lines(exp(x$post_97p5) ~ exp(x$Hg_length_mm_log), type = "l", lty = 2)
})

dev.off()

## 4) Data recovery Stan model ----

tiff(
  filename = paste0("./out_figs/", "Merc_BayesianObsvsEst_Stan", ".tif"),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 1))

log_axes <- log(c(0.005, 0.05, 0.50, 5, 15))

with(fits, plot(Hg_value_log ~ fits$stan_posterior_q50,
  pch = 21, bg = "lightgrey",
  xlim = c(min(log_axes), max(log_axes)), ylim = c(min(log_axes), max(log_axes)),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
))

points(fits$stan_posterior_q50, fits$Hg_value_log,
  pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5), cex = 0.75
)

## Confidence interval on the relationship
min(fits$stan_posterior_q50)
max(fits$stan_posterior_q50)

pred_capac_full_stan <- lm(Hg_value_log ~ stan_posterior_q50, data = fits)

pred_capac_seq <- seq(-4.5, 1.5, by = 0.1)
pred_capac_prediction <- predict(pred_capac_full_stan,
  newdata = data.frame(stan_posterior_q50 = pred_capac_seq),
  level = 0.95, interval = "confidence"
)
pred_capac_prediction <- as.data.frame(pred_capac_prediction)
with(pred_capac_prediction, lines(fit ~ pred_capac_seq, lwd = 2))
summary(pred_capac_prediction)

axis(1, at = log_axes, labels = F, tck = 0.025, lwd.ticks = 1)
x_labs <- formatC(exp(log_axes), format = "f", digits = 3)
axis(1, at = log_axes, labels = x_labs, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75)

axis(2, at = log_axes, labels = F, tck = 0.025, lwd.ticks = 1, pos = log_axes[1])
y.labs <- formatC(exp(log_axes), format = "f", digits = 3)
axis(2, at = log_axes, labels = x_labs, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = log_axes[1], las = 2)

title(xlab = expression("Predicted Hg (" * mu * "g g"^-1 * ")"), line = 1.75)
title(ylab = expression("Observed Hg (" * mu * "g g"^-1 * ")"), line = 1.75)

dev.off()

## 5) Data recovery INLA model ----

tiff(
  filename = paste0("./out_figs/", "Merc_BayesianObsvsEst_INLA", ".tif"),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 1))

log_axes <- log(c(0.005, 0.05, 0.50, 5, 15))

with(fits, plot(Hg_value_log ~ inla_posterior_q50,
  pch = 21, bg = "lightgrey",
  xlim = c(min(log_axes), max(log_axes)), ylim = c(min(log_axes), max(log_axes)),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
))

points(fits$inla_posterior_q50, fits$Hg_value_log,
  pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5), cex = 0.75
)

## Confidence interval on the relationship
min(fits$inla_posterior_q50)
max(fits$inla_posterior_q50)

pred_capac_full_inla <- lm(Hg_value_log ~ inla_posterior_q50, data = fits)

pred_capac_seq <- seq(-4.5, 1.5, by = 0.1)
pred_capac_prediction <- predict(pred_capac_full_inla,
  newdata = data.frame(inla_posterior_q50 = pred_capac_seq),
  level = 0.95, interval = "confidence"
)
pred_capac_prediction <- as_data_frame(pred_capac_prediction)
with(pred_capac_prediction, lines(fit ~ pred_capac_seq, lwd = 2))
summary(pred_capac_prediction)

axis(1, at = log_axes, labels = F, tck = 0.025, lwd.ticks = 1)
x.labs <- formatC(exp(log_axes), format = "f", digits = 3)
axis(1, at = log_axes, labels = x.labs, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75)

axis(2, at = log_axes, labels = F, tck = 0.025, lwd.ticks = 1, pos = log_axes[1])
y.labs <- formatC(exp(log_axes), format = "f", digits = 3)
axis(2, at = log_axes, labels = x.labs, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = log_axes[1], las = 2)

title(xlab = expression("Predicted Hg (" * mu * "g g"^-1 * ")"), line = 1.75)
title(ylab = expression("Observed Hg (" * mu * "g g"^-1 * ")"), line = 1.75)

dev.off()

## 6) Stan vs. INLA regression coefficient comparison ----

inla_coefs <- mod_full_length_inla_summary
stan_coefs <- mod_full_length_stan_summary

inla_coefs[1:28, ]
inla_coefs[2, ]

inla_ints <- inla_coefs[1:28, ]

head(inla_ints)

inla_ints[, c(2, 4, 5, 6, 7)] <- lapply(inla_ints[, c(2, 4, 5, 6, 7)], function(x) {
  x[-c(1, 2)] <- x[-c(1, 2)] + x[1]
  x
})

names(stan_coefs)

stan_ints <- stan_coefs[1:28, ]
# stan.ints[,c(1:3)] <- lapply(stan.ints[,c(1:3)], function(x){
#
#  x[-c(1,2)] <- x[-c(1,2)] + x[1]
#  x
# })

tiff(
  filename = paste0("./out_figs/", "Merc_FishMercModel_SpeciesIntercepts", ".tif"),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(8, 3, 3, 1))
inla_ints

plot(1:28, inla_ints$`0.5quant`[1:28],
  pch = 21, bg = "lightgrey",
  xlim = c(0, 30), ylim = c(-40, 25),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
)

arrows(
  x0 = 1:28 - 0.125, y0 = inla_ints$`0.025quant`[1:28],
  x1 = 1:28 - 0.125, y1 = inla_ints$`0.975quant`[1:28], angle = 0, length = 0
)
points(1:28 - 0.125, inla_ints$`0.5quant`[1:28], pch = 16, cex = 0.5)

arrows(
  x0 = 1:28 + 0.125, y0 = stan_ints$Q2.5[1:28],
  x1 = 1:28 + 0.125, y1 = stan_ints$Q97.5[1:28], angle = 0, col = "red", length = 0
)
points((1:28) + 0.125, stan_ints$Q50[1:28], pch = 16, col = "red", cex = 0.5)

x_labels <- gsub("Hg_species_name", "", mod_full_length_stan_summary$param[1:28])
x_labels[2] <- "Length"
x_labels[1] <- "Intercept"
axis(1, at = 1:28, labels = F, tck = 0.025, lwd.ticks = 1, las = 2)
axis(1, at = 1:28, labels = x_labels, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, las = 2)

y_labels <- formatC((seq(-40, 25, by = 5)), digits = 3)
y_labels
axis(2, at = seq(-40, 25, by = 5), labels = F, tck = 0.025, lwd.ticks = 1, pos = 0)
axis(2, at = seq(-40, 25, by = 5), labels = y_labels, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = 0, las = 2)

title(xlab = expression("Species"), line = 7)
title(ylab = expression("Intercept (log (" * mu * "g g"^-1 * ")"), line = 1.75)

legend("top",
  legend = c("Stan", "R-INLA"), bty = "n", pch = 16, col = c("black", "red"),
  horiz = TRUE
)

dev.off()

## 7) Stan vs. INLA interaction regression comparison ----

inla_coefs <- mod_full_length_inla_summary
stan_coefs <- mod_full_length_stan_summary

inla_coefs$ID[29:54]
stan_coefs$param[29:54]

## Select those estimated by each program
inla_coefs_sub <- inla_coefs[29:54, ][inla_coefs$ID[29:54] %in% stan_coefs$param, ]
stan_coefs_sub <- stan_coefs[29:54, ][stan_coefs$param[29:54] %in% inla_coefs_sub$ID, ]

stan_coefs_sub$param
inla_coefs_sub$ID

identical(stan_coefs_sub$param, inla_coefs_sub$ID)

tiff(
  filename = paste0("./out_figs/", "Merc_FishMercModel_SpeciesLengthInteractions", ".tif"),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(8, 3, 3, 1))

plot(1:25, inla_coefs_sub$`0.5quant`,
  pch = 21, bg = "lightgrey",
  xlim = c(0, 26), ylim = c(-5, 5),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
)

lines(seq(1:25), y = rep(0, 25), lty = 2)

arrows(
  x0 = 1:25 - 0.125, y0 = inla_coefs_sub$`0.025quant`,
  x1 = 1:25 - 0.125, y1 = inla_coefs_sub$`0.975quant`, angle = 0, length = 0
)
points(1:25 - 0.125, inla_coefs_sub$`0.5quant`, pch = 16, cex = 0.5)

arrows(
  x0 = 1:25 + 0.125, y0 = stan_coefs_sub$Q2_5,
  x1 = 1:25 + 0.125, y1 = stan_coefs_sub$Q97_5, angle = 0, col = "red", length = 0
)
points((1:25) + 0.125, stan_coefs_sub$Q50, pch = 16, col = "red", cex = 0.5)

x_labels <- gsub("Hg_length_mm_log:Hg_species_name", "", stan_coefs_sub$param)
axis(1, at = 1:25, labels = F, tck = 0.025, lwd.ticks = 1, las = 2)
axis(1, at = 1:25, labels = x_labels, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, las = 2)

y_labels <- formatC((seq(-5, 5, by = 2.5)), digits = 3)
y_labels
axis(2, at = seq(-5, 5, by = 2.5), labels = F, tck = 0.025, lwd.ticks = 1, pos = 0)
axis(2, at = seq(-5, 5, by = 2.5), labels = y_labels, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = 0, las = 2)

title(xlab = expression("Species"), line = 7)
title(ylab = expression("Interaction coefficient (slope + __)"), line = 1.75)

legend("top",
  legend = c("STAN", "R-INLA"), bty = "n", pch = 16, col = c("black", "red"),
  horiz = TRUE
)

dev.off()

## 8) Stan vs INLA random prediction estimates ----

tiff(
  filename = paste0("./out_figs/", "Merc_FishMercModel_StanVsINLAPredictions", ".tif"),
  width = 30, height = 15, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(8, 5, 3, 1))

## INLA
boxplot(predictions$`50%` ~ predictions$Hg_species_name,
  las = 2,
  ylab = expression("Predicted Hg (" * mu * "g g"^-1 * ")"),
  xlab = "", cex.axis = 0.75,
  main = "Stan", ylim = c(-12, 5)
)

boxplot(predictions$`0.5quant` ~ predictions$Hg_species_name,
  las = 2,
  ylab = expression("Predicted Hg (" * mu * "g g"^-1 * ")"),
  xlab = "", cex.axis = 0.75,
  main = "INLA", ylim = c(-12, 5)
)

dev.off()

## 9) Spatial predictions by species ----

## Load administrative boundary layer and subset by Ontario

st_layers("D:/CEON/spatialraw/canvec_1M_CA_Admin_fgdb/canvec_1M_CA_Admin.gdb")

admin <- st_read(
  dsn = file.path(
    sp_dir,
    "spatialraw/canvec_1M_CA_Admin_fgdb/canvec_1M_CA_Admin.gdb"
  ),
  layer = "geo_political_region_2"
)
ont <- subset(admin, jurisdiction == 100)
ont <- st_transform(ont, st_crs(m1))
g <- st_graticule(ont)

## 9a) Northern Pike ----

tiff(paste0(
  "./out_figs/",
  "Merc_HgPredictions_NorthernPike",
  Sys.Date(), ".tif"
),
height = 15, width = 15, units = "cm",
compression = "lzw", res = 300
)

plot(st_geometry(predictions_spatial_np),
  pch = 16, cex = 1,
  col = "white"
)
plot(st_geometry(ont),
  col = NA,
  border = "lightgrey", lwd = 3, axes = FALSE, bty = "n", add = T
)
plot(g, add = T, col = "lightgrey")
plot(st_geometry(predictions_spatial_np),
  pch = 16, cex = 1,
  col = adjustcolor(predictions_spatial_np$jenks_col, alpha.f = 0.75), add = T
)
legend("top",
  legend = formatC(jenks_breaks$brks, digits = 2),
  pch = 16, col = viridis(6), bty = "n", horiz = T
)

dev.off()

## 9b) Walleye ----

tiff(paste0(
  "./out_figs/",
  "Merc_HgPredictions_Walleye",
  Sys.Date(), ".tif"
),
height = 15, width = 15, units = "cm",
compression = "lzw", res = 300
)

plot(st_geometry(predictions_spatial_wa),
  pch = 16, cex = 1,
  col = "white"
)
plot(st_geometry(ont),
  col = NA,
  border = "lightgrey", lwd = 3, axes = FALSE, bty = "n", add = T
)
plot(g, add = T, col = "lightgrey")
plot(st_geometry(predictions_spatial_wa),
  pch = 16, cex = 1,
  col = adjustcolor(predictions_spatial_wa$jenks_col, alpha.f = 0.75), add = T
)
legend("top",
  legend = formatC(jenks_breaks$brks, digits = 2),
  pch = 16, col = viridis(6), bty = "n", horiz = T
)

dev.off()

## 9c) Smallmouth Bass ----

tiff(paste0(
  "./out_figs/",
  "Merc_HgPredictions_SmallmouthBass",
  Sys.Date(), ".tif"
),
height = 15, width = 15, units = "cm",
compression = "lzw", res = 300
)

plot(st_geometry(predictions_spatial_sb),
  pch = 16, cex = 1,
  col = "white"
)
plot(st_geometry(ont),
  col = NA,
  border = "lightgrey", lwd = 3, axes = FALSE, bty = "n", add = T
)
plot(g, add = T, col = "lightgrey")
plot(st_geometry(predictions_spatial_sb),
  pch = 16, cex = 1,
  col = adjustcolor(predictions_spatial_sb$jenks_col, alpha.f = 0.75), add = T
)
legend("top",
  legend = formatC(jenks_breaks$brks, digits = 2),
  pch = 16, col = viridis(6), bty = "n", horiz = T
)

dev.off()

## 9d) Walleye and Northern Pike combined ----

tiff(paste0(
  "./out_figs/",
  "Merc_HgPredictions_NorthernPike",
  Sys.Date(), ".tif"
),
height = 15, width = 15, units = "cm",
compression = "lzw", res = 300
)

par(mfrow = c(2, 1))
plot(st_geometry(predictions_spatial_np),
  pch = 16, cex = 1,
  col = adjustcolor(predictions_spatial_np$jenks_col, alpha.f = 0.75)
)
plot(st_geometry(ont),
  col = adjustcolor("lightgrey", alpha = 0.75),
  border = "lightgrey", lwd = 3, axes = FALSE, bty = "n", add = T
)
plot(g, add = T, col = "lightgrey")


plot(st_geometry(predictions_spatial_wa),
  pch = 16, cex = 1,
  col = adjustcolor(predictions_spatial_wa$jenks_col, alpha.f = 0.75)
)
plot(st_geometry(ont),
  col = adjustcolor("lightgrey", alpha = 0.75),
  border = "lightgrey", lwd = 3, axes = FALSE, bty = "n", add = T
)
plot(g, add = T, col = "lightgrey")
legend("top",
  legend = formatC(jenks.breaks$brks, digits = 2),
  pch = 16, col = viridis(6), bty = "n", horiz = T
)




dev.off()



## 10) MEM importance ----

## 10a) Northern Pike ----
tiff(
  filename = paste0(
    "./out_figs/", "Merc_MEMVariableImportance_NorthernPike",
    Sys.Date(), ".tif"
  ),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 1))

plot(rf2$importance[, 1],
  pch = 21, bg = "lightgrey",
  xlim = c(-10, 140), ylim = c(0, 0.003),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
)

points(rf2$importance[, 1], type = "l", lwd = 2, col = "lightgrey")

axis(1, at = seq(-10, 140, by = 10), labels = F, tck = 0.025, lwd.ticks = 1)
axis(1, at = seq(0, 140, by = 10), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75)

axis(2, at = seq(0, 0.003, by = 0.001), labels = F, tck = 0.025, lwd.ticks = 1, pos = -10)
axis(2, at = seq(0, 0.003, by = 0.001), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = -10)

title(xlab = expression("MEM"), line = 1.75)
title(ylab = expression("% Increase in MSE"), line = 1.75)

dev.off()

## 10b) Walleye ----
tiff(
  filename = paste0(
    "./out_figs/", "Merc_MEMVariableImportance_Walleye_",
    Sys.Date(), ".tif"
  ),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 1))

plot(rf2$importance[, 1],
  pch = 21, bg = "lightgrey",
  xlim = c(-10, 140), ylim = c(0, 0.006),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
)

points(rf2$importance[, 1], type = "l", lwd = 2, col = "lightgrey")

axis(1, at = seq(-10, 140, by = 10), labels = F, tck = 0.025, lwd.ticks = 1)
axis(1, at = seq(0, 140, by = 10), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75)

axis(2, at = seq(0, 0.006, by = 0.001), labels = F, tck = 0.025, lwd.ticks = 1, pos = -10)
axis(2, at = seq(0, 0.006, by = 0.001), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = -10)

title(xlab = expression("MEM"), line = 1.75)
title(ylab = expression("% Increase in MSE"), line = 1.75)

dev.off()

## 10b) SmallmouthBass ----
tiff(
  filename = paste0(
    "./out_figs/", "Merc_MEMVariableImportance_SmallmouthBass_",
    Sys.Date(), ".tif"
  ),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 1))

plot(rf2$importance[, 1],
  pch = 21, bg = "lightgrey",
  xlim = c(-10, 140), ylim = c(0, 0.006),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
)

points(rf2$importance[, 1], type = "l", lwd = 2, col = "lightgrey")

axis(1, at = seq(-10, 140, by = 10), labels = F, tck = 0.025, lwd.ticks = 1)
axis(1, at = seq(0, 140, by = 10), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75)

axis(2, at = seq(0, 0.006, by = 0.001), labels = F, tck = 0.025, lwd.ticks = 1, pos = -10)
axis(2, at = seq(0, 0.006, by = 0.001), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = -10)

title(xlab = expression("MEM"), line = 1.75)
title(ylab = expression("% Increase in MSE"), line = 1.75)

dev.off()

## 11) Ontario lake distribution ----

lake_points <- group_by(predictions_spatial, Hg_waterbody_code)
lake_points <- summarise(predictions_spatial, unique(Hg_waterbody_code))

tiff(paste0(
  "./out_figs/",
  "Merc_LakeDistribution",
  Sys.Date(), ".tif"
),
height = 15, width = 15, units = "cm",
compression = "lzw", res = 300
)

plot(st_geometry(lake_points), pch = 16, cex = 1, col = "white")
plot(st_geometry(ont),
  col = adjustcolor("lightgrey", alpha = 0.75),
  border = "lightgrey", lwd = 3, axes = FALSE, bty = "n", add = T
)
plot(g, add = T, col = "lightgrey")
plot(st_geometry(lake_points),
  pch = 16, cex = 1, add = T,
  col = adjustcolor("black", alpha.f = 0.5)
)

dev.off()




## POTENTIAL JUNK ----

## Plot: Species Intercepts by rank ----

mod.full.sims <- as.data.frame(as.matrix(mod.full))

sel.cols <- colnames(mod.full.sims)[grep("Hg_species_name", colnames(mod.full.sims))]
sel.cols <- sel.cols[!grepl("length_mm", sel.cols)]
sel.cols

mod.full.sims

a_quant <- apply(
  X = mod.full.sims[sel.cols],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2.5", "Q50", "Q97.5")
a_quant$species <- gsub("Hg_species_name", "", sel.cols)

a_quant <- a_quant[with(a_quant, order(Q50)), ]

a_quant$signif <- sapply(1:nrow(a_quant), function(x) {
  (a_quant[x, "Q2.5"] > 0) & (a_quant[x, "Q97.5"] > 0) |
    (a_quant[x, "Q2.5"] < 0) & (a_quant[x, "Q97.5"] < 0)
})

min(as.matrix(a_quant[, 1:3]))
max(as.matrix(a_quant[, 1:3]))

tiff(
  filename = paste0("./out_figs/", "Merc_FishMercModel_SpeciesIntercepts", ".tif"),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(8, 3, 3, 1))

with(a_quant, plot(Q50,
  pch = 21, bg = "lightgrey",
  ylim = c(-8, 7),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
))

arrows(
  x0 = seq(1:nrow(a_quant)),
  x1 = seq(1:nrow(a_quant)),
  y0 = a_quant$Q2.5,
  y1 = a_quant$Q97.5, angle = 0
)

points(a_quant$Q50, pch = 16, cex = 1.5, col = ifelse(a_quant$signif, "red", "black"))


axis(1, at = seq(1:nrow(a_quant)), labels = F, tck = 0.025, lwd.ticks = 1)
axis(1, at = seq(1:nrow(a_quant)), labels = a_quant$species, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, las = 2)

axis(2, at = seq(-8, 7, by = 1), labels = F, tck = 0.025, lwd.ticks = 1, pos = 0)
axis(2, at = seq(-8, 7, by = 1), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = 0, las = 2)

title(xlab = expression("Species"), line = 7)
title(ylab = expression(Delta ~ "Intercept (log (" * mu * "g g"^-1 * ")"), line = 1.75)

dev.off()


expression("Observed Hg (" * mu * "g g"^-1 * ")")

## Plot: Species interactions by rank ----

sel.cols <- colnames(mod.full.sims)[grep("Hg_species_name", colnames(mod.full.sims))]
sel.cols <- sel.cols[grepl("length_mm", sel.cols)]
sel.cols

a_quant <- apply(
  X = mod.full.sims[sel.cols],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.50, 0.975)
)
a_quant <- data.frame(t(a_quant))
names(a_quant) <- c("Q2.5", "Q50", "Q97.5")
a_quant$species <- gsub("Hg_length_mm_log:Hg_species_name", "", sel.cols)

a_quant <- a_quant[with(a_quant, order(Q50)), ]

a_quant$signif <- sapply(1:nrow(a_quant), function(x) {
  (a_quant[x, "Q2.5"] > 0) & (a_quant[x, "Q97.5"] > 0) |
    (a_quant[x, "Q2.5"] < 0) & (a_quant[x, "Q97.5"] < 0)
})

min(as.matrix(a_quant[, 1:3]))
max(as.matrix(a_quant[, 1:3]))

tiff(
  filename = paste0("./out_figs/", "Merc_FishMercModel_SpeciesInterceptsInteractions", ".tif"),
  width = 15, height = 10, units = "cm", res = 300,
  compression = "lzw+p", type = "cairo"
)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(8, 3, 3, 1))

with(a_quant, plot(Q50,
  pch = 21, bg = "lightgrey",
  ylim = c(-1.5, 1.5),
  type = "n", xaxt = "n", yaxt = "n", yaxs = "i",
  xlab = "", ylab = "", main = "", frame.plot = F
))

arrows(
  x0 = seq(1:nrow(a_quant)),
  x1 = seq(1:nrow(a_quant)),
  y0 = a_quant$Q2.5,
  y1 = a_quant$Q97.5, angle = 0
)

points(a_quant$Q50, pch = 16, cex = 1.5, col = ifelse(a_quant$signif, "red", "black"))


axis(1, at = seq(1:nrow(a_quant)), labels = F, tck = 0.025, lwd.ticks = 1)
axis(1, at = seq(1:nrow(a_quant)), labels = a_quant$species, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, las = 2)

axis(2, at = seq(-1.5, 1.5, by = 0.5), labels = F, tck = 0.025, lwd.ticks = 1, pos = 0)
axis(2, at = seq(-1.5, 1.5, by = 1), labels = T, tick = F, lwd = 0, cex.axis = 0.75, line = -0.75, pos = 0, las = 2)



title(xlab = expression("Species"), line = 7)
title(ylab = expression("Intercept"), line = 1.75)

dev.off()
















## ************
## Plot spatial
## ************
names(m1)

m2_length <- merge(sbf_sport_preds_length, m1[, c(
  "waterbody_code", "sample_year",
  "species_code", "latitude_ddmmss",
  "longitude_ddmmss", "species_name"
)],
by = c("waterbody_code", "sample_year", "species_code"), all.x = TRUE
)

m2_mass <- merge(sbf_sport_preds_mass, m1[, c(
  "waterbody_code", "sample_year",
  "species_code", "latitude_ddmmss",
  "longitude_ddmmss", "species_name"
)],
by = c("waterbody_code", "sample_year", "species_code"), all.x = TRUE
)

m2_cond <- merge(sbf_sport_preds_cond, m1[, c(
  "waterbody_code", "sample_year",
  "species_code", "latitude_ddmmss",
  "longitude_ddmmss", "species_name"
)],
by = c("waterbody_code", "sample_year", "species_code"), all.x = TRUE
)

## Remove duplicates
m2_length <- unique(m2_length)
m2_mass <- unique(m2_mass)
m2_cond <- unique(m2_cond)

## Put in list
m2 <- list(m2_length, m2_mass, m2_cond)

## Convert to sf
m2.space <- lapply(m2, function(x) {

  ## Convert degrees minutes seconds to decimal degrees for lat and long
  dd.lat <- substring(x$latitude_ddmmss, 1, 2)
  mm.lat <- substring(x$latitude_ddmmss, 3, 4)
  ss.lat <- substring(x$latitude_ddmmss, 5, 6)

  spat.lat <- paste0(dd.lat, "d", mm.lat, "'", ss.lat, "\"")

  decdeg.lat <- char2dms(spat.lat)
  decdeg.lat <- as.numeric(decdeg.lat)

  dd.long <- substring(x$longitude_ddmmss, 1, 2)
  mm.long <- substring(x$longitude_ddmmss, 3, 4)
  ss.long <- substring(x$longitude_ddmmss, 5, 6)

  spat.long <- paste0(dd.long, "d", mm.long, "'", ss.long, "\"")

  decdeg.long <- char2dms(spat.long)
  decdeg.long <- as.numeric(decdeg.long) * -1 # negative adjusts for W

  ## Assuming NAD83
  x$latitude_decdeg <- decdeg.lat
  x$longitude_decdeg <- decdeg.long

  # plot(x$latitude_decdeg ~ x$longitude_decdeg)

  x <- st_as_sf(
    x = x, coords = c("longitude_decdeg", "latitude_decdeg"),
    crs = 4617
  )
  x <- st_transform(x = x, crs = 3161)

  x$fit <- as.numeric(x$fit)

  return(x)
})

par(mfrow = c(3, 3))

plot(m2.space[[1]]["fit"], breaks = "quantile", pch = 16, pal = viridis::viridis, key.pos = 1)
plot(m2.space[[2]]["fit"], breaks = "quantile", pch = 16, pal = viridis::viridis, key.pos = 1)
plot(m2.space[[3]]["fit"], breaks = "quantile", pch = 16, pal = viridis::viridis, key.pos = 1)

hist((m2.space[[1]]$fit))
hist((m2.space[[2]]$fit))
hist((m2.space[[3]]$fit))

## Relationship between mean concentration and standardized concentration for
## all data
mean.conc <- aggregate(value ~ waterbody_code + sample_year + species_code +
  species_name, data = m1, FUN = mean)
length.conc <- aggregate(fit ~ waterbody_code + sample_year + species_code +
  species_name, data = m2.space[[1]], FUN = mean)
mass.conc <- aggregate(fit ~ waterbody_code + sample_year + species_code +
  species_name, data = m2.space[[2]], FUN = mean)
cond.conc <- aggregate(fit ~ waterbody_code + sample_year + species_code +
  species_name, data = m2.space[[3]], FUN = mean)

## Merge the datasets
conc.merged <- merge(length.conc, mean.conc)
conc.merged <- merge(conc.merged, mass.conc)

mean.length.merged <- merge(length.conc, mean.conc)