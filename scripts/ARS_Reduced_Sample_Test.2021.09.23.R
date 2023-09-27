## ***********************************************
## INLA mixed models for size standardized arsenic
## ***********************************************

## Run INLA mixed models and predict [As] values at the lake-level while removing observations to observe changes in posterior quantiles


## **********
## SETUP ----
## **********

## Load packages
library(INLA)
library(dplyr)
library(readxl)
library(glmmTMB)
library(effects)
library(sjPlot)
library(ggplot2)

#install.packages("INLA", repos = c(getOption("repos"),
#                                   INLA = "https://inla.r-inla-download.org/R/testing"), dep = TRUE)

## *******************
## ANALYSIS SETUP ----
## *******************

## Read the saved file
#setwd("")
fish <- read.csv("./data/Fish_As_2021.09.16.csv")

fish_lake <- subset(fish, System_Type == "Lake")

## Some data cleanup and variable creation 
fish_lake <- dplyr::rename(fish_lake, As_ugg = As_ug_ww)
fish_lake$Sampling_Year <- as.factor(as.numeric(sapply(strsplit(fish_lake$Sampling_Date, split = "-"), function(x){x[[1]]})))
fish_lake$As_ugg_log <- log(fish_lake$As_ugg)
fish_lake$RWT_log <- log(fish_lake$RWT)
fish_lake$TLEN_log <- log(fish_lake$TLEN)

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

with(fish_lake, plot(log(As_ugg) ~ log(TLEN)))
with(fish_lake, plot(log(As_ugg) ~ log(RWT))) 


## *****************************************
## WALLEYE MIXED MODEL REGRESSION, INLA ----
## *****************************************

## Slightly more complicated to set up but very similar results and better estimates on the variance components 
## which I feel is an important diagnostic 

## Walleye
WE <- subset(fish_lake, Taxon == "WALL") 
WE$WATERBODY_CODE <- factor(WE$Waterbody) # for inla model
WE$WATERBODY_CODE1 <- as.integer(WE$WATERBODY_CODE) # for inla model
WE$WATERBODY_CODE2 <- WE$WATERBODY_CODE1 + max(WE$WATERBODY_CODE1) # for inla model
WE_n_waterbody <- dplyr::n_distinct(WE$WATERBODY_CODE) # for inla model

## Set up precision --> standard deviation formula; Bayesian models use precision (tau) where sd = 1/sqrt(tau) 
MySqrt <- function(x) {
        1 / sqrt(x)
}


## **********
## + RWT ----
## **********

## Run model, same model as GLMMTMB as above
INLA_WE_RWT <- inla(As_ugg_log ~ RWT_log + 
                            f(WATERBODY_CODE1, n = 2 * WE_n_waterbody, model = "iid2d") +  
                            f(WATERBODY_CODE2, RWT_log, copy = "WATERBODY_CODE1"),
                    data = WE, 
                    control.predictor = list(
                            compute = TRUE, 
                            quantiles = c(0.025, 0.5, 0.975)
                    ),
                    control.compute = list(
                            cpo = TRUE
                    )
)

summary(INLA_WE_RWT)

## There is a lot of information in the INLA objects 
## Below, I show how to extract the pertinent information for the models 

## Fitted model
INLA_WE_RWT_fits <- data.frame(WATERBODY_CODE = WE$WATERBODY_CODE, 
                               As_ugg_log = WE$As_ugg_log)

## Alternative comparison to fitted and residuals
INLA_WE_RWT_fits$inla_posterior_q50 <- INLA_WE_RWT$summary.fitted.values[, "0.5quant"]
INLA_WE_RWT_fits$inla_posterior_q2p5 <- INLA_WE_RWT$summary.fitted.values[, "0.025quant"]
INLA_WE_RWT_fits$inla_posterior_q97p5 <- INLA_WE_RWT$summary.fitted.values[, "0.975quant"]
INLA_WE_RWT_fits$resid_inla <- INLA_WE_RWT_fits$As_ugg_log - INLA_WE_RWT_fits$inla_posterior_q50

par(mfrow = c(1,3))
## Residuals vs. fit on linear scale
lw1_inla <- loess(resid_inla ~ inla_posterior_q50, data = INLA_WE_RWT_fits)
j <- order(INLA_WE_RWT_fits$inla_posterior_q50)
plot(resid_inla ~ inla_posterior_q50, data = INLA_WE_RWT_fits, pch = 21, bg = adjustcolor("lightgrey", alpha.f = 0.5))
lines(lw1_inla$fitted[j] ~ INLA_WE_RWT_fits$inla_posterior_q50[j], lwd = 2)
## Looks about the same as the other models but some poor fit at higher values

## R2
plot(INLA_WE_RWT_fits$As_ugg_log ~ INLA_WE_RWT_fits$inla_posterior_q50, pch = 16) 
fitted_test <- lm(INLA_WE_RWT_fits$As_ugg_log ~ INLA_WE_RWT_fits$inla_posterior_q50)
summary(fitted_test) ## model explains about 75% of the variation in As 
abline(fitted_test)

## Histogram of residuals
hist(INLA_WE_RWT_fits$resid_inla) # approximately normally distributed

## Generate parameter estimates from posterior distribution for fixed and random effects and manipulate them into a single dataframe
INLA_WE_RWT_fixed <- data.frame(
        ID = rownames(INLA_WE_RWT$summary.fixed),
        INLA_WE_RWT$summary.fixed, stringsAsFactors = FALSE
)
names(INLA_WE_RWT_fixed) <- c("ID", names(INLA_WE_RWT$summary.fixed))
INLA_WE_RWT_fixed$Type <- "Fixed"
head(INLA_WE_RWT_fixed)

INLA_WE_RWT_random_intercept <- INLA_WE_RWT$summary.random$WATERBODY_CODE1
INLA_WE_RWT_random_intercept$ID <- as.character(INLA_WE_RWT_random_intercept$ID)
INLA_WE_RWT_random_intercept$WATERBODY_CODE1 <- INLA_WE_RWT_random_intercept$ID
INLA_WE_RWT_random_intercept <- merge(INLA_WE_RWT_random_intercept, WE[,c("WATERBODY_CODE1", "Waterbody")], no.dups = TRUE)
INLA_WE_RWT_random_intercept <- INLA_WE_RWT_random_intercept[!duplicated(INLA_WE_RWT_random_intercept),]
INLA_WE_RWT_random_intercept$Type <- "Random Intercept - Waterbody"
head(INLA_WE_RWT_random_intercept)

INLA_WE_RWT_random_slope <- INLA_WE_RWT$summary.random$WATERBODY_CODE2
INLA_WE_RWT_random_slope$ID <- as.character(INLA_WE_RWT_random_slope$ID)
INLA_WE_RWT_random_slope$WATERBODY_CODE2 <- INLA_WE_RWT_random_slope$ID
INLA_WE_RWT_random_slope <- merge(INLA_WE_RWT_random_slope, WE[,c("WATERBODY_CODE2", "Waterbody")])
INLA_WE_RWT_random_slope <- INLA_WE_RWT_random_slope[!duplicated(INLA_WE_RWT_random_slope),]
INLA_WE_RWT_random_slope$Type <- "Random Slope - Waterbody"
head(INLA_WE_RWT_random_slope)

INLA_WE_RWT_summary <- list(
        INLA_WE_RWT_fixed,
        INLA_WE_RWT_random_intercept,
        INLA_WE_RWT_random_slope
)

INLA_WE_RWT_summary <- dplyr::bind_rows(INLA_WE_RWT_summary)
head(INLA_WE_RWT_summary)
tail(INLA_WE_RWT_summary)

## The estimated coefficient for the grand intercept is 
## -4.28 (-4.68, -3.89) using the 50th percentile from the posterior distribution and the 0.025 and 0.975 quantiles for the error 

## The estimated slope coefficient is 
## 0.14 (0.08, 0.20) using the same values 

## Variance parameters, difficult to read in this form. 

## Calculating the sd as in GLMMTMB (sigma) from the precision (tau) values 
INLA_WE_RWT$marginals.hyperpar$`Precision for the Gaussian observations`
INLA_WE_RWT$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
INLA_WE_RWT$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`

tau_WATERBODY_CODE1 <- INLA_WE_RWT$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 1)`
tau_WATERBODY_CODE1_RSLOPE <- INLA_WE_RWT$marginals.hyperpar$`Precision for WATERBODY_CODE1 (component 2)`
tau_residual <- INLA_WE_RWT$marginals.hyperpar$`Precision for the Gaussian observations`

## Generate standard deviations as in TMB
(sigma_WATERBODY_CODE1 <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1)) # not clear on how this function works?
(sigma_WATERBODY_CODE1_RSLOPE <- inla.emarginal(MySqrt, tau_WATERBODY_CODE1_RSLOPE))
(sigma_residual <- inla.emarginal(MySqrt, tau_residual))

## Generate variance components 
(WE_variance_components <- c(sigma_residual, sigma_WATERBODY_CODE1, sigma_WATERBODY_CODE1_RSLOPE)^2)


## ****************************
## WALLEYE PREDICTION TEST ----
## ****************************

## Predictions in INLA are made by appending your prediction information to the dataframe used to run the model 
## Except, the As_ugg values are set to NA
## Then you extract the posterior information 

WE_complete <- WE[complete.cases(WE), ]

WE_predict <- data.frame(Waterbody = unique(WE[,c("Waterbody")]))
WE_predict$TLEN_log <- median(WE_complete$TLEN_log)
WE_predict$RWT_log <- median(WE_complete$RWT_log)
WE_predict <- WE_predict[complete.cases(WE_predict),]

WE_predict$WATERBODY_CODE <- factor(WE_predict$Waterbody) 
WE_predict$WATERBODY_CODE1 <- as.integer(WE_predict$WATERBODY_CODE)   
WE_predict$WATERBODY_CODE2 <- WE_predict$WATERBODY_CODE1 + max(WE_predict$WATERBODY_CODE1) 

head(WE_predict) 


# WE_predict_inla contains one row for each waterbody, with As_ugg_log set to 'NA' - this will stay the same
WE_predict_inla <- merge(WE_predict, WE[,c("WATERBODY_CODE1", "WATERBODY_CODE2")], sort = FALSE)
WE_predict_inla <- WE_predict_inla[!duplicated(WE_predict_inla),]
WE_predict_inla$As_ugg_log <- NA
nrow(WE_predict)
nrow(WE_predict_inla)


## Make sure these are identical so that we can easily line up with the predictions from TMB
identical(WE_predict$Waterbody, WE_predict_inla$Waterbody)

# WE_inla_prediction is all walleye [As] observations which is used to inform the model - the model will then be used to predict lake-level As 
# values from 'WE_predict_inla'. 

WE_inla_prediction <- WE[names(WE_predict_inla)]

filter(as.data.frame(table(WE_inla_prediction$Waterbody)), Freq >= 10)

## Create loop that uses lakes with high n (>= 10), then progressively remove 1 observation at a time while rerunning the model and 
## appending the results in a single dataframe for comparison

## Select waterbody names with a sample size >= 10
test.lk <- filter(as.data.frame(table(WE_inla_prediction$Waterbody)), Freq >= 10)
test.lk <- as.character(test.lk$Var1)

## Run prediction models for lakes with highest sample sizes
WE_inla_test <- foreach(j = 1:length(test.lk)) %do% {
        
        WE_inla_test_lk <- foreach(i = 1:nrow(filter(WE_inla_prediction, Waterbody == test.lk[j]))-1) %do% {
        
                ## Randomly remove rows from the specified waterbody (repeating the process for 1:n-1 iterations)
                if(i >= 1){
                        lk.rm <- WE_inla_prediction[-sample(which(WE_inla_prediction$Waterbody == test.lk[j]), i), ] # subset does not work here when i = 0
                } else {
                        lk.rm <- WE_inla_prediction
                }
        
                ## Append lakes to be predicted (i.e., WE_predict_inla)
                WE_lk.rm <- rbind(lk.rm, WE_predict_inla)
                row.names(WE_lk.rm) <- 1:nrow(WE_lk.rm) # have to reset row names for indexing to work properly when extracting posterior information
        
                ## Run INLA model
                message("Running INLA prediction model for ", test.lk[j], " (lake ", j, " of ", nrow(test.lk), ") with ", i, " of ", nrow(filter(WE_inla_prediction, Waterbody == test.lk[j])), " observations removed")
        
                INLA_WE_RWT_prediction_model <- inla(As_ugg_log ~ RWT_log + 
                                                        f(WATERBODY_CODE1, n = 2 * WE_n_waterbody, model = "iid2d") +   
                                                        f(WATERBODY_CODE2, RWT_log, copy = "WATERBODY_CODE1"),
                                                data = WE_lk.rm, 
                                                control.predictor = list(
                                                        compute = TRUE, 
                                                        quantiles = c(0.025, 0.5, 0.975)
                                                ),
                                                control.compute = list(
                                                        cpo = TRUE
                                                )
                )
        
                ## Retain only predictions for specified lake
                lk.row <- subset(WE_lk.rm, Waterbody == test.lk[j])
        
                INLA_WE_RWT_predict_posteriors <- INLA_WE_RWT_prediction_model$summary.fitted.values[
                        max(as.numeric(row.names(lk.row))), # pulls a single row for each iteration; appended row to be predicted will always be the highest row number, therefore using max() function to index
                        c("0.025quant", "0.5quant", "0.975quant")
                ]
        
                ## Convert from log
                INLA_WE_RWT_predict_posteriors <- exp(INLA_WE_RWT_predict_posteriors)
        
                ## Combine predictions with lake info
                WE_lk_predict <- cbind(WE_lk.rm[max(as.numeric(row.names(lk.row))), ], INLA_WE_RWT_predict_posteriors) 
                WE_lk_predict <- WE_lk_predict[!duplicated(WE_lk_predict),]
        
                ## Add column to identify sample size
                WE_lk_predict <- mutate(WE_lk_predict, n = nrow(filter(WE_inla_prediction, Waterbody == test.lk[j])) - i)
        
                return(WE_lk_predict)
        
        }

        return(WE_inla_test_lk)

}

WE_inla_test <- bind_rows(WE_inla_test)

## Write csv for test results
#write.csv(WE_inla_test, "./data/INLA_size_standardized_data/reduced_sample_predicts/WE_inla_SampleSize_test.2021.09.22.csv", row.names = FALSE)


## *****************
## PLOT RESULTS ----
## *****************

#WE_inla_test <- read.csv("./data/INLA_size_standardized_data/reduced_sample_predicts/WE_inla_SampleSize_test.2021.09.22.csv", check.names = FALSE)

max.n <- filter(as.data.frame(table(WE_inla_prediction$Waterbody)), Freq >= 10)
max.n <- max(max.n$Freq)

## Plot all lakes
WE_test_plot <- ggplot(WE_inla_test, aes(x = n, y = `0.5quant`, colour = Waterbody, shape = Waterbody)) +
        geom_point(size = 4) +
        geom_line() +  
#        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +        
                
        scale_shape_manual(values = 1:12) +
        scale_x_continuous("Sample size", labels = as.character(1:max.n), breaks = 1:max.n) +
        
        labs(title = "Walleye INLA Prediction Test", y = "[As] ug/g prediciton (0.5 quantile)") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.title = element_text(size = 16)) + 
        theme(axis.text = element_text(size = 10))

#ggsave("./figures/INLA_SampleSize_Test/WE_AllLakes.2021.09.22.jpg", WE_test_plot, height = 8, width = 12, dpi = 300)

## Kenogamsisis
WE_keno_test_plot <- ggplot(subset(WE_inla_test, Waterbody == "Kenogamisis Lake" ), aes(x = n, y = `0.5quant`)) +
        geom_point(colour = "steelblue3", size = 4) +
        geom_line(colour = "steelblue1") +  
        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), color = "black", width = 0.2) +        
        
        scale_x_continuous("Sample size", labels = as.character(1:max.n), breaks = 1:max.n) +
        
        labs(title = "Walleye INLA Prediction Test; Kenogamisis Lake", y = "[As] ug/g prediciton (0.5 quantile)") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.title = element_text(size = 16)) + 
        theme(axis.text = element_text(size = 10))

#ggsave("./figures/INLA_SampleSize_Test/WE_Kenogamisis.2021.09.22.jpg", WE_keno_test_plot, height = 8, width = 12, dpi = 300)

## Long lake (Sudbury)
WE_long_test_plot <- ggplot(subset(WE_inla_test, Waterbody == "Long Lake (Sudbury)" ), aes(x = n, y = `0.5quant`)) +
        geom_point(colour = "firebrick3", size = 4) +
        geom_line(colour = "firebrick1") +  
        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), color = "black", width = 0.2) +        
        
        scale_x_continuous("Sample size", labels = as.character(1:32), breaks = 1:32) +
        
        labs(title = "Walleye INLA Prediction Test; Long Lake (Sudbury)", y = "[As] ug/g prediciton (0.5 quantile)") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.title = element_text(size = 16)) + 
        theme(axis.text = element_text(size = 10)) +
        coord_cartesian(y = c(0,0.4))

#ggsave("./figures/INLA_SampleSize_Test/WE_LongSudbury.2021.09.22.jpg", WE_long_test_plot, height = 8, width = 12, dpi = 300)

## Attawapiskat Lake
WE_atta_test_plot <- ggplot(subset(WE_inla_test, Waterbody == "Attawapiskat Lake" ), aes(x = n, y = `0.5quant`)) +
        geom_point(colour = "forestgreen", size = 4) +
        geom_line(colour = "green3") +  
        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), color = "black", width = 0.2) +        
        
        scale_x_continuous("Sample size", labels = as.character(1:28), breaks = 1:28) +
        
        labs(title = "Walleye INLA Prediction Test; Attawapiskat Lake", y = "[As] ug/g prediciton (0.5 quantile)") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.title = element_text(size = 16)) + 
        theme(axis.text = element_text(size = 10)) +
        coord_cartesian(y = c(0,0.4))

#ggsave("./figures/INLA_SampleSize_Test/WE_Attawapiskat.2021.09.22.jpg", WE_atta_test_plot, height = 8, width = 12, dpi = 300)


## Identify and plot lakes that deviate most from grand intercept (while still having relatively large sample size)
filter(INLA_WE_RWT_summary, Type == "Random Intercept - Waterbody" & Waterbody %in% test.lk) %>% arrange(desc(abs(`0.5quant`)))

## Bigwood Lake
(WE_bigwood_test_plot <- ggplot(subset(WE_inla_test, Waterbody == "Bigwood Lake" ), aes(x = n, y = `0.5quant`)) +
        geom_point(colour = "orchid4", size = 4) +
        geom_line(colour = "orchid2") +  
        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), color = "black", width = 0.2) +        
        
        scale_x_continuous("Sample size", labels = as.character(1:28), breaks = 1:28) +
        
        labs(title = "Walleye INLA Prediction Test; Bigwood Lake", y = "[As] ug/g prediciton (0.5 quantile)") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) +
        theme(axis.title = element_text(size = 16)) + 
        theme(axis.text = element_text(size = 10)) +
        coord_cartesian(y = c(0,0.4)))

#ggsave("./figures/INLA_SampleSize_Test/WE_Bigwood.2021.09.22.jpg", WE_bigwood_test_plot, height = 8, width = 12, dpi = 300)


## Lang Lake
(WE_lang_test_plot <- ggplot(subset(WE_inla_test, Waterbody == "Lang-Lake" ), aes(x = n, y = `0.5quant`)) +
                geom_point(colour = "salmon4", size = 4) +
                geom_line(colour = "salmon2") +  
                geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), color = "black", width = 0.2) +        
                
                scale_x_continuous("Sample size", labels = as.character(1:28), breaks = 1:28) +
                
                labs(title = "Walleye INLA Prediction Test; Lang Lake", y = "[As] ug/g prediciton (0.5 quantile)") +
                theme(plot.title = element_text(hjust = 0.5, size = 20)) +
                theme(axis.title = element_text(size = 16)) + 
                theme(axis.text = element_text(size = 10)) +
                coord_cartesian(y = c(0,0.4)))

#ggsave("./figures/INLA_SampleSize_Test/WE_Lang.2021.09.22.jpg", WE_lang_test_plot, height = 8, width = 12, dpi = 300)





















