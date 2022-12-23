# General terminology

## Model types

1. No item effects (Ind). Models are fitted with participant-effects but not item-effects. In SDT models, criteria c (1:Number_Thresholds) for individuals are sampled from a parent distribution with N(mu_cr_c,sigma_cr_c)
2. No Item effects, Random Crit (Indgam). Models are fitted with participant=effects but not item-effects. In SDT models, criteria c (1:Number_Thresholds) for individuals are implemented as participant-specific deviations for mu_cr_c (i.e., part of the random effects structure for distributional parameters)
3. Item effects, Random Crit (Item).  Models are fitted with participant-effects and item-effects. In SDT models, criteria c (1:Number_Thresholds) for individuals are implemented as participant-specific deviations for mu_cr_c (i.e., part of the random effects structure for distributional parameters)

## Model comparison/generalizability

1. Deviance (Full/Deviance). Model fitted to full data set, comparison on basis of in-sample deviamce
2. K(=10)-Fold crossvalidation (KFCV). Model fitted to 90% of data. 10% held-out for each participant, for ~10% held-out. Comparison on basis of out-of-sample deviance (repeated 10 times, summed across folds)
3. Leave-Participant-Out (LOP). Model fitted to 90% of data. 10% of participants held out. Comparison on basis of out-of-sample deviance (repeated 10 times, summed across folds) where predictions for individual-level deviation from grand mean for parameters are done on the basis of group-level parameters.

## Models

- Equal-variance Gaussian signal detection model (EVSDT)
- Equal-variance Gumbel signal detection model (Gumbel)
- Unqual-variance Gaussian signal detection model (UVSDT)
- Dual-process signal detection model (DPSDT)
- Dual-process signal detection model with response mapping (DPRMSDT/DPRM2SDT)
- Mixture signal detection model, strength of unattended items = new (M0SDT)
- Mixture signal detection model (MASDT)
- 2 High Threshold model (2HTM)

## Available data

- all processed datasets we collected for the meta-analysis is available in RecognitionData.zip, with the Shiny-App under https://nklange.shinyapps.io/RecognitionData/ allowing some overview over individual data sets
- all fits to the full data sets (Full/Deviance), i.e., MCMC samples fit object
- all predictions/deviance following Full, KFCV and LOP fits (largely available in greater detail for exploration in shiny apps)
- Not available: MCMC fit objects for the vast majority of KFCV and LOP fits

# Data used in the meta-analysis

Datafiles are contained in ProcessedData folder. These data sets include mark-ups for Folds for KFCV and LOP and constitute the data as we used them to fit models to. In some cases, these data sets are simplified or extracted from larger data sets. The full data sets prior the processing for our meta-analysis specifically, are available in Data.zip as .csv files

# Fitting data

The itemInd-column in the data set is used to determine whether item effects are fitted or not and how the data and models are prepared. 

Where itemInd = "Ind": no item effects  
where no item effects but random crit itemInd = "Indgam"  
where itemInd = "Item": item effects  

- MakeStanDataModel.R: prepare stan-friendly data from .rds data files for all models
- FitModel.R: fit models with stan

Resulting fits are in Fits folder.

# Make predictions

Different files to make predictions for different model types. Generally same functions with some variations to account for model variations

- LOPTheta.R: functions for models with item effects
- LOPTheta_Ind.R: functions for models without item effects
- LOPTheta_Ingam.R: functions for models without items effects + random criteria

- PredictDeviance.R: flow to predict Deviance (full fits) (as well as posterior p-values and LOOIC)
- PredictKFCV.E: flow to predict and out-of-sample K-Fold deviance
- PredictLOP.R: flow to predict out-of-sample Leave-Participant-Out deviance

Resulting files are in Predictions folder

# Aggregate predictions for visualizations

AggregateInd.R  
AggregateIndgam.R  
AggregateItem.R  

Resulting files are in ProcessedPredictions folder / as aggregated .rds files in Predictions

# Make visualizations 

MakeVisualizations.R: visualizations for manuscript for waffle plots and some distributional plots

# Extract Parameters

Extract parameters from fit files for visualiations in shiny apps. Resulting files in ParameterSummary folder

ExtractParameters.R  
ExtractParameters_Item.R  

# Shiny apps

## Data

distribution of all collected and tidied data (92 data sets): https://nklange.shinyapps.io/RecognitionData/ (Folder ShinyApp_Data)

## Fitted models with no item effects (79 data sets)

tbc

## Fitted models for all model types (10 data sets)

tbc



