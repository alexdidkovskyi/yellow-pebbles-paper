---
title: "README"
output: html_document
---

# Yellow Pebbles: River Bed Load Transport Analysis

This repository contains the R code used in the analysis for the paper "Statistical Modelling of Sediment Transport Dynamics in Gravel-Bed Rivers Using RFID Tracer Observations" published in the Journal of Geophysical Research: Earth Surface.

## Overview

This codebase implements statistical approaches for analyzing sediment transport dynamics in gravel-bed rivers using Radio Frequency Identification (RFID) tracer data. The analysis includes:

1. Markov Chain analysis of tracer mobility states
2. Assessment of flood history effects on tracer mobility
3. XGBoost regression models for predicting tracer displacement and virtual velocity
4. Estimation and validation of bed load volumes

## Data

The analysis is based on a dataset of nearly 3,500 tracer displacement measurements collected during 27 sediment-mobilizing events in a Pre-Alpine reach in Italy (Caldone River). The dataset includes:

- Tracer characteristics (dimensions, weight)
- Tracer locations before and after flood events
- Flow rate and duration of mobilizing events
- Morphological units in the river reach

## Directory Structure

- `data/`: Data directory for input/output files
- `src/`: R code files organized by functionality
- `scripts/`: Analysis scripts that use the R functions

## Setup and Usage

### Requirements

The code requires R (version 4.0.0 or higher) and the following packages:

```R
# Core data processing and analysis
library(tidyverse)
library(magrittr)

# Statistical modeling
library(glmnet)
library(xgboost)
library(markovchain)

# Geospatial analysis
library(sp)
library(sf)
library(proj4)
library(raster)
library(geosphere)
library(igraph)
library(RTriangle)

# Visualization and output
library(ggplot2)
library(plotly)
library(openxlsx)
```

### Running the Analysis

1. Clone this repository
2. Run the script in the `scripts/` directory 

## Main Components

### Data Preprocessing

The preprocessing steps include:

- Loading and cleaning tracer locations and characteristics

- Calculating graph distances using triangulation

- Feature engineering and dimensionality reduction

### Markov Chain Analysis

Analysis of state transitions in:
- Two-state model (moved/not moved)

- Three-state model for displacement (long/typical/zero)

- Three-state model for virtual velocity (high/typical/zero)

### XGBoost Modeling

Regression modeling for:

- Typical displacement

- Long displacement

- Typical virtual velocity

- High virtual velocity

### Validation

Three-step validation approach:
1. Cross-validation with error metrics (RMSE, MAE, Wasserstein distance)
2. Prediction of movement during non-consecutive events
3. Comparison with sediment volumes in retention basin


## License

This code is provided for research purposes only. Please contact the authors for usage permissions.

## Contact

For questions or issues, please contact the authors of the manuscript