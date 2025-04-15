# Bayesian Spatio-Temporal Dendritic River Modeling

This repository contains the code and resources for Bayesian spatio-temporal modeling of dendritic river networks. The project is implemented in R and focuses on analyzing and predicting spatial and temporal patterns in river networks using Bayesian statistical methods.

## Features

- **Bayesian modeling**: Implements probabilistic models for spatio-temporal data.
- **River networks**: Handles dendritic (tree-like) river network structures.
- **Spatio-temporal analysis**: Supports analysis and prediction over both space and time.

## Repository Structure

### Core Files

The following are the key files in the `GitHub_folder` directory:

- [`bayes_model_no23.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/bayes_model_no23.R): A script for Bayesian modeling (version without need for the SSN package).
- [`bayes_model_no23_noSSN.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/bayes_model_no23_noSSN.R): Another version of the Bayesian modeling script.
- [`functions.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/functions.R): Contains helper functions used across the project.
- [`get_weights.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/get_weights.R): Script to calculate weights for the model.
- [`slicy_dicey.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/slicy_dicey.R): Code for the Slice Sampler Used.

### Subdirectories

- [`dist_and_weights_mats`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/tree/main/GitHub_folder/dist_and_weights_mats): Contains distance and weight matrices.

## Requirements

To run the scripts in this repository, you need the following:

- **R** (version 4.3.3 or higher)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod.git
   cd bayesian_spatio_temp_dendritic_river_mod/GitHub_folder
