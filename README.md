# Bayesian Spatio-Temporal Dendritic River Modeling
This is the code for my master's project for the MS in Statistics at BYU. I tried to do well commenting everything, below I'll go into detail about what each file does.

## How to use

Originally, this data was formatted for the old SSN package that was depricated around R-version 4.3.3. The needed wieghts and distances that the SSN package is usually needed to calculate are saved as R objects so as to not need to learn how to install the old package by finding the old package and figuring out how to get R to go back to the old version. The data that I used for the project is not included, as I don't hold the rights to it. However, there is code to create simulated data over the river network that was used to explore.

The file structure as it now stands should be what you need to run everything. Note the subdirectories that contain weights needed for the tail-up correlation calculations.



### Core Files

The following are the key files in the `GitHub_folder` directory: The first contains the main code.

The main folder to use and explore the model is
- [`bayes_model_no23_noSSN.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/bayes_model_no23_noSSN.R):
     This is the file that can be used without the SSN package. The data are not included as previously mentioned, but later down the code one can simulate data using the stream network used in the project. The rest of the files are supplementary to support this one and don't need to be directly opened-though the functions.R file has helpful comments about the functions used within this file. Requires R (Version >= 4.3.3)

How to use without data: (For those just wanting to simulate river network data)
* Start from the top of the code, and comment out the sections that require the code (Should be marked with a comment or include any sort of reference to data). Skip over any lines that run in data. 
* Follow the comments and refer to functions.R for any unfamiliar functions used.
* Refer to paper about any math used in the functions. 

- [`bayes_model_no23.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/bayes_model_no23.R):
     This file is the original used for the analysis, it requires the SSN package-not usable without it. (Needs R version 4.3.3 or at minimum slightly older versions))

- [`functions.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/functions.R): Contains helper functions used across the project.
- [`get_weights.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/get_weights.R): Script to calculate weights, distances, and other helpful matrices for the model.
- [`slicy_dicey.R`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/blob/main/GitHub_folder/slicy_dicey.R): Code for the Slice Sampler Used.

### Subdirectories

- [`dist_and_weights_mats`](https://github.com/thomas-wit/bayesian_spatio_temp_dendritic_river_mod/tree/main/GitHub_folder/dist_and_weights_mats): Contains distance and weight matrices needed for calculations.


