# Human Gloss Perception: Figure Generation and Data Processing
This repository contains the source code, example data, and scripts required to process behavioral and computational model data, and to generate all main and supplementary figures for the manuscript: Takuma Morimoto, Arash Akbarinia, Katherine Storrs, Jacob R. Cheeseman, Hannah E. Smithson, Karl R. Gegenfurtner and Roland W. Fleming, “Human gloss perception reproduced by tiny neural networks”. bioRxiv. https://doi.org/10.1101/2025.05.09.653112

################## (1) System Requirements ##################

Operating System:

Windows, macOS, Ubuntu (tested on macOS ver 15.3.2)

MATLAB Version:

MATLAB R2020a or newer (tested on R2020a, R2024b)

Toolboxes/Dependencies:
Image Processing Toolbox (tested on Version 11.1, Version 24.2)
Statistics and Machine Learning Toolbox (tesetd on Version 11.7, Version 24.2)
Curve Fitting Toolbox (tested on Version 3.5.11, Version 24.2)
Computer Vision Toolbox (tested on Version 9.2, Version 24.2)
Optimization Toolbox (tested on Version 8.5, Version 24.2)

No non-standard hardware required; standard desktop/laptop computer

################## (2) Installation Guide ##################

Instructions:

Download and unzip the zip file at https://github.com/takuma929/gloss_tinynetworks or clone the repository:

	git clone gloss_tinynetworks

Open MATLAB and set the root directory as your working directory.

Add all subfolders to the MATLAB path.

################## (3) Demo ##################

To generate all figures and process the data as in the manuscript:

Start MATLAB and navigate to the main folder.

Run the following command in the MATLAB Command Window:

    	main

All figures will be saved to the figs/ directory, and cleaned data files will be saved in the data/ directory.

Expected output:

PDF and PNG files for all figures (main + supplementary)

Cleaned MATLAB data files (e.g., onlineData.mat)

Typical install time:

Setup: Less than 1 minute (copy/unzip files and MATLAB path setup)

Data download: The full repository, including all source code and required example datasets, is > 6 GB.

Downloading the data files may take for a while, depending on your internet speed

Expected run time: ~ 5 – 10 minutes on a standard computer

################## (4) Instructions for Use ##################

To run individual figure scripts or customize analyses:

Each figure has its own script (e.g., fig3_model_comparison.m for Figure 3).

You can run any script individually after running process_onlinedata.m to prepare the data.

Reproducing results in the manuscript:

Running the main.m script will reproduce all results and figures reported in the paper.

Scripts are modular; users can adapt or extend analysis as needed.

################## (5) License ##################

This code is distributed under the MIT License.

################## (6) Open Source Repository ##################

GitHub: https://github.com/takuma929/gloss_tinynetworks

DOI: Data will be uploaded to Zenodo or a similar repository upon publication.

################## (7) Software Description and Documentation ##################

Key operations:
The software processes behavioral and computational model response data, applies exclusion criteria, computes summary statistics, and generates all figures in the manuscript.

Fundamental tasks:
Data cleaning, model evaluation, analysis, figure generation for publication.

Algorithms and approach:
Includes standard psychophysical data processing, correlation analysis, kernel fitting, t-SNE visualization, and comparison of model predictions with human perception.

Dependencies:
Standard MATLAB toolboxes only. Other dependencies are included in 'function' directory in the Github repository.