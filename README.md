Human Gloss Perception and Tiny Neural Networks: Figure Generation and Data Processing
This repository contains the source code, example data, and scripts required to process behavioral and computational model data, and to generate all main and supplementary figures for the manuscript:

Takuma Morimoto, Arash Akbarinia, Katherine Storrs, Jacob R. Cheeseman, Hannah E. Smithson, Karl R. Gegenfurtner, and Roland W. Fleming
“Human gloss perception reproduced by tiny neural networks”. bioRxiv.
https://doi.org/10.1101/2025.05.09.653112

1. System Requirements
Operating System:
Windows, macOS, or Ubuntu (tested on macOS ver 15.3.2)

MATLAB Version:
R2020a or newer (tested on R2020a, R2024b)

Required MATLAB Toolboxes:

Image Processing Toolbox (tested on v11.1, v24.2)

Statistics and Machine Learning Toolbox (tested on v11.7, v24.2)

Curve Fitting Toolbox (tested on v3.5.11, v24.2)

Computer Vision Toolbox (tested on v9.2, v24.2)

Optimization Toolbox (tested on v8.5, v24.2)

Hardware:
Standard desktop or laptop computer (no special requirements)

2. Installation Guide
Download the repository:

Download ZIP or clone via terminal:

bash
Copy
Edit
git clone https://github.com/takuma929/gloss_tinynetworks.git
Open MATLAB, set the repository root directory as your working directory.

Add all subfolders to the MATLAB path.

3. Demo: Generate All Figures
To generate all figures and process the data as in the manuscript:

Start MATLAB and navigate to the main project folder.

Run the following command in the MATLAB Command Window:

matlab
Copy
Edit
main
All figures will be saved to the figs/ directory, and cleaned data files will be saved in the data/ directory.

Expected Output:

PDF and PNG files for all main and supplementary figures

Cleaned MATLAB data files (e.g., onlineData.mat)

Typical Install & Run Times:

Setup: < 1 minute (copy/unzip files, add to path)

Data download: Repository (with datasets) is > 6 GB; download time depends on your internet connection

Expected run time: ~5–10 minutes on a standard computer

4. Instructions for Use
Run Individual Figures or Custom Analyses:

Each figure has its own script (e.g., fig3_model_comparison.m for Figure 3)

Prepare the data first by running process_onlinedata.m

Scripts are modular and can be adapted or extended for your analysis

Full Reproducibility:

Running main.m will reproduce all results and figures as reported in the manuscript

5. License
Distributed under the MIT License.

6. Open Source Repository
GitHub: https://github.com/takuma929/gloss_tinynetworks

DOI: Data will be uploaded to Zenodo (or similar) upon publication.

7. Software Description & Documentation
Key Operations:
Processes behavioral and computational model response data, applies exclusion criteria, computes summary statistics, and generates all figures.

Fundamental Tasks:
Data cleaning, model evaluation, analysis, and publication-ready figure generation.

Algorithms & Approach:
Standard psychophysical data processing, correlation analysis, kernel fitting, t-SNE visualization, and comparison of model predictions with human perception.

Dependencies:
Standard MATLAB toolboxes only. Other code dependencies are included in the function directory.