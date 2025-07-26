%--------------------------------------------------------------------------
% main.m — Master Script for Data Processing and Figure Generation
%
% This script serves as the main pipeline for the project. It executes all
% necessary data processing and figure generation steps for the manuscript.
%
% Functionality:
%   - Processes raw experimental and model data
%   - Applies exclusion criteria and data cleaning routines
%   - Runs all scripts needed to reproduce main and supplementary figures
%
% Workflow:
%   1. Clean workspace and close figures
%   2. Process raw online experimental data and save cleaned data structures
%   3. Sequentially generate all manuscript figures (main and supplementary)
%   4. Each figure script is modular and can be run independently if needed
%
% Output:
%   - Cleaned data (e.g., onlineData.mat)
%   - Manuscript figures (PDF, PNG) in the 'figs' directory
%
% Usage:
%   Simply run this script from the project root folder.
%
% Note:
%   - Ensure that this repository is present in the MATLAB path
%   - Some figure scripts may take several minutes to complete depending on
%     computational resources.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

clearvars;close all;clc %% cleaning

tic

save_fig_parameters % save figure parameters (e.g. fontsize) in "data" directory.

process_onlinedata % generate "onlineData.mat" in "data" directory.

fig2_online_exp_results % generate figure 2

fig3_model_comparison % generate figure 3

fig4b_24kernels % generate figure 4b

fig4c_4d_9c_analyze_kernel % generate figure 4c, 4d and 9c

fig4e_1kernel_fittingResults  % generate figure 4e

fit_kernel_gauss_ridges % perform fitting for 24 kernels from one-layer models

fig4f_24kernel_fittingResults % generate figure 4f

fig5b_tSNEplot % generate figure 5b

fig6_highlight_manipulation % generate figure 6

fig7_Serrano_dataset % generate figure 7

fig8_real_photographs % generate figure 8

fig9_texture_objects % generate figure 9

figS1_lighting_colordistribution % generate figure S1 (color distribution)

figS1_make_thumbnails % generate figure S1 (thumbnail)

figS2AND3_scatter_lighting_shape % generate figure S2 and S3

figS4_offlinevsonline % generate figure S4

figS5_effect_background % generate figure S5

toc