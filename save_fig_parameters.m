%--------------------------------------------------------------------------
% Figure Parameter Initialization Script
%
% This script initializes and saves commonly used figure parameters for 
% figures. It sets up figure dimensions, font sizes, and font styles, 
% and creates directories for storing figures and data. 
% These parameters can be loaded by other scripts to ensure consistent 
% formatting across all figures in the project.
%
% Outputs:
%   - Directory structure: 'data/' and 'figs/' folders created (if not present)
%   - Figure parameter file: 'data/fig_parameters.mat'
%
% Author: TM, 2025
%--------------------------------------------------------------------------

clearvars; close all;

disp('Saving figure parameters...')

% Create directories for saving figures and data if they do not already exist
if ~exist('figs', 'dir')
    mkdir('figs')
end

% Set figure width (in centimeters) for two-column and one-column formats
figp.twocolumn = 18.5;              % Standard two-column width
figp.onecolumn = figp.twocolumn/2;  % Standard one-column width

% Set font sizes for general text and axes
figp.fontsize      = 6;   % General font size
figp.fontsize_axis = 7;   % Axis label font size

% Set font name for all figures
figp.fontname = 'Arial';

% Save the figure parameter structure for use in other scripts
save(fullfile('data','fig_parameters'), 'figp')

disp('Done.')
