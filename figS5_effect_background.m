%--------------------------------------------------------------------------
% This script generates Supplementary Figure S5 for the manuscript.
% It compares model performance with and without background cues by:
% - Loading correlation coefficients for networks trained with/without image backgrounds,
% - Plotting paired observer results and means for each network type (one-layer and three-layer),
% - Formatting figures showing correlation to human judgments as a function of kernel number,
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clean workspace and load data
clearvars; close all;

disp('Generating figure S5...')

% Load figure parameters and correlation data
load(fullfile('data', 'fig_parameters'))
load(fullfile('data', 'corr_bg_nobg'), 'corr_nobg', 'corr_bg')

% Visualization parameters
gap  = 0.2;
navy = [0.2855, 0.3294, 0.3890];

%% Loop through network types and plot
for area = {'onelayer', 'threelayer'}
    % Extract correlation data
    bg_corr   = corr_bg.(area{1}).human;
    nobg_corr = corr_nobg.(area{1}).human;

    % Create figure
    fig = figure;
    hold on

    % Plot paired lines and scatter means
    for N = 1:size(bg_corr, 2)
        for n = 1:24
            line([N-gap, N+gap], [nobg_corr(n, N), bg_corr(n, N)], ...
                'LineWidth', 0.3, 'Color', [0.7 0.7 0.7]);
        end
        scatter(N-gap, mean(nobg_corr(:, N)), 25, navy, ...
            'o', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.2);
        scatter(N+gap, mean(bg_corr(:, N)), 25, navy, ...
            'd', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.2);
    end

    % Axes formatting
    ax = gca;
    ax.Units       = 'centimeters';
    ax.Position    = [1.0 0.9 2.4 3.2];
    ax.FontName    = 'Arial';
    ax.FontSize    = figp.fontsize;
    ax.Color       = [0.97 0.97 0.97];
    ax.XColor      = 'k';
    ax.YColor      = 'k';
    ax.TickDir     = 'out';
    ax.LineWidth   = 0.5;
    ax.XTick       = 1:4;
    ax.YTick       = 0.5:0.25:1;
    ax.YTickLabel  = {'0.50','0.75','1.00'};
    xlim([0.4 4.6]);
    ylim([0.5 1]);

    % Custom X-axis labels depending on area type
    if strcmp(area{1}, 'onelayer')
        ax.XTickLabel = {'1','2','4','9'};
    elseif strcmp(area{1}, 'threelayer')
        ax.XTickLabel = {'9','16','32','64'};
    end

    % Labels
    xlabel('kernel N', 'FontWeight', 'Bold');
    ylabel('Correlation to human', 'FontWeight', 'Bold');

    % Figure formatting
    fig.Units          = 'centimeters';
    fig.Position       = [10, 10, figp.twocolumn/5, figp.twocolumn/4];
    fig.Color          = 'w';
    fig.InvertHardcopy = 'off';

    ticklengthcm(ax, 0.1);
    grid off
    box off

    % Export to PDF
    exportgraphics(fig, fullfile('figs', ['figS5(bgvsnobg_', area{1}, ').pdf']), ...
        'ContentType', 'vector');
end


%% Statistical comparison: background vs. no-background
[h_two, p_two, ci_two, stats_two] = ...
    ttest(mean(corr_bg.threelayer.human, 2), mean(corr_nobg.threelayer.human, 2), 'tail', 'both');

[h_one, p_one, ci_one, stats_one] = ...
    ttest(mean(corr_bg.onelayer.human, 2), mean(corr_nobg.onelayer.human, 2), 'tail', 'both');

% Display results
disp('--- Three-layer stats ---');
disp(stats_two)
disp(['p = ', num2str(p_two)])

disp('--- One-layer stats ---');
disp(stats_one)
disp(['p = ', num2str(p_one)])

disp('Done.')
close all
