%--------------------------------------------------------------------------
% This script generates Figure 7c/d for the manuscript.
% It shows model predictions against the Serrano et al. validation dataset by:
% - Loading predictions from one-layer and three-layer models,
% - Extracting lighting and shape categories for groupwise analysis,
% - Computing Pearson correlations between model predictions and human
%   response for each lighting and shape condition,
% - Plotting bar graphs of correlation coefficients for each model and grouping,
% - Annotating mean values and category labels,
% - Saving the figures in the 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clean workspace
clearvars; close all;

disp('Generating figure 7...')

%% Load data
load(fullfile('data', 'fig_parameters'))
load(fullfile('data', 'onlineData'))

%% Load Serrano dataset predictions
table = readtable(fullfile('data', 'pred_validation', 'prediction_Serrano_dataset_nobg.csv'));
imagename = table2cell(table(:,1));

% Extract lighting and shape categories
[list.Lighting, list.Shape] = getlist(imagename);

% Store predictions
resall.onelayer    = table2array(table(:,2));
resall.threelayer  = table2array(table(:,3));
resall.Serrano     = table2array(table(:,4));

%% Compute correlations grouped by lighting/shape
for type = {'Lighting', 'Shape'}
    cnt = 0;
    for target = list.(type{1})'
        cnt = cnt + 1;
        idx = contains(imagename, target{1});
        res.onelayer   = table2array(table(idx,2));
        res.threelayer = table2array(table(idx,3));
        res.Serrano    = table2array(table(idx,4));

        for area = {'onelayer', 'threelayer'}
            corrCoeff.(area{1}).(type{1})(:,cnt) = corr(res.(area{1}), res.Serrano, 'type', 'Pearson');
        end
    end
end

%% Prepare for plotting
cmap.onelayer   = [91 105 124]/255;
cmap.threelayer = [91 97 246]/255;

% Original data order for labeling
textList.threelayer.Shape    = 1:9;
textList.onelayer.Shape      = [1 4 2 9 6 3 8 5 7];
textList.threelayer.Lighting = 1:9;
textList.onelayer.Lighting   = [2 9 8 4 7 1 6 5 3];

% Align ordering between one- and three-layer models
[~, shapeIdx] = ismember(1:9, textList.onelayer.Shape);
textList.onelayer.Shape    = 1:9;
textList.threelayer.Shape  = textList.threelayer.Shape(shapeIdx);

[~, lightIdx] = ismember(1:9, textList.onelayer.Lighting);
textList.onelayer.Lighting   = 1:9;
textList.threelayer.Lighting = textList.threelayer.Lighting(lightIdx);

%% Plot correlation bar graphs
for area = {'onelayer', 'threelayer'}
    for type = {'Shape', 'Lighting'}
        name = list.(type{1});
        corr_list = corrCoeff.(area{1}).(type{1});
        corr_list = corr_list(~isnan(corr_list));
        [corr_sorted, ~] = sort(corr_list, 'descend');

        fig = figure; hold on
        b = bar(corr_sorted, 'EdgeColor', 'none', 'FaceColor', cmap.(area{1}));

        % Mean line
        meany = mean(corr_list);
        line([0, 100], [meany, meany], 'LineStyle', ':', 'LineWidth', 1, 'Color', [1 0.5 1])
        text(8.5, meany*1.08, num2str(round(meany*100)/100), 'FontSize', 7, 'FontName', 'Arial', 'Color', 'm')

        % Axis settings
        ax = gca;
        ax.FontName = 'Arial';
        ax.FontSize = figp.fontsize;
        ax.Color = ones(1,3) * 0.97;
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.LineWidth = 0.5;
        ax.Units = 'centimeters';

        xlim([0, length(corr_list)+1])
        ylim([-0.12, 0.82])
        ax.Position = [0.95, 0.37, 3.3, 3.9];
        yticks([-0.1, 0:0.25:1])
        ax.YTickLabel = {'-0.10','0.00','0.25','0.50','0.75'};

        xticks(0:0.25:1)
        ax.XTickLabel = {};

        xlabel(type{1})
        ylabel('Correlation coefficient')

        % Add numbers
        for n = 1:9
            num = textList.(area{1}).(type{1})(n);
            if corr_sorted(n) > 0
                text(n-0.26, -0.04, num2str(num), 'FontSize', figp.fontsize, 'FontName', figp.fontname)
            end
        end

        grid on; ax.XGrid = 'off'; box off
        fig.Units = 'centimeters';
        fig.Color = 'w';
        fig.InvertHardcopy = 'off';
        fig.Position = [10, 10, figp.twocolumn/4, figp.twocolumn/4];

        ticklengthcm(ax, 0.0)

        % Export figure
        filename = sprintf('fig7%s(%s_%s_Serrano).pdf', ...
            char('c' + strcmp(area{1}, 'threelayer')), area{1}, type{1});
        exportgraphics(fig, fullfile('figs', filename), 'ContentType', 'vector')
    end
end

disp('Done.')
close all

%% Function to extract lighting and shape labels
function [lightingEnvironments, shapes] = getlist(imageNames)
    lightingEnvironments = cell(size(imageNames));
    shapes = cell(size(imageNames));

    for i = 1:length(imageNames)
        tokens = regexp(imageNames{i}, '\[(.*?)\]', 'tokens');
        lightingEnvironments{i} = tokens{1}{1};
        shapes{i} = tokens{2}{1};
    end

    lightingEnvironments = unique(lightingEnvironments);
    shapes = unique(shapes);
end