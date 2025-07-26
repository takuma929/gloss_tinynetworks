%--------------------------------------------------------------------------
% This script generates Supplementary Figure S4 for the manuscript.
% It compares lab-based (offline) and online data by:
% - Computing human–ground-truth correlations per observer and grouping responses,
% - Loading and aligning data from both online and offline experiments,
% - Calculating within-observer, and group-level correlations,
% - Plotting scatter plots and histograms comparing offline (lab) and online results,
% - Saving all figures in the 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clear workspace
clearvars; close all;

disp('Generating figure S4...')

%% Load data
load(fullfile('data','onlineData'))
load(fullfile('data','imageStats_corrCoeff'))
load(fullfile('data','imgStats_multiRegression_corrCoeff'))
load(fullfile('data','fig_parameters'))

%% Compute GT vs Human correlation per observer
groupN_summary = zeros(54,1);
for N = 1:length(data)
    corrCoeff.gtvshuman(N) = corr(gt(N).Pellacini_c, mean(data(N).response_Pellacini_c, 2));
    onlineData(:,:,N) = data(N).response_Pellacini_c;
    groupN_summary(N) = data(N).groupN;
end

%% Organize response by group
for groupN = 1:54
    idx = find(groupN_summary == groupN);
    for i = 1:length(idx)
        response_group.(['group', num2str(groupN)])(:,i) = mean(data(idx(i)).response_Pellacini_c(1:84,:), 2);
        if i == 1
            gt_group.(['group', num2str(groupN)]) = gt(idx(i)).Pellacini_c(1:84);
        end
    end
end

%% Load and prepare offline data
clear offline_group corrCoeff
load(fullfile('data','offlineData.mat'))
offlineData = data; clear data
obsList = fieldnames(offlineData);

%% Online group means
online_group1 = mean(response_group.group1(1:72,:), 2);
online_group2 = mean(response_group.group2(1:72,:), 2);

%% Apply offline corrections (identity scaling for now)
for obs = obsList'
    data = offlineData.(obs{1});
    for imgN = 1:72
        for sessionN = 1:2
            % No actual scaling applied (scale = 0)
            data.Pellacini_c.group1(imgN, sessionN) = data.Pellacini_c.group1(imgN, sessionN);
            data.Pellacini_c.group2(imgN, sessionN) = data.Pellacini_c.group2(imgN, sessionN);
        end
    end
    offlineData.(obs{1}) = data;
end

%% Compute within-observer correlations (session 1 vs 2)
cnt = 0;
for obs = obsList'
    cnt = cnt + 1;
    data = offlineData.(obs{1});

    if cnt == 1
        [n_img, n_session] = size(data.Pellacini_c.group1);
        offline_group.group1 = zeros(n_img, n_session, numel(obsList));
        offline_group.group2 = zeros(n_img, n_session, numel(obsList));
    end

    offline_group.group1(:,:,cnt) = data.Pellacini_c.group1;
    offline_group.group2(:,:,cnt) = data.Pellacini_c.group2;

    all_resp = [data.Pellacini_c.group1; data.Pellacini_c.group2];
    corrCoeff.withinObs(cnt,1) = corr(all_resp(:,1), all_resp(:,2));
end

%% Correlation: average online vs offline
for group = {'group1','group2'}
    online_avg  = mean(response_group.(group{1})(1:72,:), 2);
    offline_avg = mean(squeeze(mean(offline_group.(group{1}), 2)), 2);
    corrCoeff.offlinevsonline_mean.(group{1}) = corr(online_avg, offline_avg);
end

%% Correlation: all combinations between online and offline observers
for group = {'group1','group2'}
    online  = response_group.(group{1})(1:72,:);
    offline = squeeze(mean(offline_group.(group{1}), 2));
    corr_list = [];

    for i = 1:size(offline,2)
        for j = 1:size(online,2)
            corr_list(end+1,1) = corr(online(:,j), offline(:,i));
        end
    end

    corrCoeff.offlinevsonline.(group{1}) = corr_list;
end

%% Scatter plots: online vs offline
navy = [91 105 124]/255 * 0.8;

for group = {'group1', 'group2'}
    offline_avg = mean(squeeze(mean(offline_group.(group{1}), 2)), 2);
    online_avg  = mean(response_group.(group{1})(1:72,:), 2);
    
    offline_se = std(squeeze(mean(offline_group.(group{1}), 2)), [], 2) / sqrt(size(offline_group.(group{1}), 3));
    online_se  = std(response_group.(group{1})(1:72,:), [], 2) / sqrt(size(response_group.(group{1}), 2));

    fig = figure; hold on
    line([-1 1], [-1 1], 'Color', 'k', 'LineWidth', 0.5)
    errorbar(online_avg, offline_avg, offline_se, offline_se, online_se, online_se, ...
        'LineStyle', 'none', 'Color', [0 0 0], 'LineWidth', 0.1, 'CapSize', 0)
    scatter(online_avg, offline_avg, 20, navy, 'filled', ...
        'MarkerEdgeColor', [.9 .9 .9], 'LineWidth', 0.5)

    xlabel('Observer setting - online', 'FontWeight', 'Bold')
    ylabel('Observer setting - lab', 'FontWeight', 'Bold')
    xlim([-0.01, 0.15]); ylim([-0.01, 0.15])
    xticks(0:0.05:0.15); yticks(0:0.05:0.15)

    ax = gca;
    ax.Units = 'centimeters';
    ax.FontName = 'Arial';
    ax.FontSize = figp.fontsize;
    ax.Color = [.97 .97 .97];
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 0.5;
    ax.Position = [0.9 0.9 3.3 3.3];
    ticklengthcm(ax, 0.0)
    axis square; grid off; box off

    fig.Units = 'centimeters';
    fig.Position = [10, 10, figp.twocolumn/4, figp.twocolumn/4];
    fig.Color = 'w'; fig.InvertHardcopy = 'off';

    exportgraphics(fig, fullfile('figs', ['figS4a(scatter_offlinevsonline_', group{1}, ').pdf']), 'ContentType', 'vector')
end

%% Histogram plots
close all
gray = [.3 .3 .3];

for group = {'group1','group2'}
    fig_histogram(corrCoeff.offlinevsonline.(group{1}), ...
        calculateBins(corrCoeff.offlinevsonline.(group{1})), ...
        'Corr. lab vs. online', 'Frequency', gray, ...
        fullfile('figs', ['figS4b(hist_offlinevsonline', group{1}, ').pdf']));
end

fig_histogram(corrCoeff.withinObs, 30, ...
    'Corr. within obs', 'Frequency', gray, ...
    fullfile('figs','figS4c(hist_offline_withinObs).pdf'))

disp('Done.')
close all

%% Functions
function fig_histogram(data, nbins, x_label, y_label, color, filename)
    load(fullfile('data', 'fig_parameters'))
    fig = figure; hold on

    data_median = median(data);
    h = histogram(data, nbins, 'Normalization', 'count', ...
        'EdgeColor', [0.97 0.97 0.97], 'FaceColor', color);

    line([data_median, data_median], [-1e5 1e5], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 1)
    ymax = max(h.Values);
    
    xlabel(x_label,'FontWeight', 'Bold');
    ylabel(y_label,'FontWeight', 'Bold');
    
    ax = gca;
    ax.Units = 'centimeters';
    ax.FontName = 'Arial';
    ax.FontSize = figp.fontsize;
    ax.Color = [.97 .97 .97];
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 0.5;
    ax.Position = [1.05 0.9 2.1 3.2];
    ax.TickDir = 'out';
    ax.XLim = [0 1];
    ax.YLim = [-ymax*0.02, ymax*1.1];
    ax.XTick = 0:0.25:1;
    ax.YTick = [0,ymax];
    ax.XTickLabel = {'0.00','','0.50','','1.00'};
    ax.YTickLabel = {'0', num2str(round(ymax*100)/100)};
    ticklengthcm(ax, 0.1)
    grid off; box off

    text(data_median - 0.07, ymax * 1.15, num2str(round(data_median*100)/100), ...
        'Color', 'm', 'FontSize', figp.fontsize, 'FontName', 'Arial')

    fig.Units = 'centimeters';
    fig.Position = [10, 10, figp.twocolumn/5, figp.twocolumn/4];
    fig.Color = 'w'; fig.InvertHardcopy = 'off';

    exportgraphics(fig, filename, 'ContentType', 'vector')
end

function numBins = calculateBins(data)
    data = data(:);
    binWidth = 2 * iqr(data) * length(data)^(-1/3);
    numBins = max(1, ceil((max(data) - min(data)) / binWidth));
end
