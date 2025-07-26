%--------------------------------------------------------------------------
% This script generates Figure 2 for the manuscript.
% It analyzes human gloss perception data by:
% - Organizing participant responses and ground-truth values by group,
% - Computing inter-observer, intra-observer, and human-vs-ground-truth correlations,
% - Creating visualizations: scatter plots, histograms, and 3×3 tiled image grids
%   illustrating correspondence between human and ground-truth gloss settings.
% The script saves each figure panel in the 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

clearvars; close all; % cleaning

disp('Generating figure 2...')

%% Load Data
load(fullfile('data','onlineData'),'gt','data')
load(fullfile('data','fig_parameters'))

%% Organize responses and compute metrics
nGroups = 54;
nSubjects = length(data);

% Extract group number for each participant
groupN_summary = arrayfun(@(x) x.groupN, data);

% Collect response data for all participants
onlineData = zeros(size(data(1).response_Pellacini_c,1), size(data(1).response_Pellacini_c,2), nSubjects);
for i = 1:nSubjects
    onlineData(:,:,i) = data(i).response_Pellacini_c;
end

% Organize responses and ground truth by group
response_group = struct();
gt_group = struct();
for g = 1:nGroups
    idx = find(groupN_summary == g);
    groupResponses = zeros(size(data(1).response_Pellacini_c,1), numel(idx));
    for j = 1:numel(idx)
        groupResponses(:,j) = mean(data(idx(j)).response_Pellacini_c, 2);
    end
    response_group.(sprintf('group%d',g)) = groupResponses;
    gt_group.(sprintf('group%d',g)) = gt(idx(1)).Pellacini_c; % Representative ground truth for the group
end

% Calculate human vs. human and human vs. ground truth correlations
human_human_corr = [];
human_gt_corr = [];
for g = 1:nGroups
    groupData = response_group.(sprintf('group%d',g));
    gtData = gt_group.(sprintf('group%d',g));
    nSubj = size(groupData,2);
    for s = 1:nSubj
        obs1 = groupData(:,s);
        obs_rest = mean(groupData(:,setdiff(1:nSubj,s)),2);
        human_human_corr(end+1,1) = corr(obs1, obs_rest);
        human_gt_corr(end+1,1) = corr(obs1, gtData);
    end
end

% Correlation: ground truth vs. human (average response for each participant)
corrCoeff.gtvshuman = arrayfun(@(i) corr(gt(i).Pellacini_c, mean(data(i).response_Pellacini_c,2)), 1:nSubjects)';

% Correlation: human vs. human (across participants, last 12 images)
corrCoeff.acrossParticipant_all = [];
for i = 1:nSubjects
    for j = i+1:nSubjects
        obs1 = mean(data(i).response_Pellacini_c(75:end,:),2); % last 12 images
        obs2 = mean(data(j).response_Pellacini_c(75:end,:),2);
        corrCoeff.acrossParticipant_all(end+1,1) = corr(obs1, obs2);
    end
end

% Within-participant correlation (between the two columns in the response matrix)
corrCoeff.withinParticipant = arrayfun(@(i) corr(onlineData(:,1,i), onlineData(:,2,i)), 1:nSubjects)';

%% Scatter plot (fig2a)
load(fullfile('data','humanANDgt_Pellacini_c'))

fig = figure;
line([0 1],[0 1],'LineWidth',0.5,'Color',[0 0 0]);hold on
scatter(gt_Pellacini_c, human_Pellacini_c.mean,10,[0.5,0.5,0.5],'o','filled','LineWidth',0.5,'MarkerFaceAlpha',0.1,'MarkerEdgeColor','none');

ax = gca;
fig.Units           = 'centimeters';
fig.Position = [10,10,figp.twocolumn/4,figp.twocolumn/4];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

margin = 0.003;
xlim([0-margin 0.15+margin]);ylim([0-margin 0.15+margin])

xticks(linspace(0,0.1458,2))
yticks(linspace(0,0.1458,2))

ax.XTickLabel = {'0.00','0.15'};
ax.YTickLabel = {'0.00','0.15'};

xlabel('Ground-truth','FontWeight', 'Bold');ylabel('Observer setting','FontWeight', 'Bold');

ax.FontName = 'Arial';
ax.Color = [.97 .97 .97];
ax.FontSize = figp.fontsize;
ax.XColor = 'k';ax.YColor = 'k';
ax.TickDir = 'out';

ax.LineWidth = 0.5;
ax.Units = 'centimeters';
ax.Position = [1.0 0.9 3.2 3.2];
ticklengthcm(ax,0.1)
grid minor
box off
exportgraphics(fig,fullfile('figs','fig2a(scatter_humanvsgt).pdf'),'ContentType','vector')

%% Draw histograms (fig2b, c, d)
gray = [.3 .3 .3];

fig_histogram(corrCoeff.gtvshuman,calculateBins(corrCoeff.gtvshuman),'Corr. obs vs. gt','Frequency',gray,fullfile('figs','fig2b(hist_humanvsgt_correlation).pdf'))
fig_histogram(corrCoeff.withinParticipant,calculateBins(corrCoeff.withinParticipant),'Corr. within obs','Frequency',gray,fullfile('figs','fig2c(hist_within_participant_correlation).pdf'))
fig_histogram(human_human_corr,calculateBins(human_human_corr),'Corr. across obs','Frequency',gray,fullfile('figs','fig2d_1(hist_acrossParticipants_correlation).pdf'))
fig_histogram(corrCoeff.acrossParticipant_all,30,'Corr. across obs','Frequency',gray,fullfile('figs','fig2d_2(hist_humanvshuman_common12images).pdf'))

%% Make tiled images (fig2e)
rng(0); % for reproduciability
load(fullfile('data','humanANDgt_Pellacini_c'))

% Parameters
n_cell_rows = 3; % number of bins (human)
n_cell_cols = 3; % number of bins (gt)
n_img_per_cell = 8;
img_tile_rows = 2;
img_tile_cols = 4;

img_gap_px = 5;    % gap between images inside 2x4
panel_gap_px = 40;  % wider gap between 3x3 panels

% Define bin edges and labels
edges = linspace(0,0.1487,4);
labels = {'Low', 'Mid', 'High'};

gt_bins = discretize(gt_Pellacini_c, edges);
human_bins = discretize(human_Pellacini_c.mean, edges);

cell_images = cell(n_cell_rows, n_cell_cols);

for i = 1:n_cell_rows
    for j = 1:n_cell_cols
        idx = find(human_bins == i & gt_bins == j);
        if numel(idx) > n_img_per_cell
            idx = idx(randperm(numel(idx), n_img_per_cell));
        end
        img_list = cell(1, n_img_per_cell);
        for k = 1:n_img_per_cell
            if k <= numel(idx)
                img_path = fullfile('data','imgs','imgs3888_bg_png',sprintf('img%d.png', idx(k)));
                if exist(img_path, 'file')
                    img = imread(img_path);
                    if size(img,3) == 1, img = repmat(img,1,1,3); end
                    img_list{k} = img;
                else
                    img_list{k} = [];
                end
            else
                img_list{k} = [];
            end
        end
        img_list = fillMissingWithWhite(img_list);
        cell_images{i,j} = tileImagesWithGap(img_list, img_tile_rows, img_tile_cols, img_gap_px);
    end
end

% Determine panel (subpanel) size
cell_height = size(cell_images{1,1},1);
cell_width  = size(cell_images{1,1},2);

% Build panel-wide vertical and horizontal gaps
v_gap = uint8(255 * ones(cell_height, panel_gap_px, 3)); % vertical panel gap
h_gap = uint8(255 * ones(panel_gap_px, (cell_width+panel_gap_px)*n_cell_cols - panel_gap_px, 3)); % horizontal panel gap

% Compose rows with wide vertical gap between panels
row_blocks = cell(1, n_cell_rows);
for i = 1:n_cell_rows
    row_cells = cell(1, 2*n_cell_cols-1);
    for j = 1:n_cell_cols
        % Flip row index: high at top, low at bottom
        row_cells{2*j-1} = cell_images{n_cell_rows+1-i, j};
        if j < n_cell_cols
            row_cells{2*j} = v_gap;
        end
    end
    row_blocks{i} = cat(2, row_cells{:});
end

% Compose final image with wide horizontal gaps between panel-rows
final_blocks = cell(1, 2*n_cell_rows-1);
for i = 1:n_cell_rows
    final_blocks{2*i-1} = row_blocks{i};
    if i < n_cell_rows
        final_blocks{2*i} = h_gap;
    end
end

final_image = cat(1, final_blocks{:});

figure; imshow(final_image); axis image off
imwrite(final_image, fullfile('figs','fig2e(tiled_grid).png'));

disp('Done.')
close all

%% Functions
function out_list = fillMissingWithWhite(img_list)
    sz = [];
    for n = 1:numel(img_list)
        if ~isempty(img_list{n})
            sz = size(img_list{n});
            break;
        end
    end
    if isempty(sz), sz = [128 128 3]; end
    for n = 1:numel(img_list)
        if isempty(img_list{n})
            img_list{n} = uint8(255*ones(sz));
        end
    end
    out_list = img_list;
end

function tiled = tileImagesWithGap(img_list, nrows, ncols, gap_px)
    idx = 1;
    h = size(img_list{1},1); w = size(img_list{1},2);
    row_blocks = cell(1, nrows);
    for i = 1:nrows
        row_imgs = cell(1, 2*ncols-1);
        for j = 1:ncols
            row_imgs{2*j-1} = img_list{idx};
            idx = idx + 1;
            if j < ncols
                row_imgs{2*j} = uint8(255*ones(h, gap_px, 3));
            end
        end
        row_blocks{i} = cat(2, row_imgs{:});
    end
    for i = 2:nrows
        row_blocks{i} = cat(1, uint8(255*ones(gap_px, size(row_blocks{1},2), 3)), row_blocks{i});
    end
    tiled = cat(1, row_blocks{:});
end

%% Functions
function fig_histogram(data, nbins, x_label, y_label, color, filename)
    load(fullfile('data','fig_parameters'))
    fig = figure;
    data_median = median(data); hold on
    line([data_median,data_median],[-10000 10000],'Color','m','LineStyle',':','LineWidth',1)
    h = histogram(data, nbins, 'Normalization', 'count');
    h.EdgeColor = ones(3,1)*0.97;
    h.FaceColor = color;
    [~,peakid] = max(h.Values);
    peakval = mean([h.BinEdges(peakid),h.BinEdges(peakid+1)]);
    ax = gca;
    fig.Units = 'centimeters';
    fig.Position = [10,10,figp.twocolumn/5,figp.twocolumn/4];
    fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    ax.XLim = [0 1];
    ymax = max(h.Values);
    ax.YLim = [-2 ymax*1.1];
    xticks(0:0.25:1)
    yticks([0 ymax])
    t1 = text(data_median-0.07, ymax*1.15, num2str(round(data_median*100)/100), 'Color','m');
    t1.FontSize = figp.fontsize; t1.FontName = 'Arial';
    ax.XTickLabel = {'0.00','','0.50','','1.00'};
    ax.YTickLabel = {'0',num2str(round(ymax*100)/100)};
    xlabel(x_label,'FontWeight', 'Bold'); ylabel(y_label,'FontWeight', 'Bold');
    ax.FontName = 'Arial'; ax.Color = ones(3,1)*0.97;
    ax.FontSize = figp.fontsize; ax.XColor = 'k'; ax.YColor = 'k';
    ax.LineWidth = 0.5; ax.Units = 'centimeters';
    ax.Position = [1.05 0.9 2.1 3.2];
    ax.TickDir = 'out';
    ticklengthcm(ax,0.1)
    grid off; ax.YGrid = 'off'; box off
    exportgraphics(fig, filename, 'ContentType', 'vector')
end

function numBins = calculateBins(data)
    data = data(:);
    Q1 = quantile(data, 0.25); Q3 = quantile(data, 0.75); IQR = Q3 - Q1;
    binWidth = 2 * IQR * length(data)^(-1/3);
    numBins = ceil((max(data) - min(data)) / binWidth);
end
