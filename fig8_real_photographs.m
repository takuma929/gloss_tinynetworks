%--------------------------------------------------------------------------
% This script generates Figure 8 for the manuscript.
% It evaluates model classification performance on real-world photographs by:
%   - Loading human-provided binary gloss labels and model predictions,
%   - Determining optimal decision thresholds for both one-layer and three-layer models,
%   - Computing accuracy, d-prime, and threshold distance for each model,
%   - Plotting a scatter comparison of model predictions (Figure 8a),
%   - Generating image grids for correctly classified glossy and matte examples (Figures 8b, 8c),
%   - Reporting summary statistics to the console.
% Outputs: Figures in 'figs' directory and performance metrics for inclusion in the manuscript.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

% Script for evaluating responses of one-layer model and three-layer model
clearvars; close all;

disp('Generating figure 8...')

rng(0); % freeze seed for reprodiciability

% Load figure parameters (font sizes, colors, etc.)
load(fullfile('data','fig_parameters'))

% Load dataset containing observer responses and model predictions for real-world photographs 
T = readtable('./data/pred_validation/prediction_photographs.csv');

% Extract clean image numbers (remove .png extension)
numbers = regexprep(T.image_name, '\.png', '');

%% Classification and performance evaluation for each model/model
% Define models to compare: one-layer (single-kernel), three-layers
models = {'onelayer', 'threelayers'};
labels = T.proportion_glossy; % Binary labels (0: matte, 1:glossy)

for a = 1:length(models)
    model = models{a};
    predictions.(model) = T.(['glossiness_' model]); % Model output
    
    % Find best decision threshold to maximize accuracy
    [bestThresh.(model), ~] = find_best_threshold(predictions.(model), labels);

    % Convert continuous predictions to binary class using best threshold
    pred_classes.(model) = predictions.(model) > bestThresh.(model);

    % Compute classification accuracy
    accuracy.(model) = mean(pred_classes.(model) == labels) * 100;

    % --- Signal detection theory (d-prime) calculation ---
    % "1" is glossy, "0" is matte
    actual_positive = labels == 1;
    actual_negative = labels == 0;
    predicted_positive = pred_classes.(model) == 1;

    hit_rate = sum(predicted_positive & actual_positive) / sum(actual_positive);
    fa_rate  = sum(predicted_positive & actual_negative) / sum(actual_negative);

    % Correction for extreme values (avoid infinite d-prime)
    n_pos = sum(actual_positive); n_neg = sum(actual_negative);
    hit_rate = max(min(hit_rate, 1 - 1/(2*n_pos)), 1/(2*n_pos));
    fa_rate  = max(min(fa_rate , 1 - 1/(2*n_neg)), 1/(2*n_neg));

    dprime.(model) = norminv(hit_rate) - norminv(fa_rate); % Sensitivity index
end

%% Display image tile grids for correctly/incorrectly classified images (onelayer model)

% Find image indices for each classification outcome
correct_glossy_idx   = find((pred_classes.onelayer == 1) & (labels == 1));
correct_matte_idx    = find((pred_classes.onelayer == 0) & (labels == 0));
incorrect_glossy_idx = find((pred_classes.onelayer == 1) & (labels == 0));
incorrect_matte_idx  = find((pred_classes.onelayer == 0) & (labels == 1));

% Shuffle image order within each group and create grid images for visualization
N = 24; % Number of images per grid

% Correct glossy predictions
shuffled_names = T.image_name(correct_glossy_idx);
shuffled_names = shuffled_names(randperm(length(shuffled_names)));
makeTileImage(shuffled_names, 'fig8b(example_glossy_imgs).png', N);

% Correct matte predictions
shuffled_names = T.image_name(correct_matte_idx);
shuffled_names = shuffled_names(randperm(length(shuffled_names)));
makeTileImage(shuffled_names, 'fig8c(example_matte_imgs).png', N);

%% Generate scatter plot comparing model predictions
fig = figure; ax = gca; hold on;

x = predictions.onelayer;  % Single-kernel model predictions
y = predictions.threelayers;  % Three-layer model predictions

% Label groups for plotting
matte_idx  = labels == 0;
glossy_idx = labels == 1;

% Scatter: blue for matte, magenta for glossy (can adjust as needed)
scatter(x(matte_idx),  y(matte_idx),  80, [18,111,16]/255,  'o', 'filled', 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.7, 'LineWidth', 0.5);
scatter(x(glossy_idx), y(glossy_idx), 80, [208,61,139]/255, 'o', 'filled', 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.7, 'LineWidth', 0.5);

% Annotate each point with image index for reference
% comment in if you want to show image numbers
for idx = find(matte_idx)'
    %t = text(x(idx)-0.001, y(idx), numbers{idx}); t.FontSize = 4; t.Color = 'w';
end
for idx = find(glossy_idx)'
    %t = text(x(idx)-0.001, y(idx), numbers{idx}); t.FontSize = 4; t.Color = 'w';
end

% Decision boundary lines (vertical/horizontal) for each model
xline(bestThresh.onelayer, ':k', 'LineWidth', 1);
yline(bestThresh.threelayers, ':k', 'LineWidth', 1);

% Figure formatting (units, size, axis, labels)
fig.Units = 'centimeters';
fig.Position = [10, 10, figp.twocolumn/2, figp.twocolumn/2];
fig.Color = 'w'; fig.InvertHardcopy = 'off';

xlim([0,0.14]); ylim([-0.1,0.6]);
xticks([0 0.07,0.14]); yticks([0,0.3,0.6]);
xlabel('Prediction (single-kernel model)');
ylabel('Prediction (three-layer model)');
ax.XTickLabel = {'0.00','0.07','0.14'};
ax.YTickLabel = {'0.00','0.30','0.60'};
ax.FontName = 'Arial'; ax.Color = [.97 .97 .97]; ax.FontSize = figp.fontsize;
ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 0.5; ax.Units = 'centimeters';
ax.Position = [1.0 0.9 7.6 7.6];
ticklengthcm(ax, 0.0); % Standardize tick mark length (helper below)
axis xy; grid off; box off;

% Save vector figure
exportgraphics(fig, fullfile('figs', 'fig8a(real_photographs_prediction).pdf'), 'ContentType', 'vector');

%% Calculate and print summary statistics for each model (accuracy, d-prime, threshold distance)
for a = 1:length(models)
    name = models{a};
    preds = predictions.(name);
    pred_class = pred_classes.(name);
    thresh = bestThresh.(name);
    label = labels;

    % Identify misclassified indices and compute their distance from threshold
    misclassified_idx = find(pred_class ~= label);
    dists = abs(preds(misclassified_idx) - thresh);
    mean_dist = mean(dists);

    % Print performance metrics to console
    fprintf('Model: %s\n', name);
    fprintf('  Accuracy: %.2f%%\n', accuracy.(name));
    fprintf('  d-prime: %.2f\n', dprime.(name));
    fprintf('  Mean distance to threshold (misclassified images): %.4f\n\n', mean_dist);
end

disp('Done.')
close all

%% ------------------------- Functions ----------------------------
function [bestThresh, bestAcc] = find_best_threshold(preds, labels)
% Returns the threshold maximizing accuracy for binary classification.
    thresholds = linspace(min(preds), max(preds), 500);
    accs = arrayfun(@(t) mean((preds > t) == labels), thresholds);
    [bestAcc, idx] = max(accs);
    bestThresh = thresholds(idx);
end

function makeTileImage(imageList, outputFilename, N)
% Creates a grid image of N images from imageList and saves to outputFilename.
    gap = 5; imagesPerRow = 6;
    totalImgs = min(numel(imageList), N);
    tileImgs = cell(1, N);

    for i = 1:totalImgs
        img = imread(fullfile('./data/imgs/realworld_photographs/', imageList{i}));
        if ndims(img) == 2, img = repmat(img, 1, 1, 3); end
        img = cropToRatio(img, 768, 512); img = imresize(img, [256, 384]);
        tileImgs{i} = img;
    end

    % Fill unused slots with white
    for i = totalImgs+1:N
        tileImgs{i} = uint8(255*ones(256,384,3));
    end

    % Layout: rows/columns
    rows = ceil(N / imagesPerRow);
    tileHeight = size(tileImgs{1},1); tileWidth = size(tileImgs{1},2);
    tileCanvas = uint8(255*ones(rows*(tileHeight+gap)+gap, imagesPerRow*(tileWidth+gap)+gap, 3));

    for idx = 1:N
        r = ceil(idx / imagesPerRow);
        c = mod(idx-1, imagesPerRow) + 1;
        row_start = gap + (r-1)*(tileHeight+gap) + 1;
        col_start = gap + (c-1)*(tileWidth+gap) + 1;
        tileCanvas(row_start:row_start+tileHeight-1, col_start:col_start+tileWidth-1, :) = tileImgs{idx};
    end
    imwrite(tileCanvas, fullfile('./figs', outputFilename));
end

function imgOut = cropToRatio(img, targetWidth, targetHeight)
% Crops the input image to the specified aspect ratio, centered
    [h, w, ~] = size(img);
    targetRatio = targetWidth / targetHeight; imgRatio = w / h;
    if imgRatio > targetRatio
        newWidth = round(targetRatio * h);
        startX = round((w - newWidth) / 2);
        imgOut = imcrop(img, [startX, 1, newWidth-1, h-1]);
    else
        newHeight = round(w / targetRatio);
        startY = round((h - newHeight) / 2);
        imgOut = imcrop(img, [1, startY, w-1, newHeight-1]);
    end
end

function [rgb, M] = XYZToSRGBPrimary(XYZ)
    M = [3.2410 -1.5374 -0.4986 ; -0.9692 1.8760 0.0416 ; 0.0556 -0.2040 1.0570];
    if (~isempty(XYZ)), rgb = M*XYZ; else, rgb = []; end
end