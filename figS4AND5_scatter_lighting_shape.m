%--------------------------------------------------------------------------
% This script generates Supplementary Figures S2 and S3 for the manuscript.
% It analyzes human-vs-ground-truth gloss judgements by:
% - Plotting scatter plots of human vs. ground-truth gloss ratings for all 36
%   lighting environments and object shapes (with corresponding thumbnails),
% - Calculating and annotating correlation coefficients for each group,
% - Assembling ranked group images for both lighting and shape categories,
% - Saving all panels and final composite images in 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clear workspace
clearvars; close all;

disp('Generating figure S2 and S3...')

%% Load data
load(fullfile('data','fig_parameters'))
load(fullfile('data','humanANDgt_Pellacini_c'))
load(fullfile('data','onlineData'))
load(fullfile('data','allgroup_imgN'))

navy = [91 105 124]/255 * 0.8;

%% Count number of observers per image
groupN_summary = zeros(54,1);
for N = 1:length(data)
    onlineData(:,:,N) = data(N).response_Pellacini_c;
    groupN_summary(N) = data(N).groupN;
end

img3888_obsCount = zeros(3888,1);
for N = 1:295
    groupN = groupN_summary(N);
    img3888_obsCount(imgN_record(:,groupN)) = img3888_obsCount(imgN_record(:,groupN)) + 1;
end

%% Load lightprobe images
lightprobe_list = [1 4 5 7 8 9 10 13 16 17 18 20 21 22 23 24 25 28 29 30 31 32 33 34 36 37 39 40 41 42 43 45 49 50 54 56];
lightprobe_dir = fullfile('data','imgs','imgs_lightprobe');

for lightprobeN = 1:36
    img_path = fullfile(lightprobe_dir, ['LightProbe', num2str(lightprobe_list(lightprobeN)), '.png']);
    I.lightprobe(:,:,:,lightprobeN) = imresize(imread(img_path), [128 256]);
end

%% Load shape images
shape_dir = fullfile('data','imgs','imgs_obj_mesh');
for objN = 1:36
    img_path = fullfile(shape_dir, ['shape', num2str(objN), '.png']);
    I.shape(:,:,:,objN) = imresize(imread(img_path), [128 128]);
end

%% Scatter plots: Human vs GT grouped by lighting environments
table = readtable(fullfile('data','onlineExp_condition_summary.csv'));

for lightprobeN = 1:36
    id = find(table.lightprobeN == lightprobeN);
    fig = figure; hold on

    % Lightprobe thumbnail
    img = I.lightprobe(:,:,:,lightprobeN);
    image([-0.005 0.055], [0.123 0.153], flipud(img)*1.2)

    % Plot data
    errorbar(gt_Pellacini_c(id), human_Pellacini_c.mean(id), ...
        human_Pellacini_c.std(id) ./ sqrt(img3888_obsCount(id)), ...
        'LineStyle', 'none', 'LineWidth', 0.3, 'Color', [0 0 0], 'Capsize', 0);
    scatter(gt_Pellacini_c(id), human_Pellacini_c.mean(id), ...
        15, navy, 'filled', 'MarkerEdgeColor', [.9 .9 .9], 'LineWidth', 0.5);

    % Identity line
    line([-400 400], [-400 400], 'Color', [0 0 0], 'LineWidth', 0.5)
    
    % Correlation coefficient
    corrCoeff_gtvshuman.lightprobe(lightprobeN) = corr(gt_Pellacini_c(id), human_Pellacini_c.mean(id));
    t = text(0.123, 0.167, sprintf('%0.2f', round(corrCoeff_gtvshuman.lightprobe(lightprobeN)*100)/100));
    t.FontSize = figp.fontsize; t.FontName = 'Arial';
    
    % Aesthetics
    formatFigure(fig, figp.twocolumn / 6 * 0.95);
    configureAxes(gca);
    exportgraphics(fig, fullfile('figs','temp', ...
        ['figS2(scatter_humanvsgt_lightprobe', num2str(lightprobeN), ').png']), ...
        'ContentType', 'image', 'Resolution', 1200);
    close all
end

%% Save lightprobe group image
[~, rankN.lightprobe] = sort(corrCoeff_gtvshuman.lightprobe, 'descend');
out.lightprobe = assembleGroupImage('lightprobe', 2, rankN.lightprobe, 6, 6);
imwrite(out.lightprobe, fullfile('figs', 'figS2(scatter_humanvsgt_groupedbylighting).png'))

%% Scatter plots: Human vs GT grouped by shape
for objN = 1:36
    id = find(table.objN == objN);
    fig = figure; hold on

    % Shape thumbnail with transparency on white
    img = double(I.shape(:,:,:,objN))/255;
    alpha_channel = ~all(img >= 1.2, 3);
    hImg = image([-0.005 0.045], [0.105 0.155], flipud(img));
    set(hImg, 'AlphaData', alpha_channel);

    % Plot data
    errorbar(gt_Pellacini_c(id), human_Pellacini_c.mean(id), ...
        human_Pellacini_c.std(id) ./ sqrt(img3888_obsCount(id)), ...
        'LineStyle', 'none', 'LineWidth', 0.3, 'Color', [0 0 0], 'Capsize', 0);
    scatter(gt_Pellacini_c(id), human_Pellacini_c.mean(id), ...
        15, navy, 'filled', 'MarkerEdgeColor', [.9 .9 .9], 'LineWidth', 0.5);

    % Identity line
    line([-400 400], [-400 400], 'Color', [0 0 0], 'LineWidth', 0.5)

    % Correlation coefficient
    corrCoeff_gtvshuman.shape(objN) = corr(gt_Pellacini_c(id), human_Pellacini_c.mean(id));
    t = text(0.123, 0.167, sprintf('%0.2f', round(corrCoeff_gtvshuman.shape(objN)*100)/100));
    t.FontSize = figp.fontsize; t.FontName = 'Arial';

    formatFigure(fig, figp.twocolumn / 6 * 0.95);
    configureAxes(gca);
    exportgraphics(fig, fullfile('figs','temp', ...
        ['figS3(scatter_humanvsgt_shape', num2str(objN), ').png']), ...
        'ContentType', 'image', 'Resolution', 1200);
    close all
end

%% Save shape group image
[~, rankN.shape] = sort(corrCoeff_gtvshuman.shape, 'descend');
out.shape = assembleGroupImage('shape', 3, rankN.shape, 6, 6);
imwrite(out.shape, './figs/figS3(scatter_humanvsgt_groupedbyshape).png');

%% save rank information
save(fullfile('data','rankN_lightprobeANDshape'),'rankN')

disp('Done.')
close all

%% Functions
function configureAxes(ax)
    ax.XLim = [-0.008, 0.155];
    ax.YLim = [-0.008, 0.155];
    ax.XTick = [0, 0.15]; ax.YTick = [0, 0.15];
    ax.XTickLabel = {'',''}; ax.YTickLabel = {'',''};
    ax.Units = 'centimeters';
    ax.FontName = 'Arial';
    ax.FontSize = 7;
    ax.LineWidth = 0.5;
    ax.Color = [.95 .95 .95]/0.95;
    ax.XColor = 'k'; ax.YColor = 'k';
    ax.Position = [0.3 0.15 2.3 2.3];
    axis xy; grid off; box off
    ticklengthcm(ax, 0.0)
end

function formatFigure(fig, width_cm)
    fig.Units           = 'centimeters';
    fig.Position        = [10, 10, width_cm, width_cm];
    fig.Color           = 'w';
    fig.InvertHardcopy  = 'off';
    fig.PaperPosition   = [0, 10, 8.45, 8.45];
end

function outImg = assembleGroupImage(type, imgN, rankOrder, w, h)
    cnt = 0; clear figs
    for idx = rankOrder
        cnt = cnt + 1;
        img_path = fullfile('figs','temp', sprintf('figS%d(scatter_humanvsgt_%s%d).png', ...
            imgN, type, idx));
        img = imread(img_path);
        if cnt == 1
            imgsize = size(img);
        end
        figs(:,:,:,cnt) = padarray(imresize(img, imgsize(1:2)), [20,50], 255, 'both');
    end
    s = size(figs);
    outImg = reshape(permute(reshape(figs, [s(1),s(2),s(3),w,h]), [1 5 2 4 3]), [s(1)*h, s(2)*w, s(3)]);
end
