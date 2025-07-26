%--------------------------------------------------------------------------
% This script generates Figure 5b for the manuscript.
% It visualizes the t-SNE embeddings of image representations from a 
% three-layer network model by:
% - Loading image, mask, and human gloss rating data,
% - Applying masks to isolate objects and set backgrounds to white,
% - Displaying images or glossiness-colored scatter points in t-SNE space 
%   for different network layers and the pixel domain,
% - Saving each intermediate panel and composing a composite figure.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clear workspace
clearvars; close all;

disp('Generating figure 5b...')

%% Load data
load(fullfile('data', 'fig_parameters'))
load(fullfile('data', 'humanANDgt_Pellacini_c'))

% Normalize human responses to [0,1]
human3888 = human_Pellacini_c.mean;
human3888 = human3888 / max(human3888(:));

% Load t-SNE coordinates
load(fullfile('data', 'tSNE_threelayer_model.mat'))

% Create output directory if it doesn't exist
if ~exist(fullfile('figs', 'temp'), 'dir')
    mkdir(fullfile('figs', 'temp'))
end

%% Load and mask relevant images
for imgN = img_list
    imgPath = fullfile('data', 'imgs', 'imgs3888_nobg_png', ['img', num2str(imgN), '.png']);
    maskPath = fullfile('data', 'imgs', 'obj_mask', ['img', num2str(imgN), '_mask.png']);
    
    img = double(imresize(imread(imgPath), [128 128], 'nearest')) / 255;
    mask = double(imresize(imread(maskPath), [128 128], 'nearest')) / 255;
    
    img_masked = img .* repmat(mask, 1, 1, 3);
    black_pixels = all(img_masked == 0, 3);
    img_masked(repmat(black_pixels, [1, 1, 3])) = 1; % Set black pixels to white
    
    img3888(:, :, :, imgN) = img_masked;
end

%% Plot settings
img_s = 0.05;
cmapN = 1024;
cmap = brewermap(cmapN, '*PiYG');

I_all = [];

for plot_type = {'object', 'glossiness'}
    I = [];
    
    for area = {'pixelSpace', 'conv1_max', 'conv2_max', 'pooling_out'}
        Y = scaleXY(vals_tSNE.(area{1}));
        fig = figure; hold on
        
        if strcmp(plot_type{1}, 'object')
            cnt = 0;
            for imgN = img_list
                cnt = cnt + 1;
                img = flipud(img3888(:, :, :, imgN));
                
                % Create alpha channel for white background
                white_pixels = all(img >= 1, 3);
                alpha_channel = ones(size(img,1), size(img,2));
                alpha_channel(white_pixels) = 0;
                
                hImg = image(img, 'XData', [Y(cnt,1)-img_s, Y(cnt,1)+img_s], ...
                                'YData', [Y(cnt,2)-img_s, Y(cnt,2)+img_s]);
                set(hImg, 'AlphaData', alpha_channel);
            end
            
        elseif strcmp(plot_type{1}, 'glossiness')
            scatterColors = cmap(floor((human3888(img_list)+0.001)*(cmapN-1)), :);
            scatter(Y(:,1), Y(:,2), 30, scatterColors, 'o', ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3, 'MarkerFaceAlpha', 1);
        end
        
        % Figure aesthetics
        ax = gca;
        axis xy square
        fig.Units = 'centimeters';
        fig.Color = 'w';
        fig.InvertHardcopy = 'off';
        fig.PaperPosition = [0, 10, 8.45, 8.45];
        fig.Position = [10, 10, figp.twocolumn/4, figp.twocolumn/4];
        
        ax.XLim = [-0.1, 1.1];
        ax.YLim = [-0.1, 1.1];
        ax.XTick = [0, 0.15];
        ax.YTick = [0, 0.15];
        ax.XTickLabel = {};
        ax.YTickLabel = {};
        
        ax.FontName = 'Arial';
        ax.FontSize = figp.fontsize;
        ax.Color = [0.95, 0.95, 0.95];
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.LineWidth = 0.5;
        ax.Units = 'centimeters';
        ax.Position = [0.2, 0.22, 4.1, 4.1];
        grid off; box off;
        ticklengthcm(ax, 0.0)
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        
        % Save figure
        exportPath = fullfile('figs', 'temp', ['tSNE_', plot_type{1}, '_', area{1}, '.png']);
        exportgraphics(fig, exportPath, 'ContentType', 'image', 'Resolution', 1200)
        close all
        
        % Append to combined image
        I = [I, imresize(imread(exportPath), [512, 512])];
    end
    I_all = [I_all; I];
end

% Save final composite image
imwrite(I_all, fullfile('figs', 'fig5b(tSNE_plot).png'))


disp('Done.')
close all

%% Function to scale X and Y to [0,1]
function Y_scaled = scaleXY(Y)
    x = Y(:, 1);
    y = Y(:, 2);
    
    x_scaled = (x - min(x)) / (max(x) - min(x));
    y_scaled = (y - min(y)) / (max(y) - min(y));
    
    Y_scaled = [x_scaled, y_scaled];
end
