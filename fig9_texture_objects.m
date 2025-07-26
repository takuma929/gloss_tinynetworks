%--------------------------------------------------------------------------
% This script generates Figure 9b and related panels for the manuscript.
% It visualizes gloss predictions for textured objects by:
% - Creating tiled montages of low- and high-contrast example images,
% - Loading and comparing predicted glossiness from single- and double-kernel models
%   (for both uniform and textured conditions),
% - Plotting scatter plots of model predictions (with mean markers and formatting),
% - Saving all figures in the 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clear workspace
clearvars; close all;

disp('Generating figure 9...')

load(fullfile('data', 'fig_parameters'))

%% Save tiled images (low vs high contrast examples)
indices = [200, 794, 1034, 1214, 1224, 1480];
img_folder_low  = fullfile('data', 'imgs', 'imgs_texture_nobg_lowcontrast');
img_folder_high = fullfile('data', 'imgs', 'imgs_texture_nobg_highcontrast');

w = 3; h = 2; % tile grid size: width x height
tile_size = 128;
gap = 5;

tiled_image_lowcontrast  = generateTiledImageFromFiles(indices, tile_size, gap, w, h, img_folder_low);
tiled_image_highcontrast = generateTiledImageFromFiles(indices, tile_size, gap, w, h, img_folder_high);

imwrite(tiled_image_lowcontrast, fullfile('figs', 'fig9b(imgs_lowcontrast).png'));
imwrite(tiled_image_highcontrast, fullfile('figs', 'fig9b(imgs_highcontrast).png'));

%% Scatter plot: predicted glossiness (uniform vs. textured)
cont_list.singlekernel = {'lowcontrast', 'highcontrast'};
cont_list.doublekernel = {'highcontrast'};

cmap.singlekernel = [90, 105, 124]/255;
cmap.doublekernel = [228, 108, 12]/255;

for model = {'singlekernel', 'doublekernel'}
    for cont = cont_list.(model{1})
        T = readtable(fullfile('data', 'pred_validation', ...
            ['prediction_textured_', model{1}, '_', cont{1}, '.csv']));
        
        fig = figure;hold on
        ax = gca;

        % Reference line y = x
        line([0,1], [0,1], 'LineWidth', 0.1, 'Color', [0.3 0.3 0.3]);

        % Scatter plot of predictions
        scatter(T.uniform, T.textured, 5, cmap.(model{1}), 'o', ...
            'filled', 'LineWidth', 0.1, 'MarkerFaceAlpha', 0.2);
        
        % Mean point
        scatter(mean(T.uniform), mean(T.textured), 70, 'bx', 'LineWidth', 1.5);

        % Figure formatting
        fig.Units = 'centimeters';
        fig.Position = [10, 10, figp.twocolumn/4, figp.twocolumn/4];
        fig.Color = 'w';
        fig.InvertHardcopy = 'off';
        
        % Axis settings
        xlim([0, 0.14]); ylim([0, 0.14]);
        xticks([0, 0.07, 0.14]);
        yticks([0, 0.07, 0.14]);
        xlabel('Predicted glossiness (uniform)');
        ylabel('Predicted glossiness (textured)');
        ax.XTickLabel = {'0.00', '0.07', '0.14'};
        ax.YTickLabel = {'0.00', '0.07', '0.14'};
        ax.FontName = 'Arial';
        ax.FontSize = figp.fontsize;
        ax.Color = [0.97, 0.97, 0.97];
        ax.XColor = 'k'; ax.YColor = 'k';
        ax.LineWidth = 0.5;
        ax.Units = 'centimeters';
        ax.Position = [0.9, 0.9, 3.3, 3.3];

        ticklengthcm(ax, 0.0)
        axis xy
        grid off
        box off

        % Save figure
        if strcmp(model{1},'singlekernel')
            filename = sprintf('fig9b(scatter_singlekernel_%s).pdf', cont{1});
        elseif strcmp(model{1},'doublekernel')
            filename = sprintf('fig9e(scatter_doublekernel_%s).pdf', cont{1});
        end
        
        exportgraphics(fig, fullfile('figs', filename), 'ContentType', 'vector')
    end
end

disp('Done.')
close all

%% Function: Generate tiled image from file indices
function tiled_image = generateTiledImageFromFiles(indices, tile_size, gap, w, h, img_folder)
% Generates a tiled image from a list of image indices.

    if numel(indices) ~= w * h
        error('Number of indices (%d) must match grid size w * h = %d', numel(indices), w * h);
    end

    tiles = cell(1, w * h);

    % Load and resize each image
    for i = 1:numel(indices)
        img_name = fullfile(img_folder, sprintf('img%d.png', indices(i)));
        if ~isfile(img_name)
            error('File not found: %s', img_name);
        end
        img = imresize(imread(img_name), [tile_size, tile_size]);
        tiles{i} = img;
    end

    % Determine output image size
    if size(tiles{1}, 3) == 1
        % Grayscale
        tiled_image = uint8(255 * ones(h * tile_size + (h - 1) * gap, ...
                                      w * tile_size + (w - 1) * gap));
    else
        % RGB
        tiled_image = uint8(255 * ones(h * tile_size + (h - 1) * gap, ...
                                      w * tile_size + (w - 1) * gap, 3));
    end

    % Assemble tiles into grid
    for row = 1:h
        for col = 1:w
            idx = (row - 1) * w + col;
            row_start = (row - 1) * (tile_size + gap) + 1;
            col_start = (col - 1) * (tile_size + gap) + 1;
            tiled_image(row_start:row_start+tile_size-1, ...
                        col_start:col_start+tile_size-1, :) = tiles{idx};
        end
    end
end
