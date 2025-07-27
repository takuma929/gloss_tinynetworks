%--------------------------------------------------------------------------
% This script generates Figure 9b and related panels for the manuscript.
% It visualizes gloss predictions for textured objects by:
% - Creating tiled montages of low- and high-contrast example images,
% - Loading and comparing predicted glossiness from single- and double-kernel models
%   (for both uniform and textured conditions),
% - Plotting scatter plots of model predictions (with mean markers and formatting),
% - Conduct and report statistical test,
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

MAE_all = struct();

for model = {'singlekernel', 'doublekernel'}
    for cont = cont_list.(model{1})
        T = readtable(fullfile('data', 'pred_validation', ...
            ['prediction_textured_', model{1}, '_', cont{1}, '.csv']));
        
        % Compute MAE between uniform and textured
        MAE_uniform_vs_textured = abs(T.textured - T.uniform);
        MAE_all.(model{1}).(cont{1}) = MAE_uniform_vs_textured;
        
        % For info: mean, std, N
        mean_MAE = mean(MAE_uniform_vs_textured);
        std_MAE = std(MAE_uniform_vs_textured);
        fprintf('\nModel: %s, Contrast: %s\n', model{1}, cont{1});
        fprintf('MAE between uniform and textured: mean = %.4f, SD = %.4f, N = %d\n', ...
            mean_MAE, std_MAE, numel(MAE_uniform_vs_textured));
        
        % Normality
        [h_norm, p_norm] = lillietest(MAE_uniform_vs_textured);
        if h_norm
            fprintf('Normality test for MAE: NOT normal (Lilliefors p = %.4g)\n', p_norm);
        else
            fprintf('Normality test for MAE: normal (Lilliefors p = %.4g)\n', p_norm);
        end
        
        % ---- Existing plotting code ----
        fig = figure;hold on
        ax = gca;
        line([0,1], [0,1], 'LineWidth', 0.1, 'Color', [0.3 0.3 0.3]);
        scatter(T.uniform, T.textured, 5, cmap.(model{1}), 'o', ...
            'filled', 'LineWidth', 0.1, 'MarkerFaceAlpha', 0.2);
        scatter(mean(T.uniform), mean(T.textured), 70, 'bx', 'LineWidth', 1.5);
        fig.Units = 'centimeters';
        fig.Position = [10, 10, figp.twocolumn/4, figp.twocolumn/4];
        fig.Color = 'w';
        fig.InvertHardcopy = 'off';
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
        if strcmp(model{1},'singlekernel')
            filename = sprintf('fig9b(scatter_singlekernel_%s).pdf', cont{1});
        elseif strcmp(model{1},'doublekernel')
            filename = sprintf('fig9e(scatter_doublekernel_%s).pdf', cont{1});
        end
        exportgraphics(fig, fullfile('figs', filename), 'ContentType', 'vector')
    end
end

% --- Compare MAE between singlekernel and doublekernel (highcontrast) ---

MAE_single = MAE_all.singlekernel.highcontrast;
MAE_double = MAE_all.doublekernel.highcontrast;
diff_vals = MAE_double - MAE_single;

% Normality test for paired difference
[h_norm, p_norm] = lillietest(diff_vals);

fprintf('\n--- MAE Comparison between models (highcontrast) ---\n');
if h_norm
    fprintf('Normality test for difference: NOT normal (Lilliefors p = %.4g)\n', p_norm);
    % Wilcoxon signed-rank test (one-sided: doublekernel < singlekernel)
    [p_wil, h_wil, stats_wil] = signrank(MAE_double, MAE_single, 'tail', 'left');
    % Effect size (r = z / sqrt(n))
    z_val = stats_wil.zval;
    n = numel(diff_vals);
    r = z_val / sqrt(n);
    fprintf(['Wilcoxon signed-rank test (doublekernel < singlekernel): ' ...
        'z = %.3f, p = %.4g, effect size r = %.3f, n = %d\n'], ...
        z_val, p_wil, r, n);
else
    fprintf('Normality test for difference: normal (Lilliefors p = %.4g)\n', p_norm);
    [~, p, ci, stats] = ttest(MAE_double, MAE_single, 'Tail', 'left');
    [~, ~, ci_two, ~] = ttest(MAE_double, MAE_single, 'Tail', 'both');
    df = stats.df;
    tstat = stats.tstat;
    d = mean(diff_vals) / std(diff_vals);
    fprintf('Paired one-sided t-test (doublekernel < singlekernel): t(%d) = %.3f, p = %.4g, Cohen''s d = %.3f\n', ...
        df, tstat, p, d);
    fprintf('95%% CI for mean difference: [%.4f, %.4f]\n', ci_two(1), ci_two(2));
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
