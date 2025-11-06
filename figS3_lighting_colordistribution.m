%--------------------------------------------------------------------------
% This script generates Supplementary Figure S1c for the manuscript.
% It visualizes the chromatic (CIELAB a*b*) distributions of 36 light probe 
% illumination fields used in the experiments by:
% - Loading and converting HDR light probe images to the CIELAB color space,
% - Plotting the chromaticity distribution of each light probe and overlaying the daylight locus,
% - Saving each panel as an individual figure,
% - Assembling all chromaticity plots into a composite grid image.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clear workspace
clearvars; close all;

disp('Generating figure S1...')

load(fullfile('data', 'fig_parameters'))
load(fullfile('data', 'humanANDgt_Pellacini_c'))
load(fullfile('data', 'rankN_lightprobeANDshape'))

%% Lightprobe information
p.lightprobeN = [1 4 5 7 8 9 10 13 16 17 18 20 21 22 23 24 25 28 29 30 31 32 ...
                 33 34 36 37 39 40 41 42 43 45 49 50 54 56];
lightprobe_dir = fullfile('data', 'imgs', 'imgs_lightprobe_hdr');
wp_d65 = [0.9504, 1.0000, 1.0888] * 100;

% a*b* coordinates for daylight locus
ab_daylight = [7.8542, 35.2322; 4.4413, 25.8058; 2.3027, 17.7861; ...
               1.0069, 10.9652; 0.2632, 5.1428; -0.1197, 0.1471; ...
              -0.2674, -4.1634; -0.2624, -7.9201; -0.1588, -11.2005; ...
               0.0068, -14.0852; 0.2103, -16.6379; 0.4356, -18.9097; ...
               0.6718, -20.9422];

%% Plot chromatic distribution for each lightprobe
for lightprobeN = 1:36
    probe_idx = p.lightprobeN(lightprobeN);
    lightprobe = hdrread(fullfile(lightprobe_dir, sprintf('LightProbe%d.hdr', probe_idx)));
    lightprobe(lightprobe < 0) = 0; % set negative values to zero
    
    [h, w, ~] = size(lightprobe);
    
    if h > 512
        scale = 512 / h;
        new_w = round(w * scale);
        lightprobe = imresize(lightprobe, [512, new_w], 'bilinear');
    end

    % Convert RGB to XYZ
    rgb = reshape(lightprobe, [], 3);
    XYZ = SRGBPrimaryToXYZ(rgb')';
    XYZ = XYZ / max(XYZ(:,2)) * 100;

    % Remove dark pixels
    valid = XYZ(:,2) > 0.1;
    rgb = rgb(valid, :);
    XYZ = XYZ(valid, :);
    Lab_lightprobe = XYZToLab(XYZ', wp_d65')';

    %% Plot chromaticity (a*b*) distribution
    fig = figure; hold on
    smpl = 10; % sampling interval
    scatter(Lab_lightprobe(1:smpl*2:end,2), Lab_lightprobe(1:smpl*2:end,3), ...
        2, power(rgb(1:smpl*2:end,:) ./ sum(rgb(1:smpl*2:end,:),2) * 1.5, 1/2.2), ...
        'o', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 0.1);

    plot(ab_daylight(:,1), ab_daylight(:,2), 'k-', 'LineWidth', 0.5)
    scatter(mean(Lab_lightprobe(1:smpl:end,2)), mean(Lab_lightprobe(1:smpl:end,3)), ...
        40, 'kx', 'LineWidth', 1, 'MarkerFaceAlpha', 0.2);

    %% Format figure
    fig.Units = 'centimeters';
    fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    fig.Position = [10, 10, figp.twocolumn / 6, figp.twocolumn / 6];

    ax = gca;
    ax.XLim = [-50, 50];
    ax.YLim = [-50, 50];
    ax.XTick = [-50, 0, 50];
    ax.YTick = [-50, 0, 50];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontName = 'Arial';
    ax.Units = 'centimeters';
    ax.FontSize = figp.fontsize;
    ax.Color = [.97 .97 .97];
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 0.2;
    ax.Position = [0.2, 0.29, 2.7, 2.7];
    ticklengthcm(ax, 0.1)
    grid on; box on

    %% Export figure
    exportgraphics(fig, fullfile('figs', 'temp', ...
        sprintf('figS1c(lightprobe%d_colordistribution.png', lightprobeN)), ...
        'ContentType', 'image', 'Resolution', 300)

    close all
end

%% Combine all chromaticity distribution plots
I = [];
for lightprobeN = 1:36
    img_path = fullfile('figs', 'temp', sprintf('figS1c(lightprobe%d_colordistribution.png', lightprobeN));
    temp = imresize(imread(img_path), [224, 224]);
    I(:,:,:,lightprobeN) = padarray(temp, [20, 80], 255, 'both');
end

I = I(:,:,:,rankN.lightprobe);
s = size(I);
I2 = reshape(permute(reshape(I, [s(1), s(2), s(3), 6, 6]), [1 5 2 4 3]), [s(1)*6, s(2)*6, 3]);

imshow(uint8(I2))
imwrite(uint8(I2), fullfile('figs', 'figS1c(lightprobe_chromaticdistribution).png'))

disp('Done.')
close all