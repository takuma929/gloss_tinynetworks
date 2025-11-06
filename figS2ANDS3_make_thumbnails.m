%--------------------------------------------------------------------------
% This script generates Supplementary Figures S1a and S1b for the manuscript.
% It creates and saves thumbnail grids of all light probes and object shapes
% used in the experiments by:
% - Loading and resizing each light probe and shape image,
% - Padding and assembling images into 6×6 grids (36 items each),
% - Flipping specific object images for consistent presentation,
% - Saving the resulting composite images in 'figs' directory.
%
% Author: TM, 2025
%--------------------------------------------------------------------------

%% Clear workspace
clearvars; close all;

load(fullfile('data', 'rankN_lightprobeANDshape'), 'rankN')

%% --- Create thumbnails for lightprobes ---
lightprobe_indices = [1 4 5 7 8 9 10 13 16 17 18 20 21 22 23 24 25 28 ...
                      29 30 31 32 33 34 36 37 39 40 41 42 43 45 49 50 54 56];

lightprobe_dir = fullfile('data', 'imgs', 'imgs_lightprobe');
cnt = 0;

for lightprobeN = rankN.lightprobe
    cnt = cnt + 1;
    probe_idx = lightprobe_indices(lightprobeN);
    img_path = fullfile(lightprobe_dir, sprintf('LightProbe%d.png', probe_idx));
    img = imresize(imread(img_path), [128 256]);
    I.lightprobe(:,:,:,cnt) = padarray(img, [10, 10], 255, 'both');
end

%% Assemble and save lightprobe thumbnails
w = 6; h = 6; % grid dimensions (6x6 = 36)
s = size(I.lightprobe);

out.lightprobe = reshape(permute(reshape(I.lightprobe, [s(1), s(2), s(3), w, h]), [1 5 2 4 3]), [s(1)*h, s(2)*w, s(3)]);
imwrite(out.lightprobe, fullfile('figs', 'figS1b(thumbnail_lighting).png'));
imshow(out.lightprobe);pause(2);
close all

%% --- Create thumbnails for shapes ---
shape_dir = fullfile('data', 'imgs', 'imgs_obj_mesh');
cnt = 0;

for shapeN = rankN.shape
    cnt = cnt + 1;
    img_path = fullfile(shape_dir, sprintf('shape%d.png', shapeN));
    img = imread(img_path);
    
    % Flip shape 21 and 25 horizontally
    if ismember(shapeN, [21, 25])
        img = flip(img, 2);
    end
    I.shape_obj(:,:,:,cnt) = img;
end

%% Assemble and save shape thumbnails
w = 6; h = 6;
s = size(I.shape_obj);

out.shape_obj = reshape(permute(reshape(I.shape_obj, [s(1), s(2), s(3), w, h]), [1 5 2 4 3]), [s(1)*h, s(2)*w, s(3)]);
imwrite(out.shape_obj, fullfile('figs', 'figS1a(thumbnail_obj).png'));
imshow(imresize(out.shape_obj, 0.5));pause(2);
close all
