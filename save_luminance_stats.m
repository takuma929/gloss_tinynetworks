%--------------------------------------------------------------------------
% This optional script processes rendered images and computes summary luminance statistics.
% This allows the reproduction of ImageStats.mat file provided in the 'data' directory. 
% For each of 3888 images, it:
%   - Loads the RGB image and corresponding object mask,
%   - Applies the mask to extract object pixels,
%   - Converts sRGB to linear RGB and then to XYZ color space,
%   - Computes luminance (Y) for object pixels,
%   - Calculates summary statistics (mean, min, max, std, median, quartiles, skewness, kurtosis, contrast)
%     for each image's luminance distribution,
%   - Stores results in the ImageStats structure for further analysis or modeling.
%--------------------------------------------------------------------------

clearvars;close all;clc

img_dir = fullfile('data','imgs','imgs3888_bg_512_png');
mask_dir = fullfile('data','imgs','obj_mask');

for imgN = 1:3888
    imgN
    I = double(imread([img_dir,'/img',num2str(imgN),'.png']))/255;
    I_mask = double(imread(fullfile(mask_dir,['img',num2str(imgN),'_mask.png'])))/255;

    obj = I.*repmat(I_mask,1,1,3);
    srgb_linear = SRGBGammaUncorrect(uint8(obj*255));
    s = size(srgb_linear);
    XYZ = SRGBPrimaryToXYZ(reshape(srgb_linear,[s(1)*s(2),s(3)])');
    
    Y = nonzeros(XYZ(2,:));
   
    % glossiness predictor
    ImageStats.mean(imgN) = mean(Y);
    ImageStats.min(imgN) = min(Y);
    ImageStats.max(imgN) = max(Y);    
    ImageStats.std(imgN) = std(Y);
    ImageStats.median(imgN) = median(Y);
    ImageStats.Q1(imgN) = prctile(Y,25);
    ImageStats.Q3(imgN) = prctile(Y,75);
    ImageStats.kurtosis(imgN) = kurtosis(Y);
    ImageStats.skewness(imgN) = skewness(Y);
    ImageStats.contrast(imgN) = (max(Y)-min(Y))/(max(Y)+min(Y));
end

%save(fullfile('data','ImageStats'),'ImageStats')
