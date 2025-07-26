function Out = gloss_spectraTotristimulusvals_400to720(HSI)

% Size of Input
imagesize = size(HSI);
% Vectorization of HyperSpectralImage
if length(size(HSI)) == 3
    HSI_reshaped = reshape(HSI,imagesize(1)*imagesize(2),imagesize(3));
elseif length(size(HSI)) < 3
    HSI_reshaped = HSI;
end

% Load Color Matching Function
load('lms_400to720.mat');
load('T_xyz1931.mat');

% White Point
%wp  = cie_whitepoint('D65');
wp = [0.9504 1.0000 1.0888];

% Parameters for CIECAMO2
Yb = 20;    % NumericScalar, relative luminance of reference white in the adapting field.  - Average across 3 or 4 angles / Average across whole scenes  
LA = 64/pi/5;   % NumericScalar, adapting field luminance (cd/m^2).   - Average across whole scenes 

% Spline xyz to make it 400nm to 720 nm by 10 nm step
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,[400,10,33])';

% Weighting for LMS cone signals
Lw = 0.689903;
Mw = 0.348322;
Sw = 0.0371597/0.0192; 

XYZ_reshaped = HSI_reshaped*T_xyz*683*10; % Spctrum to XYZ

sRGB_reshaped = XYZToSRGBPrimary(XYZ_reshaped')'; % XYZ to sRGB
sRGB_reshaped = max((sRGB_reshaped/max(max(sRGB_reshaped))),0);
sRGB_reshaped = sRGB_reshaped.^(1/2.4);

Lab_reshaped = XYZToLab(XYZ_reshaped',wp')'; % XYZ to Lab

% Reshape All Chromatic Coordinates
if length(size(HSI)) == 3
    Out.XYZ = reshape(XYZ_reshaped,imagesize(1),imagesize(2),3);
    Out.Lab = reshape(Lab_reshaped,imagesize(1),imagesize(2),3);
    Out.sRGB = reshape(sRGB_reshaped,imagesize(1),imagesize(2),3);

elseif length(size(HSI)) < 3
    Out.XYZ = XYZ_reshaped;
    Out.Lab = Lab_reshaped;
    Out.sRGB = sRGB_reshaped;
end

Out.XYZ(isnan(Out.XYZ))=0;
Out.Lab(isnan(Out.Lab))=0;
Out.sRGB(isnan(Out.sRGB))=0;
end
