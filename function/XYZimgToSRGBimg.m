function sRGB = Gloss_XYZimgToSRGBimg(xyz)
    s = size(xyz);
    srgb = XYZToSRGBPrimary(reshape(xyz,s(1)*s(2),s(3))');
    srgb_out = SRGBGammaCorrect(srgb);
    sRGB = reshape(double(srgb_out')/255,s);
    %imshow(sRGB)
end