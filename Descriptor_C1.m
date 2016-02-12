function Descriptor_C1(name, VF, cutL)
% ------------------------------------------------------------------------- 
% Function used:
% overlap_particle_calc.m
% area_stat.m
% faster_elongation.m
% Transform.m
% -------------------------------------------------------------------------
disp('Preprocessing the image. It may take several minutes.')
image = imread([name, '.tif']);  % The image used in previous image analysis, 330-60-100-0
image = double(image);

if cutL > length(image);
    cutL = length(image);
end

image = image(1:cutL,1:cutL);

img = medfilt2(image, [5 5]);

% h=fspecial('average', [ceil(cutL/10) ceil(cutL/10)]);
% h=fspecial('average', [50 50]);
% image_f = imfilter(image, h);

image_f = medfilt2(image, [50 50]);

L = length(image);

img = img - image_f;
[~, Bimg] = Transform(img, 1-VF);    
Bimg = abs( Bimg/255 - 1 );
GBimg = Bimg .* double( image );  % overlay gray scale and binary image

% -------------- filter again to guarantee smoothness ---------------------
GBimg = medfilt2(GBimg, [5 5]);
[~, Bimg] = Transform(GBimg, 1-VF);
Bimg = abs( Bimg/255 - 1 );
GBimg = Bimg .* double( GBimg ); 
% -------------------------------------------------------------------------

GBimg = uint8(GBimg);

imwrite(GBimg,[ name, '_GB_double_filter.tif'] );

