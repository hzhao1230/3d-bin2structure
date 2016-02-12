function img = invcolor(image)
% Inverse the color of image
image = double(image);
img = 255 - image;

% % Example:
% image = imread('GP_testImage20140122.tif');  % read in a image
% img = invcolor(image);                       % inverse color
% imwrite(img, 'GP_testImage_inverse.tif');    % save to folder
% imshow( uint8(img) );                        % view in MATLAB