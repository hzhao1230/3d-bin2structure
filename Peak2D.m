function [ ploc, cplot ] = Peak2D(image, plotctrl)
% -------------------------------------------------------------------------
% Find local maximas in 2D image
%
% Input:
% An image, or a matrix. 
% It is suggested to use square images.
% plotctrl: "1" gives a plot; "0" gives no plot
%
% Output:
% ploc: list of peak locations (xy coordinates)
% cplot: binary image of peaks
%
% Comments:
% Before using this function, the image should be pre-filtered (use a
% little longer filtering length)
%
% By Hongyi Xu, Northwestern U, Feb 10, 2014
% -------------------------------------------------------------------------
image = double(image);
image2 = image;

% image = medfilt2(image, [5 5]);

% % Filtering (fast average filter), 10 pixel
% h=fspecial('average', [10 10]);
% image = imfilter(image, h);

[Lx, Ly] = size(image);

mimg = zeros(Lx, Ly);  % Mark image: marking peaks

for ii = 1:1:Lx
    
    data = image(ii,:);   % Get a line
    [~,locs] = findpeaks(data);
    locs = locs(:);
    
    if data(1) > data(2)
        locs = [ locs; 1];
    end
    if data(Ly) > data(Ly-1)
        locs = [ locs; Ly ];
    end
    
    tline = zeros(1, Ly);
    tline(locs) = 1;
    
    mimg(ii,:) = mimg(ii,:) + tline;
    
end

for ii = 1:1:Ly
    
    data = image(:,ii);
    [~,locs] = findpeaks(data);
    locs = locs(:);
    
    if data(1) > data(2)
        locs = [ locs; 1];
    end
    if data(Lx) > data(Lx-1)
        locs = [ locs; Lx ];
    end
    
    tcolm = zeros(Lx, 1);
    tcolm(locs) = 1;
    
    mimg(:,ii) = mimg(:,ii) + tcolm;
    
end

cplot = zeros(Lx+4, Ly+4);
ploc = [];
for ii = 1:1:Lx
    for jj = 1:1:Ly
        
        cx = ii + 2;
        cy = jj + 2;
        if mimg(ii,jj) == 2 ...
                && sum( sum( cplot( cx-2:cx+2, cy-2:cy+2 ) ) ) == 0
%                 && cplot(cx-1, cy-1)~=1 && cplot(cx-1, cy+1)~=1 && cplot(cx+1, cy-1)~=1 && cplot(cx+1, cy+1)~=1 ...
%                 && cplot(cx-1, cy)~=1 && cplot(cx+1, cy)~=1 && cplot(cx, cy-1)~=1 && cplot(cx, cy+1)~=1
            cplot(cx, cy) = 1;
            ploc = [ ploc; ii, jj ];
        end
        
    end
end

cplot = cplot(3:Lx+2, 3:Ly+2);

% Plot part
if plotctrl ~= 0
    imagesc(image2); hold on
    plot(ploc(:,2),ploc(:,1),'k+')
    axis equal
end
