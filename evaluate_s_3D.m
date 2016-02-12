function corrf = evaluate_s_3D(img)

image = img;

pixel = length(image);
% find the locations of edges
image_ext = zeros([pixel+2 pixel+2 pixel+2]); % expanded image
image_ext(2:pixel+1,2:pixel+1,2:pixel+1) = image;
image_tmp = 6*image - image_ext(1:pixel,2:pixel+1,2:pixel+1) - image_ext(3:pixel+2,2:pixel+1,2:pixel+1)...
                    - image_ext(2:pixel+1,1:pixel,2:pixel+1) - image_ext(2:pixel+1,3:pixel+2,2:pixel+1)...
                    - image_ext(2:pixel+1,2:pixel+1,1:pixel) - image_ext(2:pixel+1,2:pixel+1,3:pixel+2);

% in image_tmp
% +: "1" edge
% -: "0" edge
% 0: white/black interior
img_1_edge = (image_tmp > 0);
img = double(img_1_edge);
image = img;
VF = sum( image(:) )/L^3;
clear image

L = length(img);

R = floor( L/2 );
count = zeros(R+1,1);
Bn = zeros(R+1,1);

img = img - mean( mean( mean(img) ) );
F = fftn(img);
c = fftshift(ifftn(F.*conj(F)));

maxV = max( max( max(c) ) );
[ic, jc, kc] = ind2sub( size(img),find(c == maxV) );


for i = 1:1:L
    for j = 1:1:L
        for k = 1:1:L
            r = round( sqrt((i-ic)^2 + (j-jc)^2 + (k-kc)^2 ) ); % r is the distances between a given point and the max value point(peak).
            if r<=R   % this constraints is to confine the counting range into a circle which has the radius of R, and centered at the peak
                Bn(r+1) = Bn(r+1) + c(i,j,k);  % add all the value of the pixels together, which are located on the same circle centered at the peak
                count(r+1) = count(r+1) +1;  % also count the number of the pixels located on the same circle. It should be monotone increasing.
            end
        end
    end
end

Bn = Bn./count; % calculate the average value of each Bn(average value of the pixels located on the same circle)
corrf = Bn./maxV; 

LL = length(corrf);
realcorr = VF^2 + (VF - VF^2) ./ ( corrf(1) - corrf(LL) ) .* ( corrf - corrf(LL) );
corrf = realcorr;
