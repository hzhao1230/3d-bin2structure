function corrf = evaluate_3D(img)
% This function use the FFT approach in two-point correlation function
% evaluation. It is enlighten by the evaluate.m from Dr. Dikin.
% Written by Hongyi Xu, Northwestern U, 4/23/2013
% 
% Input:
% img: an 3D microstructure image. It should be equilateral.
%
% Output:
% corrf: correlation function curve

% % ----------------------------- For debug ---------------------------------
% load Recon_300X300X300_April17_layers
% img = img(1:50,1:50,1:50);
% img = double(img);
% % -------------------------------------------------------------------------
image = img;
L = length(img);

R = floor( L/2 );
count = zeros(R+1,1);
Bn = zeros(R+1,1);

img = img - mean( mean( mean(img) ) );
F = fftn(double(img));
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

VF = sum( image(:) )/L^3;
LL = length(corrf);
realcorr = VF^2 + (VF - VF^2) ./ ( corrf(1) - corrf(LL) ) .* ( corrf - corrf(LL) );
corrf = realcorr;


