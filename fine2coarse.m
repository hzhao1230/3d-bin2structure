function [ img, Bimg ] = fine2coarse(image, scale)
% This function combines the voxels of a large finescale structure together
% to obtain a smaller coarsescale structure.
% Written by Hongyi Xu
% Inputs:
% image: the fine scale binary image, cube, 1 for filler, 0 for matrix
% scale: the ratio of rescaled size over original size. 320->80, you should
% put 80/320 here.
% Outputs:
% img: the rescaled greyscale image
% Bimg: the rescaled binary image
% Require the function of resize.m

image = int8(image);
L = length(image);
VF = sum( image(:) )/L^3;

LL = round( L * scale );
num = round(1/scale);
img = zeros(LL,LL,LL);
xx = 0;
yy = 0;
zz = 0;

for ii = 1:num:L
%     ii
    xx = xx+1;
    yy = 0;
    for jj = 1:num:L
        
        yy = yy+1;
        zz = 0;
        for kk = 1:num:L
            
            zz = zz+1;
            img(xx,yy,zz) = mean( mean( mean(      image(  ii: min(L,ii+num-1), jj: min(L,jj+num-1), kk:min(L,kk+num-1)  )       ) ) );
            
        end
        
    end
end

VFlist = img(:);
lower = min(VFlist);
upper = max(VFlist);
NN = length(VFlist);

for ii = lower:1:upper-1
    
    threshold1 = sum( VFlist > ii );
    threshold2 = sum( VFlist > ii + 1 );
    if threshold1 >= VF*NN && threshold2 <=VF*NN
        break;
    end
end

Bimg = img > ii;
Bimg = int8(Bimg);

img = img/(1/scale)^3;


temp_img = Bimg;
L = length(Bimg);

AN = VF*L^3 - sum(sum(sum(temp_img)));  % AN: add number

if AN > 0  % Need add "1" voxels
    
    cnt = 0;
    flg = 0;
    while 1
        
        C = find_edge_0_3D(temp_img);
        no_list = randperm( size(C,1) );
        
        for ii = 1:1:size(C,1)
            
            no = no_list(ii);
            x = C(no,1);
            y = C(no,2);
            z = C(no,3);
            temp_img(x,y,z) = 1;
            cnt = cnt + 1;
            
            if cnt >= AN
                flg = 1;
                break;
            end
            
        end
        
        if flg
            break;
        end
        
    end
    
else  % Need add "0" voxels
    
    AN = abs(AN);   
    cnt = 0;
    flg = 0;
    while 1
        
        C = find_edge_1_3D(temp_img);
        no_list = randperm( size(C,1) );
        
        for ii = 1:1:size(C,1)
            
            no = no_list(ii);
            x = C(no,1);
            y = C(no,2);
            z = C(no,3);
            temp_img(x,y,z) = 0;
            cnt = cnt + 1;
            
            if cnt >= AN
                flg = 1;
                break;
            end
        end
        
        if flg
            break;
        end

    end
    
end

Bimg = temp_img;  % img is the final product

