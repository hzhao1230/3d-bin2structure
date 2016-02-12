function MAIN(img_name, type, color, sphere, rb, rl, ru, cutL, VF, recon_length, scale, visualize_fine_recon, visualize_coarse_recon, workingdir)
%% Introduction
%
% It is the main function for 3D descriptor-based microstructure
% characterization and reconstruction.
% To use the 3D descriptor-based C&R software package, the user only need 
% to open this file and click on "Run".
%
% There are two parts that need user's operation: 
%   I. Define inputs (SEM image name, VF, size of reconstruction, size of 
%      small coarsened reconstruction, 3D plots yes/no).
%   II. Operation IMAGEJ.
%
% This package is developed by Hongyi Xu, Northwestern University
% PI: Prof. Wei Chen, Northwestern University
% Contact: hongyixu2014@u.northwestern.edu


%% I: User-defined info. User-operation required (1 out of 2)
% ------------------------ User-Define Region -----------------------------
% Three basic inputs: name of SEM image (.tif), volume fraction, reconstruction size
% % Example:
% img_name = 'GP_testImage20140122';
% type = 1;
% color = 0;
% sphere = 0; 
% cutL = 500;
% VF = 0.31; % Volume fraction
% recon_length = 300; % Side length of reconstruction (cube)
% scale = 60/300;  % Scaling factor for generate an coarsened reconstruction: from 300 to 60
% workingdir = 'C:\Users\Hongyi Xu\Working\MATLAB_working\Characterize_Morozov\Deliverable_3D_Desciptor_v3.3';
% % Whether plot the 3D reconstruction? 0-no, 1-yes
% visualize_fine_recon = 0;     % For the 300 size fine reconstruction. It may take a long time to plot the large structure
% visualize_coarse_recon = 1;   % For the 60 size coarsened reconstruction
% % -------------------------------------------------------------------------
cd(workingdir)  % go to the working directory

img_name_inv = img_name;
if color == 1  % if background color is brighter than filler (white matrix)
    image = imread([img_name,'.tif']);
    image = invcolor(image);
    imwrite(uint8(image), [img_name,'_inv.tif']);
    img_name_inv = [img_name, '_inv'];
    clear image
end

%% II: Characterization. User-operation required (2 out of 2)
if type == 0 % if it is a binary image
    
    Descriptor_C2_Binary(img_name_inv, VF, recon_length, cutL);
    
else  % If it is a greyscale image

    % Step 1: image preprocessing
    Descriptor_C1(img_name_inv, VF, cutL);

    % Step 2: 2D characterization and 3D prediction
    Descriptor_C2(img_name_inv, VF, recon_length);
    
end

%% III: Reconstruction
cd(workingdir)  % go to the working directory

load( [ img_name_inv, '_GB_double_filter_3D_results.mat' ] ); % Load characterization results, which include a variable "name" (= img_name)
clearvars -except  name  recon_length  VF  num_3D  ND3D  Predict_3D_As_mean  Predict_3D_As_var scale  visualize_fine_recon  visualize_coarse_recon  type  sphere  rb  rl  ru% only keep the useful variables

% Start microstrutcure reconstruction
name_3D = [name, '_3D_recon'];
Descriptor_Recon_3D(name_3D, recon_length, VF, num_3D, ND3D, Predict_3D_As_mean, Predict_3D_As_var, sphere, rb, rl, ru);

clearvars -except  name  name_3D  scale  visualize_fine_recon  visualize_coarse_recon  type  img_name  sphere  rb  rl  ru
% clc

disp('The reconstruction is finished! Now working on file moving, corrlation function calculation and plotting ...')

% %% IV: Rescale the large reconstruction into a small one
% % Here the 300x300x300 is rescaled into a smaller 60x60x60 structure (for FEA)
% load( [name_3D, '.mat'] )
% % clear -except  name img  scale  visualize_fine_recon  visualize_coarse_recon
% 
% [ img_coarse, Bimg_coarse ] = fine2coarse(img, scale);
% save( [ name, '_coarsen_recon.mat'], 'img_coarse', 'Bimg_coarse' )


%% V: Create a folder (the same name as the image), and move results into it
mkdir([ name, '_results' ]);
movefile( [name, '_GB_double_filter_2D_results.mat'], [ name, '_results' ])
movefile( [name, '_GB_double_filter_3D_results.mat'], [ name, '_results' ])
movefile( [name, '_3D_recon.mat'], [ name, '_results' ])
% movefile( [ name, '_coarsen_recon.mat'], [ name, '_results' ])
if type == 1
    movefile( [ name, '_GB_double_filter.tif'], [ name, '_results' ])
end
movefile( [ name_3D, '_center_list.mat'], [ name, '_results' ])
movefile( [ name_3D, '_3D_orinatation.mat'], [ name, '_results' ])
movefile( [ name_3D, '_3D_geometry.mat'], [ name, '_results' ])

%% VI: Point correlation & Surface correlation
% disp('Calculating the point and surface correlation function ...')
% corrf = evaluate_3D(img);
% corrs = evaluate_s_3D(img);
% figure('color',[1,1,1]);
% hold on
% plot( 0:1:length(corrf)-1, corrf )
% plot( 0:1:length(corrf)-1, corrs, 'r' )
% xlabel('Distance (Voxel)')
% ylabel('Correlation Function')
% box on
% legend('2-point correlation function', 'Surface correlation function')
% save( [ name, '_corr_func.mat'], 'corrf', 'corrs' )
% movefile( [ name, '_corr_func.mat'] , [ name, '_results' ])

% %% VII: Visualization
% if visualize_coarse_recon == 1
%     cord = ThreeD2Coordinate(Bimg_coarse);
%     figure('color',[1,1,1]);
%     voxel_image( cord, 1, 'g', 0.5, 'b' ); 
%     axis equal
%     box on
%     L = length(Bimg_coarse);
%     axis([1 L 1 L 1 L])
%     view([1,0.5,0.5])
% end
% 
% if visualize_fine_recon == 1
%     cord = ThreeD2Coordinate(img);
%     figure('color',[1,1,1]);
%     voxel_image( cord, 1, 'g', 0.5, 'b' ); 
%     axis equal
%     box on
%     L = length(img);
%     axis([1 L 1 L 1 L])
%     view([1,0.5,0.5])
% end
%  
% clc
% disp('The 3D descriptor-based C&R is completed!')
% 
