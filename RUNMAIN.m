%% Create 3D recon parameters in a separte folder using Hongyi's code. 

img_name = 'crop_terthiophenePGMA_2wt%'; % black matrix, white NP
type = 0;
color = 0;
sphere = 0;
cutL = 1036; 
rb = 0;
rl = 0; ru = 0;
VF = 0.01; 
recon_length = 300; 
scale = 1/4;
plot_3D = 0;
plot_coarse = 0;
workingdir = './';

main(img_name, type, color, sphere, rb,rl,ru, cutL,VF, recon_length, scale, plot_3D, plot_coarse, workingdir)

%% Combine x, y, z positions, long, short axes, y, z rotation angles to a
% single matrix. 

wdir = ['./',img_name,'_results'];
cd(wdir)

% Get center list
load([img_name,'_3D_recon_center_list']);
x = cl(:,1); y = cl(:,2); z=cl(:,3);
% Get long and short axes
load([img_name,'_3D_recon_3D_geometry']);
la = geo_mat(:,2); sa = geo_mat(:,3);
% Get orientation angle
load([img_name,'_3D_recon_3D_orinatation']);
oy = orint_ang(:,1); oz = orint_ang(:,2); 

img_para = [x,y,z,la,sa,oy,oz];

save([img_name,'_3D_structure_output'],'img_para');