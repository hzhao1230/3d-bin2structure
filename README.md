# 3d-bin2structure
converts 2D binary image (black/0 as matrix) to 3D structure (N by 7 array) for FEA input

# Reference
Hongyi Xu, et al., Descriptor-based methodology for statistical characterization and 3D reconstruction of microstructural materials, Computational Materials Science (2014)
<http://www.sciencedirect.com/science/article/pii/S0927025613008057>

Hongyi Xu <hongyixu2014@u.northwestern.edu>


# Usage
Run RUNMAIN

# Parameters 
img_name: MAT file containing binary structure

type: 0 for binary 

color: 0 for dark background matrix

sphere: 0 for elliptical clusters

cutL: side length of square to be cut from binary image. In practice, use pixel side length of binary image. 

VF: volume fraction (0.01 for 1%)

recon_length: voxel side length of 3D reconstruction cube. In practice, choose a small size (~300-500) to prevent RAM overload. 

scale: coarsening ratio (e.g., 60/300 denotes to coarsen a 300X300X300 cube to side length of 60). In practice, choose 1/4 to 1/5. 

# Output
_output: NX7 array (N as number of clusters), with each column as X, Y, Z, long axis, short axis, Y orientation, Z orientation. 
