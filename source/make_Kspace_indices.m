%%  make_Kspace_indices %%

%% compute coordinate indices for input matrix

%%inputs:
%%  obj - input matrix of the desired size to map

%%outputs:
%%  Q - X, Y, and Z coordinates of voxels corresponding to a matrix of size(obj)

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function Q = make_Kspace_indices(obj)

%construct K-space indices
if mod(size(obj,1),2)==0
ncK1 = size(obj,1)/2+1;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1-1)./n2K1;
elseif size(obj,1)==1
vec1 = 0;
ncK1 = 1;
n2K1 = 0;
else
ncK1 = (size(obj,1)+1)/2;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1)./n2K1; 
end

if  mod(size(obj,2),2)==0
ncK2 = size(obj,2)/2+1;%central pixel
n2K2 = ncK2-1;%max radius
vec2 = (-n2K2:n2K2-1)./n2K2;
elseif size(obj,2)==1
vec2 = 0;
ncK2 = 1;
n2K2 = 0;
else
ncK2 = (size(obj,2)+1)/2;%central pixel
n2K2 = ncK2-1;%max radius
vec2 = (-n2K2:n2K2)./n2K2; 
end

if  mod(size(obj,3),2)==0
ncK3 = size(obj,3)/2+1;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3-1)./n2K3;
elseif size(obj,3)==1
vec3 = 0;
ncK3 = 1;
n2K3 = 0;
else
ncK3 = (size(obj,3)+1)/2;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3)./n2K3; 
end
[Kx Ky Kz] = meshgrid(vec2,vec1,vec3);%grid of Fourier indices
Q = sqrt(Kx.^2+Ky.^2+Kz.^2);%take magnitude
end