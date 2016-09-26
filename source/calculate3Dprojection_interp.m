%%  calculate3Dprojection_interp %%

%%calculates 2D projection from 3D object using linear interpolation of
%%central slice in Fourier space
%%inputs:
%%  modelK - 3D Fourier space of object to compute projection of
%%  phi,theta,psi - Euler angles of desired projection

%%outputs:
%%  projection - result

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.



function projection = calculate3Dprojection_interp(modelK,phi,theta,psi)
%%only works for even images in current form%%
%%%%

%get dimensions and centers
[dimx, dimy, dimz] = size(modelK);
if mod(dimx,2)==1 || mod(dimz,2)==1 || mod(dimy,2)==1
   error('GENFIRE: ERROR! All array dimensions must be even!\n\n') 
end

ncy = round((dimy+1)/2); ny2 = ncy-1;
ncx = round((dimx+1)/2); nx2 = ncx-1;
ncz = round((dimz+1)/2); nz2 = ncz-1;

[Y, X, Z] = meshgrid(-ny2:ny2-1,-nx2:nx2-1,0);

%calculate rotation matrix
R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
      -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
      sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];

  
[ky kx kz ] = meshgrid(-ny2:ny2-1,-nz2:nz2-1,-nz2:nz2-1);

%rotate coordinates
rotkCoords = R'*[X(:)';Y(:)';Z(:)'];
rotKx = rotkCoords(1,:);
rotKy = rotkCoords(2,:);
rotKz = rotkCoords(3,:);

%reshape for interpolation
rotKx = reshape(rotKx,size(X));
rotKy = reshape(rotKy,size(Y));
rotKz = reshape(rotKz,size(Z));

%calculate points on central slice
% pjK = interp3(ky,kx,kz,modelK,rotKy,rotKx,rotKz,'linear');
pjK = interp3(ky,kx,kz,modelK,rotKy,rotKx,rotKz,'cubic');

%remove any nan from interpolation
pjK(isnan(pjK))=0;

%take IFFT to obtain projection
projection = real(my_ifft(pjK(:,:)));
end