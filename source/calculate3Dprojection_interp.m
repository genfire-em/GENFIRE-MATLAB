%%  calculate3Dprojection_interp %%

%%calculates 2D projection from 3D object using linear interpolation of
%%central slice in Fourier space
%%inputs:
%%  modelK - 3D Fourier space of object to compute projection of
%%  phi,theta,psi - Euler angles of desired projection

%%outputs:
%%  projection - result

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.



function projection = calculate3Dprojection_interp(modelK,phi,theta,psi)

%get dimensions and centers
[dimx, dimy, dimz] = size(modelK);

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
ncz = round((dimz+1)/2);

[Y, X, Z] = meshgrid((1:dimy) - ncy, (1:dimx) - ncx, 0);

%calculate rotation matrix
R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
      -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
      sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];

  
[ky kx kz ] = meshgrid((1:dimy) - ncy, (1:dimx) - ncx, (1:dimz) - ncz);

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