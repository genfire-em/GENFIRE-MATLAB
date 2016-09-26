%%  calculate3Dprojection %%

%%calculates 2D projection from 3D object accurately, but slowly, using DFT
%%method
%%inputs:
%%  model - 3D object to compute projection of
%%  phi,theta,psi - Euler angles of desired projection

%%outputs:
%%  projection - result

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.


function projection = calculate3Dprojection(model,phi,theta,psi)

%%only works for even images in current form%%

%get dimensions and centers
[dimx dimy dimz] = size(model);
if mod(dimx,2)==1 || mod(dimz,2)==1 || mod(dimy,2)==1
   error('GENFIRE: ERROR! All array dimensions must be even!\n\n') 
end

ncy = round((dimy+1)/2); ny2 = ncy-1;
ncx = round((dimx+1)/2); nx2 = ncx-1;
ncz = round((dimz+1)/2); nz2 = ncz-1;

%calculate rotation matrix
R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
      -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
      sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];

%initialize coordinates
[ky kx] = meshgrid(-ny2:ny2-1,-nz2:nz2-1);
kx = single(kx(:))'; ky = single(ky(:))'; %initialize coordinates of unrotate projection slice
kz = zeros(1,dimx*dimy,'single'); %0 degree rotation is a projection onto the X-Y plane, so all points have kz=0;

rotkCoords = R'*[kx;ky;kz];%rotate coordinates
rotkx = rotkCoords(1,:);
rotky = rotkCoords(2,:);
rotkz = rotkCoords(3,:);
kOut = zeros(size(rotky),'single');

%only consider nonzero indices for efficiency
ind = find(model~=0);
[x y z] = ind2sub(size(model),ind);x = x-ncx; y = y-ncy; z = z-ncz;

for jj=1:length(kOut) %compute DFT
    kOut(jj) = sum(model(ind).*exp(-2*pi*1i*(rotkx(jj)*x/dimx + rotky(jj)*y/dimy + rotkz(jj)*z/dimz)));
end

%reshape to correct dimensions
kOut = reshape(kOut,dimx,dimy);

%take IFFT to obtain projection
projection = (real(my_ifft(kOut)));
end