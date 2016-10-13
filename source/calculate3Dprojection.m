%%  calculate3Dprojection %%

%%calculates 2D projection from 3D object accurately, but slowly, using DFT
%%method
%%inputs:
%%  model - 3D object to compute projection of
%%  phi,theta,psi - Euler angles of desired projection

%%outputs:
%%  projection - result

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.


function projection = calculate3Dprojection(model,phi,theta,psi)


%get dimensions and centers
[dimx dimy dimz] = size(model);

ncy = round((dimy+1)/2);
ncx = round((dimx+1)/2); 
ncz = round((dimz+1)/2); 

%calculate rotation matrix
R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
      -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
      sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];

%initialize coordinates
[ky kx] = meshgrid((1:dimy) - ncy, (1:dimz)-ncz);
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