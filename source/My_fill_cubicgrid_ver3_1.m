% My_fill_cubicgrid
% Y. Yang, UCLA Physics & Astronomy
% First version date: 2015. 04. 30.
% output parameter: ftArray (nx x ny x ny) Fourier 3D array
% input parameter: nx_ori (original projection 1st dim size)   
%                  ny_ori (original projection 2nd dim size)                  
%                  PROJvol - measured projection, (nx x ny x P) array with P # of projections
%                  angles - 3 Euler angle array, (3 x P) array
%                  K_thresh - threshold for acceptible distance
%                  Typeind - 1 for picking one minimal distance value
%                            2 for weighted averaging all values within threshold
%                   nx (length of oversampled projection 1st dim)
%                    ny (length of oversampled projection 2nd dim) 
%                   ThreshType - 1 for fixed interpolation threshold    
%                                2 for spacial frequency dependent
%                                threshold
%
% Second version date: 2015. 7. 7. (YY)
% Change: 1. Now this code properly process arbitrary-sized input projections
%         and oversampling. nx_ori, ny_ori, nx, ny can be arbitrary positive
%         integer, and it does not matter if they are even or odd.
%         2. This code assumes that the pixels of original projections are
%         "square", i.e. the pixel resolution is the same in both
%         directions. This version of code does not work if pixels are rectangular.
%         3. Implemented spacial frequency dependent interpolation
%         threshold (turn on and off using ThreshType argument)
%
% Thrid version date: 2015. 7. 9. (YY)
% Change: 1. Modified dimensional convention (1st dim: x, 2nd dim: y, 
%                                            3rd dim: z)
%         2. Rotation matrix is now Z(phi)*X(theta)*Z(psi), instead of
%         previous version [Z(phi)*Y(theta)*Z(psi)].
%
%         3. Zero degree projection is in XY plane, not XZ plane as before.

function ftArray = My_fill_cubicgrid_ver3_1(nx_ori, ny_ori, PROJvol, angles,K_thresh,Typeind, nx, ny, ThreshType)

% if distance below minInvThresh, minInvThresh will be used
% this is to prevent division by zero
minInvThresh = 0.001;

% initialize normal vectors and rotation matrices array
normVECs = zeros(size(PROJvol,3),3);
rotMATs = zeros(3,3,size(PROJvol,3));

phis = angles(1,:);
thetas = angles(2,:);
psis = angles(3,:);

% calculate rotation matrices and normal vectors from the rotation matrices
for i=1:size(PROJvol,3)
    vector1 = [0 0 1];
    rotmat1 = MatrixQuaternionRot(vector1,phis(i));
    
    vector2 = [0 1 0];
    rotmat2 = MatrixQuaternionRot(vector2,thetas(i));

    vector3 = [0 0 1];
    rotmat3 = MatrixQuaternionRot(vector3,psis(i));

    rotMATs(:,:,i) =  (rotmat1*rotmat2*rotmat3);

    init_normvec = [0 0 1];
    normVECs(i,:) = squeeze(rotMATs(:,:,i))*init_normvec';

 R = [-sind(psis(i))*cosd(thetas(i))*sind(phis(i))+cosd(psis(i))*cosd(phis(i)), -sind(psis(i))*cosd(thetas(i))*cosd(phis(i))-cosd(psis(i))*sind(phis(i)), sind(psis(i))*sind(thetas(i));
   cosd(psis(i))*cosd(thetas(i))*sind(phis(i))+sind(psis(i))*cosd(phis(i)) , cosd(psis(i))*cosd(thetas(i))*cosd(phis(i))-sind(psis(i))*sind(phis(i)) ,-cosd(psis(i))*sind(thetas(i));
 sind(thetas(i))*sind(phis(i)), sind(thetas(i))*cosd(phis(i)),cosd(thetas(i))];
end


% initiate Fourier space indices
if mod(nx,2) == 1
    kx_half = -1*(nx-1)/2:1:0;
    nx_half = (nx-1)/2;
else
    kx_half = -1*nx/2:1:0;
    nx_half = nx/2;
end
   
% initiate Fourier space indices
if mod(ny,2) == 1
    ky = -1*(ny-1)/2:1:(ny-1)/2;
    ky_ext = ky;
else
    ky_ext = -1*ny/2:1:ny/2;
end

kz_ext = ky_ext;
kx_half_1_ori = kx_half(1);
kx_half = kx_half * ky_ext(1) / kx_half_1_ori;

% Fourier grid
% To save time, only half of kz will be interpolated
% and centrosymmetricity will be enforced later
[KY, KX, KZ] = meshgrid(ky_ext,kx_half,kz_ext);

% initialize variables
FS = single(zeros(size(KX))*(-1)); % Fourier points
Numpt = single(zeros(size(KX))); % array to store how many points found per Fourier point
Mindist = single(ones(size(KX))*(10000)); % array to store current minimum distance
invSumdist = single(zeros(size(KX))); % array to store sum of inverse distance of Fourier point

% initiate Fourier space indices
if mod(nx_ori,2) == 1
    kx_ori = -1*(nx_ori-1)/2:1:(nx_ori-1)/2;    
else
    kx_ori = -1*nx_ori/2:1:nx_ori/2-1;
end

% initiate Fourier space indices
if mod(ny_ori,2) == 1
    ky_ori = -1*(ny_ori-1)/2:1:(ny_ori-1)/2;
else
    ky_ori = -1*ny_ori/2:1:ny_ori/2-1;
end


[KY0, KX0] = meshgrid(ky_ori,kx_ori);
    
for p=1:size(PROJvol,3)
    %p
    % current projection
    curr_proj = squeeze(PROJvol(:,:,p));

    [KY, KX, KZ] = meshgrid(ky_ext,kx_half,kz_ext);
    % obtain points-to-plane distance
    D = distancePointsPlane_YY([KX(:) KY(:) KZ(:)]', normVECs(p,:));
    
    % find Fourier points within the threshold
    
    if ThreshType == 1
        % fixed K_thresh
        Dind = find(D < K_thresh);
    elseif ThreshType == 2        
        % K_thresh as a function of spacial frequency
        KX_reg = KX/ky_ext(1)*kx_half_1_ori;
        K_mag = sqrt(KX_reg.^2+KY.^2+KZ.^2);
        K_weighted_thresh = K_thresh / n * K_mag + 0.5;    
        
        Dind = find(D < K_weighted_thresh(:)');
    end
    
    % rotate the plane to zero degree
    CP = closestpoint(normVECs(p,:)',0,[KX(Dind); KY(Dind); KZ(Dind)]);
    CP_plane = (squeeze(rotMATs(:,:,p)))\CP;
    
    clear KX KY KZ
    % picked closest point which is within the projection plain, x coordinate must be zero after rotation
    if sum(abs(CP_plane(3,:)) > 0.0001) > 0
        fprintf(1,'something wrong!\n');                
    end
    
    % consider Fourier points only within the resolution circle
    Gind = Dind(abs(CP_plane(1,:)) <= nx/2*ky_ext(1)/ kx_half_1_ori & abs(CP_plane(2,:)) <= ny/2);  % good indices  
    G_CP_plane = CP_plane(:,abs(CP_plane(1,:)) <= nx/2*ky_ext(1)/kx_half_1_ori & abs(CP_plane(2,:)) <= ny/2 );  % good in-plane coordinates
    
    % determine the available memory in MB
    if ispc % in case of Windows machine
        [~, SYSTEMVIEW] = memory;
        av_memory_size = SYSTEMVIEW.PhysicalMemory.Available / 1000000;
    elseif isunix && ~ismac  % in case of linux (or unix)
        [~,out]=system('cat /proc/meminfo | grep MemFree');
        av_memory_size=sscanf(out,'MemFree:          %d kB');
        av_memory_size = av_memory_size / 1000;
    else % in case of mac (I don't have mac now, to be implemented later)
        av_memory_size = 1000;
    end
    
    memory_per_index = 40*length(curr_proj(:))/1000000;        
    
    % determine block size for vectorized calculation
    block_size = floor(av_memory_size/memory_per_index);
    if block_size < 1
        block_size = 1;
    end
    %block_size = 500;
    cutnum = floor(length(Gind)/block_size);
    cutrem = mod(length(Gind),block_size);
    
    if cutrem~=0
        cutloopnum = cutnum + 1;
    else
        cutloopnum = cutnum;
    end
    
    % loop over Fourier points within the threshold
    for i=1:cutloopnum        
        curr_indices = ((i-1)*block_size+1):(i*block_size);
        
        if i > cutnum
            curr_indices = (cutnum*block_size+1):length(Gind);
        end
        
        % DFT calculation
        Fpoints = sum(bsxfun(@times, curr_proj(:), exp(-1*1i*2*pi*(KX0(:)*G_CP_plane(1,curr_indices)/nx/ky_ext(1)*kx_half_1_ori+KY0(:)*G_CP_plane(2,curr_indices)/ny))),1);
        
        if Typeind == 1 % minimum distance point
            Minind = find(D(Gind(curr_indices)) < Mindist(Gind(curr_indices)));
            FS(Gind(curr_indices(Minind))) = Fpoints(Minind);
            Mindist(Gind(curr_indices(Minind))) = D(Gind(curr_indices(Minind)));
        elseif Typeind == 2  %weighted avearaging
            currDist = D(Gind(curr_indices));
            currDist(currDist < minInvThresh) = minInvThresh; % if distance smaller than minInvThresh, put minInvThresh (to prevent divison by zero)
            currInvDist = 1./ currDist;           % inverse distance
            
            % re-average inverse distance
            FS(Gind(curr_indices)) = FS(Gind(curr_indices)).* invSumdist(Gind(curr_indices)) + currInvDist.*Fpoints;
            invSumdist(Gind(curr_indices)) = invSumdist(Gind(curr_indices)) + currInvDist;
            FS(Gind(curr_indices)) = FS(Gind(curr_indices)) ./ invSumdist(Gind(curr_indices));
            Numpt(Gind(curr_indices)) = Numpt(Gind(curr_indices)) + 1;
        end
            
        clear Fpoints
    end
end

clear Mindist invSumdist Numpt

% enforce centrosymmetricity
ftArray = single(zeros(nx,ny,ny));

if mod(ny,2) == 1
    ftArray_half = reshape(FS,nx_half,ny,ny);
    ftArray(1:nx_half+1,:,:) = ftArray_half;
    
    if mod(nx,2) == 1
        ftArray_ahalf = ftArray_half(1:end-1,:,:);
    else
        ftArray_ahalf = ftArray_half(2:end-1,:,:);
    end
    
    ftArray_ahalf = flip(flip(ftArray_ahalf, 2),3);   
    ftArray_ahalf = flip(ftArray_ahalf,1);
    
    ftArray(nx_half+2:end,:,:) = conj(ftArray_ahalf);
else
    ftArray_half_ext = reshape(FS,nx_half+1,ny+1,ny+1);
    ftArray_half = ftArray_half_ext(:,1:end-1,1:end-1);
    ftArray(1:nx_half+1,:,:) = ftArray_half;
    
    if mod(nx,2) == 1
        ftArray_ahalf_ext = ftArray_half_ext(1:end-1,:,:);
    else
        ftArray_ahalf_ext = ftArray_half_ext(2:end-1,:,:);   
    end
    ftArray_ahalf_ext = flip(flip(ftArray_ahalf_ext, 2),3);
    ftArray_ahalf_ext = flip(ftArray_ahalf_ext,1);
    
    ftArray(nx_half+2:end,:,:) = conj(ftArray_ahalf_ext(:,1:end-1,1:end-1));
    
end

end



function D = distancePointsPlane_YY(points, normvec)
    %distancePointsPlane_YY unsigned distances betwen 3D points and a plane
    % through origin
    %
    %   D = distancePointsPlane_YY(point, normvec)
    %   Returns the euclidean distance between points and a plane going through origin with normal vector normvec,
    %   given by: 
    %   points : (3 x n) array of 3D vectors
    %   normvec : (3 x 1) or (1 x 3) array of normal vector
    %   D     : (1 x n) vector  
 
    %
    %   ---------
    %   author : Y. Yang, UCLA Physics and Astronomy 
    %   created the 05/03/2015.
    %

    % normalized plane normal
    normvec = normvec(:) / norm(normvec);
    D = abs(points(1,:)*normvec(1) + points(2,:)*normvec(2) + points(3,:)*normvec(3));

end


function x = closestpoint(n, d, p)
    % n is the vector [A,B,C] that defines the plane
    % d is the distance of the plane from the origin
    % p is the point  [P,Q,R]
    if size(p,2) == 1
        v = (d - sum(p.*n)) / sum(n.*n);
        x = p + v * n;
    else
        nr = repmat(n,[1 size(p,2)]);
        v = (d - sum(p.*nr,1)) / sum(n.*n);
        x = p + repmat(v,[3 1]) .* nr;
    end
end

function dd = MatrixQuaternionRot(vector,theta)

theta = theta*pi/180;
vector = vector/sqrt(dot(vector,vector));
w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
      2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
      2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];

dd = RotM;
end
function out = flip(in,dimension)

switch dimension
    case 1
        for i = 1:size(in,3)
            out(:,:,i) = flipud(in(:,:,i));
        end   
    case 2
        for i = 1:size(in,3)
            out(:,:,i) = fliplr(in(:,:,i));
        end   
    case 3
        for i = 1:size(in,2)
            out(:,i,:) = fliplr(squeeze(in(:,i,:)));
        end   
end
end

