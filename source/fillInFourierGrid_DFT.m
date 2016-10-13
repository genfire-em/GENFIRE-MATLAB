% fillInFourierGrid_DFT
% Y. Yang, UCLA Physics & Astronomy
% First version date: 2015. 04. 30.
% output parameter: rec (nx x ny x nx) IFFT of ftArray
%                   ftArray (nx x ny x nx) interpolated Fourier 3D array
%                   CWArray (nx x ny x nx) interpolated confidence weight
%                   SDArray (nx x ny x nx) interpolated distance
%
%
% input parameter: PROJvol - measured projection, (nx x ny x P) array with P # of projections
%                  angles - 3 Euler angle array, (3 x P) array
%                  K_thresh - threshold for acceptible distance
%                   n1 (length of oversampled projection 1st dim)
%                    n2 (length of oversampled projection 2nd dim) 
%                   FRCcoeffs (Fourier ring correlation result for
%                              confidence weighting)%
%                   CentroSymmetricity - 0 for complex-valued reconstruction
%                                        1 for real-valued reconstruction,
%                                        centrosymmetry will be enforced to
%                                        save calculation time
%                   doGPU - 0 for not using GPU
%                           1 for using GPU
%                   doCTFcorr - 0 for not doing CTF correction
%                                1 for doing CTF correction
%                   CTFparameters - CTF parameters if doing CTF correction
%                   
%
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
%
% Fourth version date: 2016. 4. 11. (YY)
% Change: 1. Cleaned up the code a little bit
%
%         2. Made switch for centrosymmetricity, in case of complex
%         reconstruction
%
%         3. CTF correction
%
%         4. inline C function for speedup
%
% Sixth version date: 2016. 6. 26. (YY)
% Change: 1. Cleaned up the code a little bit
%
%         2. Made siwtches for CTF correction
%
%         3. Wiener filter CTF correction 
%
% Seventh version date: 2016. 8. 11. (YY)
% Change: 1. the output should be nx x ny x nx array (y is the rotation
%              axis)
%         2. fixed bug for correctly determining the cutoff sphere
%
% Eighth version date: 2016. 8. 24. (YY)
% Change: 1. C function disabled because it is actually slower
%         2. Implemented GPU capability
%
% Nineth version date: 2016. 8. 25. (YY)
% Change: 1. Made the function consistent with AJ's convention
             
function [rec, ftArray, CWArray, SDArray] = fillInFourierGrid_DFT(PROJvol, angles, K_thresh, n1, n2, FRCcoeffs, CentroSymmetricity, doGPU, doCTFcorr, CTFparameters)

tic

% original projection dimensions
n1_ori = size(PROJvol,1);
n2_ori = size(PROJvol,2);

% if distance below minInvThresh, minInvThresh will be used
% this is to prevent division by zero
minInvThresh = 0.00001;

% initialize normal vectors and rotation matrices array
normVECs = zeros(size(PROJvol,3),3);
rotMATs = zeros(3,3,size(PROJvol,3));

phis = angles(:,1);
thetas = angles(:,2);
psis = angles(:,3);

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
end


% initiate Fourier space indices
if CentroSymmetricity
    k1 = (-1*floor(n1/2):1:0) ;
    n_k1 = floor(n1/2);    
    
    k2 = (-1*floor(n2/2):1:floor(n2/2)) ;  
    k3 = (-1*floor(n1/2):1:floor(n1/2)) ;    
    
else
    k1 = (-1*ceil((n1-1)/2):1:floor((n1-1)/2)) / n1;
    
    k2 = (-1*ceil((n2-1)/2):1:floor((n2-1)/2)) / n2;
    
    k3 = k1;
end


% Fourier grid
% in case of centrosymmetry, only half of k1 will be interpolated
% and centrosymmetricity will be enforced later
[~, K1, ~] = meshgrid(k2,k1,k3);

% initialize variables
FS = zeros(size(K1))*(-1); % Fourier points
CW = zeros(size(K1))*(-1); % Fourier points

if nargout>3
    SD = zeros(size(K1))*(-1); % Fourier points
end


Numpt = zeros(size(K1)); % array to store how many points found per Fourier point
invSumTotWeight = zeros(size(K1)); % array to store sum of weights of Fourier point

% initiate Fourier space indices

k1_ori = (-1*ceil((n1_ori-1)/2):1:floor((n1_ori-1)/2)) ;    
k2_ori = (-1*ceil((n2_ori-1)/2):1:floor((n2_ori-1)/2)) ;    

[K20, K10] = meshgrid(k2_ori,k1_ori);

if doGPU
    K10G = gpuArray(K10(:));
    K20G = gpuArray(K20(:));
end
    
for p=1:size(PROJvol,3)
    % current projection
    curr_proj = squeeze(PROJvol(:,:,p));
        
    [K2, K1, K3] = meshgrid(k2,k1,k3);
    
    % obtain points-to-plane distance
    D = distancePointsPlane_YY([K1(:) K2(:) K3(:)]', normVECs(p,:));
    
    % find Fourier points within the threshold
    Dind = find(D < K_thresh);
    
    %tic
    % rotate the plane to zero degree
    CP = closestpoint(normVECs(p,:)',0,[K1(Dind); K2(Dind); K3(Dind)]);

    %toc
    
    CP_plane = (squeeze(rotMATs(:,:,p)))\CP;
    
    clear KX KY KZ
    % picked closest point which is within the projection plain, x coordinate must be zero after rotation
    if sum(abs(CP_plane(3,:)) > 0.0001) > 0
        fprintf(1,'something wrong!\n');                
    end
    
    % consider Fourier points only within the resolution circle
    Gind = Dind(abs(CP_plane(1,:)) <= n1/2 & abs(CP_plane(2,:)) <= n2/2);  % good indices  
    G_CP_plane = CP_plane(:,abs(CP_plane(1,:)) <= n1/2 & abs(CP_plane(2,:)) <= n2/2 );  % good in-plane coordinates
    
    %determine the available memory in MB
    if doGPU
        GPUinfo = gpuDevice();
        av_memory_size = round(GPUinfo.AvailableMemory/1000000);
    else
        
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
    
    %if Cscript
    %    FS = DFTinterpolate_C(K10,K20,G_CP_plane,n1,n2,k2,k1_half_1_ori,Mindist,D,Gind);
    %tic
    % loop over Fourier points within the threshold
    for i=1:cutloopnum 
        curr_indices = ((i-1)*block_size+1):(i*block_size);
        
        if i > cutnum
            curr_indices = (cutnum*block_size+1):length(Gind);
        end
         
        % CTF correction
        if doCTFcorr        
            [CTFF, gamma, ED, EA] = ctf_calc_DFT_YY(G_CP_plane(:,curr_indices),CTFparameters(p).Voltage,CTFparameters(p).pix_size,n1,n2,CTFparameters(p).Defocus,CTFparameters(p).DefStdev,CTFparameters(p).Alpha);%get CTF
            %[CTF, gamma, EnvD,EnvA] = ctf_calc_YY(Voltage,pix_size,Isize,Defocus,DefStdev,Alpha,ignore_first_peak)
            
            CTF = CTFF.*ED.*EA;           
           
                
            % in case of Wiener filtering
            if CTFparameters(p).WF
                
                corrCoeffs = CTFparameters(p).corrCoeffs;
                FSCbinEdges = CTFparameters(p).FSCbinEdges;
                
                % histogram count, F_mag < 0 will have bin index 1, and 
                % F_mag > last bin edge (this should be sqrt(2)*halfsize)
                % will have the last index
                [~, Fbin] = histc(sqrt(sum(G_CP_plane(:,curr_indices).^2,1)),[-inf FSCbinEdges inf]);
                
                if sum(Fbin==1) > 0 
                    disp('FSC bin out of bound (negative K_mag), something wrong!')
                else
                    corrCoeffs(end+1) = 0; % for Fourier magnitude outside sqrt(2)*halfsize
                                           % (which hopefully does not
                                           % exist), assign zero FSC.
                    
                    % adjust index (1st index is negative F_mag, which should not exist)                       
                    Fbin = Fbin - 1;
                    
                    % get FSC values from the histogram count result
                    FSCresults = corrCoeffs(Fbin);
                    if sum(FSCresults<0) > 0
                        fprintf(1,'negative # of FSC: %d\n',sum(FSCresults<0));
                    end
                    FSCresults(FSCresults<0) = 0;
                    
                    % get SSNR from FSC
                    SSNR = 1./(1-FSCresults);
                end                    

                throwOutThreshhold=CTFparameters(p).throwOutThreshhold;
                CTFcorr = ones(1,length(curr_indices));
                
                % using Threshold, throw out datapoints whose CTF is below
                % the threshold
                CTFcorr(abs(CTF)>=throwOutThreshhold) = CTF(abs(CTF)>=throwOutThreshhold);
                CTFcorr(abs(CTF)<throwOutThreshhold) = [];
                curr_indices(abs(CTF)<throwOutThreshhold) = [];
                SSNR(abs(CTF)<throwOutThreshhold) = [];
                
                % Wiener filter
                CTFcorr = conj(CTFcorr) ./ (abs(CTFcorr).^2 + 1./SSNR);
             
            % in case of only phase flipping
            elseif CTFparameters(p).onlyPF
                CTFcorr = ones(1,length(curr_indices));
                CTFcorr(CTF<0) = -1;
                
            % direct correction with Thresholding
            else                
                throwOutThreshhold=CTFparameters(p).throwOutThreshhold;
                CTFcorr = ones(1,length(curr_indices));
                CTFcorr(abs(CTF)>=throwOutThreshhold) = 1./CTF(abs(CTF)>=throwOutThreshhold);
                CTFcorr(abs(CTF)<throwOutThreshhold) = [];
                curr_indices(abs(CTF)<throwOutThreshhold) = [];
            end
            
        % no CTF correction
        else
            CTFcorr = ones(1,length(curr_indices));
        end               
        
        if doGPU
            G_CP_plane1_GPU_n = gpuArray(G_CP_plane(1,curr_indices)/n1);
            G_CP_plane2_GPU_n = gpuArray(G_CP_plane(2,curr_indices)/n2);
            curr_proj_GPU = gpuArray(curr_proj(:));
        
        
            % DFT calculation
            FpointsG = sum(bsxfun(@times, curr_proj_GPU, exp(-1*1i*2*pi*(K10G*G_CP_plane1_GPU_n+K20G*G_CP_plane2_GPU_n))),1);

            Fpoints = gather(FpointsG);     
            Fpoints = CTFcorr.*Fpoints;
            
            clear G_CP_plane1_GPU_n G_CP_plane2_GPU_n curr_proj_GPU FpointsG     
        else
            Fpoints = CTFcorr.*sum(bsxfun(@times, curr_proj(:), exp(-1*1i*2*pi*(K10(:)*G_CP_plane(1,curr_indices)/n1+K20(:)*G_CP_plane(2,curr_indices)/n2))),1);
        end
        
        
        %weighted avearaging
        CIND = Gind(curr_indices);

        currDist = D(CIND);
        
        Q = sqrt( (G_CP_plane(1,curr_indices)/n1).^2 + (G_CP_plane(2,curr_indices)/n2).^2 );
        
        currentConfidenceWeights = zeros(size(Q));
        currentConfidenceWeights(:) = interp1(linspace(0,1+1e-10,size(FRCcoeffs,2)),FRCcoeffs(p,:),Q(:),'linear');%make weighting map from average FRC
        currentConfidenceWeights(isnan(currentConfidenceWeights))=0;%extrapolated values will be NaN and are outside the resolution circle anyway

        currDist(currDist < minInvThresh) = minInvThresh; % if distance smaller than minInvThresh, put minInvThresh (to prevent divison by zero)
        currInvDist = 1./ currDist;           % inverse distance
        currTotWeight = currentConfidenceWeights.*currInvDist;

        % re-average inverse distance
        FS(CIND) = FS(CIND).* invSumTotWeight(CIND) + currTotWeight.*Fpoints;
        CW(CIND) = CW(CIND).* invSumTotWeight(CIND) + currTotWeight.*currentConfidenceWeights;
        invSumTotWeight(CIND) = invSumTotWeight(CIND) + currTotWeight;
        FS(CIND) = FS(CIND) ./ invSumTotWeight(CIND);
        CW(CIND) = CW(CIND) ./ invSumTotWeight(CIND);
        
        if nargout>3
            SD(CIND) = SD(CIND).* invSumTotWeight(CIND) + currTotWeight.*currDist;
            SD(CIND) = SD(CIND) ./ invSumTotWeight(CIND);
        end
        Numpt(CIND) = Numpt(CIND) + 1;     
        
        clear Fpoints
    end
    
    
    
    %toc
end

clear Mindist invSumdist Numpt

% enforce centrosymmetricity
%CentroSymmetricity
if CentroSymmetricity
    ftArray = single(zeros(n1,n2,n1));
    CWArray = single(zeros(n1,n2,n1));
    
    if nargout > 3
        SDArray = single(zeros(n1,n2,n1));
    end

    if mod(n2,2) == 1
        

        if mod(n1,2) == 1
            ftArray_half = reshape(FS,n_k1+1,n2,n1);        
            ftArray(1:n_k1+1,:,:) = ftArray_half;
            
            ftArray_ahalf = ftArray_half(1:end-1,:,:);            
            ftArray_ahalf = flip(flip(ftArray_ahalf, 2),3);   
            ftArray_ahalf = flip(ftArray_ahalf,1);
            ftArray(n_k1+2:end,:,:) = conj(ftArray_ahalf);
            
            clear ftArray_half ftArray_ahalf

            CWArray_half = reshape(CW,n_k1+1,n2,n1);        
            CWArray(1:n_k1+1,:,:) = CWArray_half;
            
            CWArray_ahalf = CWArray_half(1:end-1,:,:);            
            CWArray_ahalf = flip(flip(CWArray_ahalf, 2),3);   
            CWArray_ahalf = flip(CWArray_ahalf,1);
            CWArray(n_k1+2:end,:,:) = conj(CWArray_ahalf);
            
            clear CWArray_half CWArray_ahalf
            
            if nargout > 3
                SDArray_half = reshape(SD,n_k1+1,n2,n1);        
                SDArray(1:n_k1+1,:,:) = SDArray_half;

                SDArray_ahalf = SDArray_half(1:end-1,:,:);            
                SDArray_ahalf = flip(flip(SDArray_ahalf, 2),3);   
                SDArray_ahalf = flip(SDArray_ahalf,1);
                SDArray(n_k1+2:end,:,:) = conj(SDArray_ahalf);
                
                clear CWArray_half CWArray_ahalf           
            end            
        
        else
            ftArray_half_ext = reshape(FS,n_k1+1,n2,n1+1);   
            ftArray_half = ftArray_half_ext(:,:,1:end-1);
            ftArray(1:n_k1+1,:,:) = ftArray_half;
            
            ftArray_ahalf_ext = ftArray_half_ext(2:end-1,:,:);            
            ftArray_ahalf_ext = flip(flip(ftArray_ahalf_ext, 2),3);   
            ftArray_ahalf_ext = flip(ftArray_ahalf_ext,1);
            ftArray(n_k1+2:end,:,:) = conj(ftArray_ahalf_ext(:,:,1:end-1));
            
            clear ftArray_half ftArray_half_ext ftArray_ahalf_ext

            
            CWArray_half_ext = reshape(CW,n_k1+1,n2,n1+1);   
            CWArray_half = CWArray_half_ext(:,:,1:end-1);
            CWArray(1:n_k1+1,:,:) = CWArray_half;
            
            CWArray_ahalf_ext = CWArray_half_ext(2:end-1,:,:);            
            CWArray_ahalf_ext = flip(flip(CWArray_ahalf_ext, 2),3);   
            CWArray_ahalf_ext = flip(CWArray_ahalf_ext,1);
            CWArray(n_k1+2:end,:,:) = conj(CWArray_ahalf_ext(:,:,1:end-1));
            
            clear CWArray_half CWArray_half_ext CWArray_ahalf_ext
            
            if nargout>3
                SDArray_half_ext = reshape(SD,n_k1+1,n2,n1+1);   
                SDArray_half = SDArray_half_ext(:,:,1:end-1);
                SDArray(1:n_k1+1,:,:) = SDArray_half;

                SDArray_ahalf_ext = SDArray_half_ext(2:end-1,:,:);            
                SDArray_ahalf_ext = flip(flip(SDArray_ahalf_ext, 2),3);   
                SDArray_ahalf_ext = flip(SDArray_ahalf_ext,1);
                SDArray(n_k1+2:end,:,:) = conj(SDArray_ahalf_ext(:,:,1:end-1));

                clear SDArray_half SDArray_half_ext SDArray_ahalf_ext
            end
        end
        
        

        
    else


        if mod(n1,2) == 1
            ftArray_half_ext = reshape(FS,n_k1+1,n2+1,n1);
            ftArray_half = ftArray_half_ext(:,1:end-1,:);
            ftArray(1:n_k1+1,:,:) = ftArray_half;
        
            ftArray_ahalf_ext = ftArray_half_ext(1:end-1,:,:);
            ftArray_ahalf_ext = flip(flip(ftArray_ahalf_ext, 2),3);
            ftArray_ahalf_ext = flip(ftArray_ahalf_ext,1);

            ftArray(n_k1+2:end,:,:) = conj(ftArray_ahalf_ext(:,1:end-1,:));
            
            clear ftArray_half ftArray_half_ext ftArray_ahalf_ext
            
            CWArray_half_ext = reshape(CW,n_k1+1,n2+1,n1);
            CWArray_half = CWArray_half_ext(:,1:end-1,:);
            CWArray(1:n_k1+1,:,:) = CWArray_half;
        
            CWArray_ahalf_ext = CWArray_half_ext(1:end-1,:,:);
            CWArray_ahalf_ext = flip(flip(CWArray_ahalf_ext, 2),3);
            CWArray_ahalf_ext = flip(CWArray_ahalf_ext,1);

            CWArray(n_k1+2:end,:,:) = conj(CWArray_ahalf_ext(:,1:end-1,:));           
            
            clear CWArray_half CWArray_half_ext CWArray_ahalf_ext
            
            if nargout>3

                SDArray_half_ext = reshape(SD,n_k1+1,n2+1,n1);
                SDArray_half = SDArray_half_ext(:,1:end-1,:);
                SDArray(1:n_k1+1,:,:) = SDArray_half;

                SDArray_ahalf_ext = SDArray_half_ext(1:end-1,:,:);
                SDArray_ahalf_ext = flip(flip(SDArray_ahalf_ext, 2),3);
                SDArray_ahalf_ext = flip(SDArray_ahalf_ext,1);

                SDArray(n_k1+2:end,:,:) = conj(SDArray_ahalf_ext(:,1:end-1,:));           

                clear SDArray_half SDArray_half_ext SDArray_ahalf_ext
            end

        else
            ftArray_half_ext = reshape(FS,n_k1+1,n2+1,n1+1);
            ftArray_half = ftArray_half_ext(:,1:end-1,1:end-1);
            ftArray(1:n_k1+1,:,:) = ftArray_half;
            
            ftArray_ahalf_ext = ftArray_half_ext(2:end-1,:,:);   
            ftArray_ahalf_ext = flip(flip(ftArray_ahalf_ext, 2),3);
            ftArray_ahalf_ext = flip(ftArray_ahalf_ext,1);

            ftArray(n_k1+2:end,:,:) = conj(ftArray_ahalf_ext(:,1:end-1,1:end-1));
            
            clear ftArray_half ftArray_half_ext ftArray_ahalf_ext
            
            CWArray_half_ext = reshape(CW,n_k1+1,n2+1,n1+1);
            CWArray_half = CWArray_half_ext(:,1:end-1,1:end-1);
            CWArray(1:n_k1+1,:,:) = CWArray_half;
            
            CWArray_ahalf_ext = CWArray_half_ext(2:end-1,:,:);   
            CWArray_ahalf_ext = flip(flip(CWArray_ahalf_ext, 2),3);
            CWArray_ahalf_ext = flip(CWArray_ahalf_ext,1);

            CWArray(n_k1+2:end,:,:) = conj(CWArray_ahalf_ext(:,1:end-1,1:end-1));
            
            clear CWArray_half CWArray_half_ext CWArray_ahalf_ext
            
            if nargout>3
                SDArray_half_ext = reshape(SD,n_k1+1,n2+1,n1+1);
                SDArray_half = SDArray_half_ext(:,1:end-1,1:end-1);
                SDArray(1:n_k1+1,:,:) = SDArray_half;

                SDArray_ahalf_ext = SDArray_half_ext(2:end-1,:,:);   
                SDArray_ahalf_ext = flip(flip(SDArray_ahalf_ext, 2),3);
                SDArray_ahalf_ext = flip(SDArray_ahalf_ext,1);

                SDArray(n_k1+2:end,:,:) = conj(SDArray_ahalf_ext(:,1:end-1,1:end-1));

                clear SDArray_half SDArray_half_ext SDArray_ahalf_ext
            end
            
        end


    end
else
    ftArray = reshape(FS,n1,n2,n1);
    CWArray = reshape(CW,n1,n2,n1);
    if nargout> 3
        SDArray = reshape(SD,n1,n2,n1);
    end
end

rec = real(my_ifft(ftArray));

timeTakenToFillInGrid = toc;
timeTakenToFillInGrid = round(10*timeTakenToFillInGrid)./10;
fprintf('GENFIRE: Fourier grid assembled in %.12g seconds.\n\n',timeTakenToFillInGrid);

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

