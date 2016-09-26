%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                         %%
%%                        Welcome to GENFIRE!                              %%
%%           GENeralized Fourier Iterative REconstruction                  %%
%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENFIRE is a robust, Fourier-based reconstruction algorithm that is
%% capable of using a limited set of input projections to generate a 3D reconstruction
%% while also partially retrieving missing projection information. It does this by iterating 
%% between real and reciprocal space and applying simple constraints in each to find an optimal solution that  
%% correlates with the input projections while simultaneously obeying real space conditions
%% such as positivity (the assumption that there is no negative electron density)
%% and support (we require that the 3D object resulting from zero padded projections
%% exists entirely within some smaller region). The result is a more consistent 
%% and faithful reconstruction with superior contrast and, in many cases, resolution when
%% compared with more traditional 3D reconstruction algorithms such as
%% Filtered Back-Projection.  

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectionBottom= 1;
projectionTopIncrement = 49;


addpath ./source/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   User Parameters   %%%
filename_Projections = './data/tomoProjections.mat';%%filename of projections, which should be size NxNxN_projections where N_projections is the number of projections

filename_Angles = './data/tomoAngles.mat';%%angles can be either a 1xN_projections array containing a single tilt series, or

%%a 3xN_projections array containing 3 Euler angles for each projections in the form [phi;theta;psi]
filename_Support = './data/support60.mat'; %% NxNxN binary array specifying a region of 1's in which the reconstruction can exist 

% filename_InitialModel = '.\models\betagal\model.mat';

filename_results = './results/GENFIRE_rec_allBotHalf.mat';


global numIterations 
numIterations = 50; 

global pixelSize
pixelSize = .5; 

oversamplingRatioX =3; %%The code will zero-pad projections for you to the inputted oversampling ratio. If your projections are already oversampled
%%then set this to 1.
oversamplingRatioY =1.0; %%The code will zero-pad projections for you to the inputted oversampling ratio. If your projections are already oversampled
%%then set this to 1.

griddingMethod = 1; %% 1) Use fastest FFT method with mex-compiled function weightVals.cpp (see INSTALL_NOTES.txt for more info). 2) Use FFT method. 3) Use DFT method, which exactly calculates the closest measured value to each grid point rather than using the nearest FFT pixel. This is the most accurate but also slowest method

interpolationCutoffDistance =.7; %%radius of sphere (in pixels) within which to include measured datapoints 
%%when assembling the 3D Fourier grid

% particleWindowSize = 46; %size of window to crop out of projections (before oversampling). Comment out to just use the size of the input projections This is especially useful if you are doing CTF correction, as the CTF correction
%% is more accurate on larger arrays, but the portion of the micrograph containing your particle may be smaller than the input images. This allows you to CTF correct
%% the larger projections, then crop a window out and pad that with zeros for reconstruction.

ComputeFourierShellCorrelation =0; %%set to 1 to divide dataset in half, independently reconstruct, and compute Fourier Shell Correlation (FSC) between two halves.

numBins = 50; %number of bins for FRC averaging
%% If you do not need FRC, set ComputeFourierShellCorrelation to 0 for speed as the FRC calculation requires reconstructing everything twice




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%     Parameters that are unlikely to need to changing        %%%
targetCorrelation = 0.75;%correlation coefficient value above which is considered "strongly correlated" when computing backprojected FRC to input projections. 
Rfree_complex_target = 0.55;%value between 0 and 1 indicating the maximum acceptable convergence value for Rfree_complex for complex constraint cutoff (I recommend values between 0.5 to 0.65)
percentValuesForRfree = 0.05;
phaseErrorSigmaTolerance = 75*pi/180;%the standard deviation of the phases of the constituent datapoints is calculated for each Fourier grid point. If this standard deviation is above this tolerance threshhold, that grid point is ignored, as the phase is uncertain. This generally happens at the boundary of speckles, where the phase changes abruptly.

runAutoResolutionCutoffFinder =1;%%set to 1 to automatically find a suitable resolution cutoff for enforceable data
%to override the auto-search and set your own cutoffs for resolution to enforce, set runAutoResolutionCutoffFinder = 0 then input
%your desired cutoff here
complexCutoff = 1;

doCTFcorrection = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%   Begin Reconstruction   %%%
clc
angles = single(importdata(filename_Angles));%casting to single saves memory
projections = single(importdata(filename_Projections));
projections = projections(:,projectionBottom:projectionBottom+projectionTopIncrement,:);

if exist('filename_InitialModel','var')
    initialObject = single(importdata(filename_InitialModel));
else
    initialObject = [];
end


global support %make support variable globally accessable to avoid passing copies of large arrays around to different functions
support = single(importdata(filename_Support));
support = support(:,projectionBottom:projectionBottom+projectionTopIncrement,:);

% if ~exist('particleWindowSize','var')
%     particleWindowSizeY = size(support,2);
%     particleWindowSizeX = size(support,1);
% else
%     if particleWindowSize ~= size(support,2)
%         error('GENFIRE: ERROR! The size of your projections (set by particleWindowSize) does not match the dimensions of your support!')
%     end
% end
particleWindowSizeY = size(support,2);
particleWindowSizeX = size(support,1);
n2Y = particleWindowSizeY/2;%radius of smaller array
newDimY = particleWindowSizeY*oversamplingRatioY;%size of oversampled array
paddingY = round((newDimY-particleWindowSizeY)/2);%how many zeros to add

n2X = particleWindowSizeX/2;%radius of smaller array
newDimX = particleWindowSizeX*oversamplingRatioX;%size of oversampled array
paddingX = round((newDimX-particleWindowSizeX)/2);%how many zeros to add

averageStrongCorrelationCutoff = complexCutoff; %average resolution across all projections where backprojection FRC drops below targetCorrelation

%zero pad projections to user-inputted oversampling ratio
numProj = size(projections,3);

support = padarray(support,[paddingX paddingY paddingX]); %%zero pad support
if exist('filename_InitialModel','var')
    initialObject =  padarray(initialObject,[paddingX paddingY paddingX]); 
end

if griddingMethod==3
   error('GENFIRE: DFT gridding currently not supported in v1.5, update in near future. Change griddingMethod to 1 or 2.') 
end
%If there is only one tilt angle (such as in tomography), fill in other
%Euler angles with zeros
if size(angles,1)>3
    error('The dimension of the angles is incorrect.')
end
if size(angles,1) ==1 
angles = [zeros(1,length(angles));angles;zeros(1,length(angles))];%tomography tilt is the theta angle
end

%construct K-space indices
if mod(size(support,1),2)==0
ncK1 = size(support,1)/2+1;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1-1)./n2K1;
elseif size(support,1)==1
vec1 = 0;
ncK1 = 1;
n2K1 = 0;
else
ncK1 = (size(support,1)+1)/2;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1)./n2K1; 
end

if  mod(size(support,2),2)==0
ncK2 = size(support,2)/2+1;%central pixel
n2K2 = ncK2-1;%max radius
vec2 = (-n2K2:n2K2-1)./n2K2;
elseif size(support,2)==1
vec2 = 0;
ncK2 = 1;
n2K2 = 0;
else
ncK2 = (size(support,2)+1)/2;%central pixel
n2K2 = ncK2-1;%max radius
vec2 = (-n2K2:n2K2)./n2K2; 
end

if  mod(size(support,3),2)==0
ncK3 = size(support,3)/2+1;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3-1)./n2K3;
elseif size(support,3)==1
vec3 = 0;
ncK3 = 1;
n2K3 = 0;
else
ncK3 = (size(support,3)+1)/2;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3)./n2K3; 
end
[Kx Ky Kz] = meshgrid(vec2,vec1,vec3);%grid of Fourier indices
Kmags = sqrt(Kx.^2+Ky.^2+Kz.^2);%take magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if runAutoResolutionCutoffFinder 
    foundGoodWeights = 0;%flag to indicate when we are done
    averageStrongCorrelationCutoff = 1;%initialize to 100% for search
else
    fprintf('GENFIRE: Automatic resolution finder overwritten, continuing with user-inputted complex enforcement cutoff at %.12g%%\n\n',100*averageStrongCorrelationCutoff)
    foundGoodWeights = 1;
    
end
        


resRange = -0.05;%thickness of resolution ring to use for removal of datapoints for Rfree test
reconstructionIdentifierNumber = 0;%reconstructionIdentifierNumbers number of times the resoluton cutoff finding loop has run
if ComputeFourierShellCorrelation
    fprintf('GENFIRE: Dividing dataset in half...\n\n')
    pjCutoff1 = round((numProj+1)/2);%divide dataset in half
    pjCutoff2 = numProj;
    fprintf('GENFIRE: Assembling Fourier grids...\n\n');
       
    %%INTERPOLATE TO GRID
    %%fillInFourierGrid_tomo_C contains the c++ code "weightVals.cpp" that must be
    %%compiled on your local machine to run. If you cannot get it to work with
    %%the command "mex weightVals.cpp" then use the version fillInFourierGrid
    %%which does the same thing but slower.
    initWeights = ones(size(projections,3),numBins);
    switch griddingMethod
        case 1
        [rec1 measuredK1] = fillInFourierGrid_tomo_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),oversamplingRatioX,oversamplingRatioY,interpolationCutoffDistance,initWeights,doCTFcorrection);%interpolate projections to Fourier grid
        [rec2 measuredK2] = fillInFourierGrid_tomo_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),particleWindowSizeY,oversamplingRatio,interpolationCutoffDistance,initWeights,doCTFcorrection);
       
%         [rec1 measuredK1] = fillInFourierGrid_tomo_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),particleWindowSize,oversamplingRatio,interpolationCutoffDistance,initWeights,doCTFcorrection);%interpolate projections to Fourier grid
%         [rec2 measuredK2] = fillInFourierGrid_tomo_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),particleWindowSize,oversamplingRatio,interpolationCutoffDistance,initWeights,doCTFcorrection);
        case 2
        [rec1 measuredK1] = fillInFourierGrid(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),particleWindowSize,oversamplingRatio,interpolationCutoffDistance,initWeights,doCTFcorrection);%interpolate projections to Fourier grid
        [rec2 measuredK2] = fillInFourierGrid(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),particleWindowSize,oversamplingRatio,interpolationCutoffDistance,initWeights,doCTFcorrection);
        case 3
        ori_projections = single(importdata(filename_Projections));%load unpadded projections (function format requires this)
        Typeind = 2;
        tic
        measuredK1 = My_fill_cubicgrid_ver3_1(dim, dim, ori_projections(:,:,1:pjCutoff1), angles(:,1:pjCutoff1), interpolationCutoffDistance, Typeind, newDim, newDim,1);
        measuredK2 = My_fill_cubicgrid_ver3_1(dim, dim, ori_projections(:,:,(pjCutoff1+1):pjCutoff2), angles(:,(pjCutoff1+1):pjCutoff2), interpolationCutoffDistance, Typeind, newDim, newDim,1);
        interpTime = toc
        rec1 = real(my_ifft(measuredK1));
        rec2 = real(my_ifft(measuredK2));    
    end

    if runAutoResolutionCutoffFinder %automatically find suitable resolution limit to enforce
        foundGoodWeights = 0;
        averageStrongCorrelationCutoff = 1;
    else
        foundGoodWeights = 1;
    end
    
    resolutionWeights1 = zeros(size(Kmags)); %% first reconstruction uses resolution extension/suppression, where lower resolution information is enforced
    %%initially and the maximum enforced resolution increases. This is followed by resolution where suppression, which is the same process run backwards, so at the final iteration 
    %%only the lowest resolution is being enforced again.
    resolutionWeights1(measuredK1~=0) = 1-Kmags(measuredK1~=0);%make lower resolution have higher confidence. Higher confidence means enforced earlier in the reconstruction and for longer overall than low confidence
    constraintConfidenceWeights1 = resolutionWeights1;
    
    resolutionWeights2 = zeros(size(Kmags)); %% first reconstruction uses resolution extension/suppression, where lower resolution information is enforced
    %%initially and the maximum enforced resolution increases. This is followed by resolution where suppression, which is the same process run backwards, so at the final iteration 
    %%only the lowest resolution is being enforced again.
    resolutionWeights2(measuredK2~=0) = 1-Kmags(measuredK2~=0);%make lower resolution have higher confidence. Higher confidence means enforced earlier in the reconstruction and for longer overall than low confidence
    constraintConfidenceWeights2 = resolutionWeights2;
    clear resolutionWeights;
   
    constraintEnforcementDelayWeights = [.95:-.05:-.15  -.15:.05:.95];%This means the reconstruction will start by only enforcing data above confidence 0.95, which corresponds to the lowest 5% of resolution based on how we just constructed the resolutionWeights. 
    
    
    while ~foundGoodWeights
        if reconstructionIdentifierNumber==0
            %%remove the values for Rfree
            measuredPointInd_complex = find(measuredK1~=0&Kmags>=(averageStrongCorrelationCutoff+resRange)&Kmags<averageStrongCorrelationCutoff); %candidate values for Rfree_complex
            P = randperm(numel(measuredPointInd_complex)); %shuffle values
            measuredPointInd_complex = measuredPointInd_complex(P);
            cutoffInd_complex = floor(numel(measuredPointInd_complex).*percentValuesForRfree); %take 5% of values
            R_freeInd_complex = measuredPointInd_complex(1:cutoffInd_complex);

            %perform reconstruction
            fprintf('GENFIRE: Reconstructing... \n\n');
            if isempty(initialObject) %if no initial object was provided, just use the IFFT of the constraints as initial object
                [recHalf1, errK1,Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,measuredK1,resolutionWeights1,constraintEnforcementDelayWeights);
                [recHalf2, errK2] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,measuredK2,resolutionWeights2,constraintEnforcementDelayWeights);
            else
                [recHalf1, errK1,Rfree_complex] = GENFIRE_iterate(numIterations,initialObject,support,measuredK1,resolutionWeights1,constraintEnforcementDelayWeights);
                [recHalf2, errK2] = GENFIRE_iterate(numIterations,initialObject,support,measuredK2,resolutionWeights2,constraintEnforcementDelayWeights); 
            end
            
        else
            
        end
        
        GENFIRE_recSmall1 = recHalf1(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);%extract the original sized result from the oversampled reconstruction
        GENFIRE_recSmall2 = recHalf2(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);%extract the original sized result from the oversampled reconstruction

%             GENFIRE_recSmall1 = recHalf1(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction
%             GENFIRE_recSmall2 = recHalf2(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction

%             paddingForBackprojection = round(max(size(GENFIRE_recSmall))/2);%padding to prepare for backprojection calculation. If there is not enough 0 padding the projections are extremely inaccurate
%             GENFIRE_recSmall1 = padarray(GENFIRE_recSmall1,[paddingForBackprojection paddingForBackprojection paddingForBackprojection]);
%             GENFIRE_recSmall2 = padarray(GENFIRE_recSmall2,[paddingForBackprojection paddingForBackprojection paddingForBackprojection]);

%             unpaddedProjections = single(importdata(filename_Projections));
%             unpaddedProjections = unpaddedProjections(:,:,:);
%             centralPixel = size(projections,2)/2+1;
%             halfWindowSize = particleWindowSize/2;
        %backproject to same angles used for reconstruction and compare to
        %input projections
        if reconstructionIdentifierNumber==0%only do backprojection comparison once
            doneFlag = 1;
            break
            fprintf('GENFIRE: Computing FRC between backprojections of reconstruction and input projections...\n\n')
            [meanFRCcoeff, spatialFrequency,FRCcoeffs1] = backprojectAndCompareFRC(padarray(unpaddedProjections(:,:,1:pjCutoff1),[paddingX paddingY 0]),angles(:,1:pjCutoff1),GENFIRE_recSmall1,numBins,pixelSize,averageStrongCorrelationCutoff,targetCorrelation);
            [meanFRCcoeff, spatialFrequency,FRCcoeffs2] = backprojectAndCompareFRC(padarray(unpaddedProjections(:,:,(pjCutoff1+1):pjCutoff2),[paddingX paddingY 0]),angles(:,(pjCutoff1+1):pjCutoff2),GENFIRE_recSmall2,numBins,pixelSize,averageStrongCorrelationCutoff,targetCorrelation);
%             constraintEnforcementDelayWeights =[.95:-.1:.15 .15 .15:.1:.95]; %%These new weights correspond to enforcing based on FRC correlation coefficients

        end
        convergenceStatus = testConvergence(Rfree_complex);%check that the Rfree converged 

            %if individual cutoffs are being used for each projection, then
            %convergence is found when the average change between cutoffs
            %is very small. The condition for Rfree_complex is the same as
            %in the first case
            if reconstructionIdentifierNumber==1;% reconstructionIdentifierNumber = 0 is the first step where resolution extension is used, reconstructionIdentifierNumber = 1 is the correlation weight step (the final step)
                foundGoodWeights = 1;%ready to exit the loop
            doneFlag = 1;    
                break
%                 averageStrongCorrelationCutoff = max([newAverageStrongCorrelationCutoff averageStrongCorrelationCutoff]);%take highest cutoff, this is just an average to report to the user. Individual weights are used projection by projection for actual reconstruction
%                 fprintf('GENFIRE: Found good constraint weights with average highly-correlated resolution out to %4.4g Angstrom (%.12g%% of maximum resolution)\n\n',1/(round(100*(mean([averageStrongCorrelationCutoff])))./2./pixelSize./100),round(10000*mean([averageStrongCorrelationCutoff]))./100)
            else
                reconstructionIdentifierNumber = reconstructionIdentifierNumber+1;%track how many times the cutoff search has run
% 
%                 if reconstructionIdentifierNumber>10
%                     error('GENFIRE: Cutoff search has stagnated; exiting.\n\n')
%                 end
                fprintf('GENFIRE: Assembling newest Fourier grid...\n\n');
%                 projectionCutoffValues1 = newProjectionCutoffValues1;%update resolution cutoffs
%                 projectionCutoffValues2 = newProjectionCutoffValues2;
%                 averageStrongCorrelationCutoff = mean([projectionCutoffValues1 projectionCutoffValues2]);
                switch griddingMethod
                    case 1
                        [rec1, measuredK1,constraintConfidenceWeights1] = fillInFourierGrid_tomo_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),oversamplingRatioX,oversamplingRatioY,interpolationCutoffDistance,FRCcoeffs1);%reassemble grid with final cutoff values
                        [rec2, measuredK2,constraintConfidenceWeights2] = fillInFourierGrid_tomo_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),oversamplingRatioX,oversamplingRatioY,interpolationCutoffDistance,FRCcoeffs2);
                   
%                         [rec1, measuredK1,constraintConfidenceWeights1] = fillInFourierGrid_tomo_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,FRCcoeffs1);%reassemble grid with final cutoff values
%                         [rec2, measuredK2,constraintConfidenceWeights2] = fillInFourierGrid_tomo_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,FRCcoeffs2);
                    case 2
                        [rec1 measuredK1,constraintConfidenceWeights1] = fillInFourierGrid(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,FRCcoeffs1);%reassemble grid with final cutoff values
                        [rec2 measuredK2,constraintConfidenceWeights2] = fillInFourierGrid(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,FRCcoeffs2);
                    case 3
                        ori_projections = single(importdata(filename_Projections));%load unpadded projections (function format requires this)
                        Typeind = 2;
                        tic
                        measuredK1 = My_fill_cubicgrid_ver3_1(dim, dim, ori_projections(:,:,1:pjCutoff1), angles(:,1:pjCutoff1), interpolationCutoffDistance, Typeind, newDim, newDim,1);%reassemble grid with final cutoff values
                        measuredK2 = My_fill_cubicgrid_ver3_1(dim, dim, ori_projections(:,:,(pjCutoff1+1):pjCutoff2), angles(:,(pjCutoff1+1):pjCutoff2), interpolationCutoffDistance, Typeind, newDim, newDim,1);
                        [null1 null2,constraintConfidenceWeights1] = fillInFourierGrid_tomo_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,initWeights);%currently the DFT code does not compute the confidence weights so still have to use this function for that purpose
                        [null1 null2,constraintConfidenceWeights2] = fillInFourierGrid_tomo_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,initWeights);
                        interpTime = toc
                        rec1 = real(my_ifft(measuredK1));
                        rec2 = real(my_ifft(measuredK2));
                end

            end
    end
    if ~runAutoResolutionCutoffFinder%if auto resolution searching is not turned on, then do the half reconstructions now
        constraintInd_complex1 = find(measuredK1~=0&Kmags<=averageStrongCorrelationCutoff);%indices of measured Fourier points to enforce with complex values
        constraintInd_complex2 = find(measuredK2~=0&Kmags<=averageStrongCorrelationCutoff);%indices of measured Fourier points to enforce with complex values

    
    %To use a constant resolution cutoff, we just need to adjust the
    %weights to be 1 where we want to enforce and 0 everywhere else
    tmpWeights1 = zeros(size(constraintConfidenceWeights1));
    tmpWeights1(constraintInd_complex1) = constraintConfidenceWeights1(constraintInd_complex1);
    tmpWeights2 = zeros(size(constraintConfidenceWeights2));
    tmpWeights2(constraintInd_complex2) = constraintConfidenceWeights2(constraintInd_complex2);
        if isempty(initialObject) %if no initial object was provided, just use the IFFT of the constraints as initial object
%             [recHalf1, errK1] = GENFIRE_iterate(numIterations,smooth3D(initObject1,2*averageStrongCorrelationCutoff),support,measuredK1,constraintConfidenceWeights1);
%             [recHalf2, errK2] = GENFIRE_iterate(numIterations,smooth3D(initObject2,2*averageStrongCorrelationCutoff),support,measuredK2,constraintConfidenceWeights2);
            [recHalf1, errK1] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,measuredK1,tmpWeights1,constraintEnforcementDelayWeights);
            [recHalf2, errK2] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,measuredK2,tmpWeights2,constraintEnforcementDelayWeights);
       
        else
            [recHalf1, errK1] = GENFIRE_iterate(numIterations,initialObject,support,measuredK1,constraintConfidenceWeights1,constraintEnforcementDelayWeights);
            [recHalf2, errK2] = GENFIRE_iterate(numIterations,initialObject,support,measuredK2,constraintConfidenceWeights2,constraintEnforcementDelayWeights); 
        end
    else
        foundGoodWeights = 0;
    end
        softMask1 = makeLooseMask(recHalf1);%making a unique loose mask prevents artificial biasing between two halves 
        softMask2 = makeLooseMask(recHalf2);
        %perform Fourier Shell correlation
        [corrCoeffsIteration spatialFrequency meanIntensity] = FourierShellCorrelate(recHalf1.*softMask1,recHalf2.*softMask2,numBins,pixelSize);%title('iterative FSC')
        figure, plot(spatialFrequency,corrCoeffsIteration,'b'),%title('FSC between 2 half datasets, iterative method')
        softMask1 = makeLooseMask(rec1.*support);
        softMask2 = makeLooseMask(rec2.*support);
        [corrCoeffsInversion spatialFrequency meanIntensity] = FourierShellCorrelate(rec1.*softMask1,rec2.*softMask2,numBins,pixelSize);%title('inversion FRC')
        hold on, plot(spatialFrequency,corrCoeffsInversion,'r'),
        title('FSC between 2 half datasets')
        ylabel('Correlation')
        xlabel('Spatial Frequency (A^-^1)')
        legend('iterative method','inversion method')

        rec = (rec1+rec2)./2; %combine two halves of non-iterative reconstruction (for the iterative one we must reconstruct the Fourier grid)
        fprintf('GENFIRE: Independent half-reconstructions complete. Combining all data for final reconstruction.\n\n')
        %%deallocate memory
        clear recHalf1
        clear recHalf2
        clear rec1
        clear rec2
        clear measuredK1
        clear measuredK2
        clear softMask1
        clear softMask2

        %%finished with dividing in half and doing independent
        %%reconstructions for FSC comparison%%
end



initWeights = ones(size(projections,3),numBins);%initially try full resolution of every projection
FRCcoeffs = initWeights;
%%INTERPOLATE TO GRID
%%fillInFourierGrid_tomo_C contains the c++ code "weightVals.cpp" that must be
%%compiled on your local machine to run. If you cannot get it to work with
%%the command "mex weightVals.cpp" then use the version fillInFourierGrid
%%which does the same thing but slower.


fprintf('GENFIRE: Assembling Fourier grid...\n\n');
switch griddingMethod
    case 1
        [recIFFT measuredK,constraintConfidenceWeights,weightedDistances, sigmaPhases] = fillInFourierGrid_tomo_C(projections,angles,oversamplingRatioX,oversamplingRatioY,interpolationCutoffDistance,initWeights,doCTFcorrection);%interpolate projections to Fourier grid
    case 2
        [recIFFT, measuredK,constraintConfidenceWeights,weightedDistances] = fillInFourierGrid(projections,angles,oversamplingRatioX,oversamplingRatioY,interpolationCutoffDistance,initWeights,doCTFcorrection);%interpolate projections to Fourier grid
    case 3
        ori_projections = single(importdata(filename_Projections)); ori_projections = ori_projections(:,:,1:10);
        Typeind = 2;
        tic
        measuredK = My_fill_cubicgrid_ver3_1(dim, dim, ori_projections, angles, interpolationCutoffDistance, Typeind, newDim, newDim,1);
        interpTime = toc
        [null1 null2,constraintConfidenceWeights,weightedDistances] = fillInFourierGrid_tomo_C(projections,angles,interpolationCutoffDistance,initWeights);%interpolate projections to Fourier grid
        recIFFT = real(my_ifft(measuredK));
end
if exist('sigmaPhases','var')
    measuredK(sigmaPhases>phaseErrorSigmaTolerance) = 0;
end
resolutionWeights = zeros(size(Kmags)); %% first reconstruction uses resolution extension/suppression, where lower resolution information is enforced
%%initially and the maximum enforced resolution increases. This is followed by resolution where suppression, which is the same process run backwards, so at the final iteration 
%%only the lowest resolution is being enforced again.
resolutionWeights(measuredK~=0) = 1-Kmags(measuredK~=0);%make lower resolution have higher confidence. Higher confidence means enforced earlier in the reconstruction and for longer overall than low confidence
constraintConfidenceWeights = resolutionWeights;
constraintEnforcementDelayWeights = [.95:-.05:-.15  -.15:.05:.95];%This means the reconstruction will start by only enforcing data above confidence 0.95, which corresponds to the lowest 5% of resolution based on how we just constructed the resolutionWeights. 

recIFFTsmall = recIFFT(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);%take back original sized array of the initial IFFT to compare before/after iteration
% recIFFTsmall = recIFFT(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%take back original sized array of the initial IFFT to compare before/after iteration
tic%start clock
reconstructionIdentifierNumber = 0;%The reconstruction is performed twice, this tracks which reconstruction we are currently on

doneFlag = 0;%Once convergence of the resolution cutoffs is found, the final reconstruction is already completed, so this flag is to skip re-doing the final reconstruction
while ~foundGoodWeights
    
    if reconstructionIdentifierNumber == 0
        %%remove values for R_free calculation on first reconstruction
        measuredPointInd_complex = find(measuredK~=0&Kmags>=(averageStrongCorrelationCutoff+resRange)&Kmags<averageStrongCorrelationCutoff); %candidate values for Rfree_complex
        P = randperm(numel(measuredPointInd_complex)); %shuffle values
        measuredPointInd_complex = measuredPointInd_complex(P); %apply shuffle
        cutoffInd_complex = floor(numel(measuredPointInd_complex).*percentValuesForRfree); %take indices for 5% of measured data
        R_freeInd_complex = measuredPointInd_complex(1:cutoffInd_complex);%take complex value for 5% of measured data

        %now create a temporary set of constraints that have this 5% of
        %datapoints removed
        tmpMeasuredK = measuredK;
        tmpMeasuredK(R_freeInd_complex) = 0;
        R_freeVals_complex = measuredK(R_freeInd_complex);

        fprintf('GENFIRE: Reconstructing... \n\n');
        if isempty(initialObject)
            [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,tmpMeasuredK,resolutionWeights,constraintEnforcementDelayWeights,R_freeInd_complex,R_freeVals_complex);   
        else
            [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,initialObject,support,tmpMeasuredK,resolutionWeights,constraintEnforcementDelayWeights,R_freeInd_complex,R_freeVals_complex);
        end
    else
        fprintf('GENFIRE: Reconstructing... \n\n');
        if isempty(initialObject)
            [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,measuredK.*constraintConfidenceWeights,resolutionWeights,constraintEnforcementDelayWeights,R_freeInd_complex,R_freeVals_complex);   
        else
            [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,initialObject,support,measuredK.*constraintConfidenceWeights,resolutionWeights,constraintEnforcementDelayWeights,R_freeInd_complex,R_freeVals_complex);
        end
    end

    GENFIRE_recSmall = GENFIRE_rec(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);%extract the original sized result from the oversampled reconstruction
%     paddingForBackprojection = round(max(size(GENFIRE_recSmall))/2);%padding to prepare for backprojection calculation. If there is not enough 0 padding the projections are extremely inaccurate
%     GENFIRE_recSmall = padarray(GENFIRE_recSmall,[paddingForBackprojection paddingForBackprojection paddingForBackprojection]);
%     unpaddedProjections = single(importdata(filename_Projections));
%     centralPixel = size(projections,2)/2+1;
%     halfWindowSize = particleWindowSize/2;
    %backproject to same angles used for reconstruction and compare FRC to
    %input projections
    if reconstructionIdentifierNumber==0%only do backprojection comparison once
        GENFIRE_rec_resolutionExtensionSuppression = GENFIRE_rec(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);
        [percentFourierGridFilled, spatialFrequency] = percentageFourierGridFilledIn(measuredK,numBins,pixelSize);
        save(filename_results,'GENFIRE_rec_resolutionExtensionSuppression','percentFourierGridFilled','spatialFrequency')
%         save(filename_results,'GENFIRE_rec_resolutionExtensionSuppression');
        doneFlag = 1;
        break
        fprintf('GENFIRE: Computing FRC between backprojections of reconstruction and input projections...\n\n')
         projectionsTmp = projections(centralPixel-halfWindowSize:centralPixel+halfWindowSize-1,centralPixel-halfWindowSize:centralPixel+halfWindowSize-1,:);
         projectionsTmp = padarray(projectionsTmp,[paddingForBackprojection paddingForBackprojection 0]);%zero pad projections
        [meanFRCcoeff, spatialFrequency,newFRCcoeffs] = backprojectAndCompareFRC(projectionsTmp,angles,GENFIRE_recSmall,numBins,pixelSize,averageStrongCorrelationCutoff,targetCorrelation);

    end
    convergenceStatus = testConvergence(Rfree_complex);%check that the error converged 
    finalRfree_complex = Rfree_complex(end); %record final value of Rfree_complex
    
    if reconstructionIdentifierNumber==1;% reconstructionIdentifierNumber = 0 is the first step where resolution extension is used, reconstructionIdentifierNumber = 1 is the correlation weight step (the final step)
        foundGoodWeights = 1;
        doneFlag = 1;
        break
    else
        reconstructionIdentifierNumber = reconstructionIdentifierNumber+1;%update to reflect we are moving to the second reconstruction
        FRCcoeffs = newFRCcoeffs;
        fprintf('GENFIRE: Assembling newest Fourier grid...\n\n');

        switch griddingMethod%%reassemble Fourier grid and weights with the coefficients FRCcoeffs from the backprojection operation
            case 1
                [recIFFT, measuredK,constraintConfidenceWeights] = fillInFourierGrid_tomo_C(projections,angles,oversamplingRatioX,oversamplingRatioY,interpolationCutoffDistance,FRCcoeffs,doCTFcorrection);
            case 2
                [recIFFT measuredK,constraintConfidenceWeights] = fillInFourierGrid(projections,angles,particleWindowSize,oversamplingRatio,interpolationCutoffDistance,FRCcoeffs,doCTFcorrection);
            case 3
                [null1, null2,constraintConfidenceWeights] = fillInFourierGrid_tomo_C(projections,angles,particleWindowSize,oversamplingRatio,interpolationCutoffDistance,FRCcoeffs,doCTFcorrection);%reassemble just weights
        end
    end

end
if ~doneFlag%This loop should only run if the auto resolution finder was overridden and the user wants to use a single sharp cutoff for complex enforcement across all projections
    if averageStrongCorrelationCutoff == 1;
        averageStrongCorrelationCutoff = 1.5; %% if all of the information is enforceable, I choose to also enforce super resolution as simulations 
        %% indicated doing so produced better reco nstructions
    end
    if  ~runAutoResolutionCutoffFinder
        constraintInd_complex = find(measuredK~=0&Kmags<=averageStrongCorrelationCutoff);%indices of measured Fourier points to enforce with complex values
    else
        constraintInd_complex = find(measuredK~=0);%indices of measured Fourier points to enforce with complex values
    end
    
    %To use a constant resolution cutoff, we just need to adjust the
    %weights to be 1 where we want to enforce and 0 everywhere else
    tmpWeights = zeros(size(constraintConfidenceWeights));
    tmpWeights(constraintInd_complex) = constraintConfidenceWeights(constraintInd_complex);
    constraintConfidenceWeights = tmpWeights;
    clear tmpWeights
    recIFFT = real(my_ifft(constraintConfidenceWeights.*measuredK));%this is the IFFT of the constraints we are about to enforce, the "weights" are just 1 or 0 here
    recIFFTsmall = recIFFT(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);%take back original sized array of the initial IFFT to compare before/after iteration
    %%object if none was inputted) to the resolution that we can enforce.
    fprintf('GENFIRE: Reconstructing... \n\n');

    if isempty(initialObject) 
        [GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,zeros(size(measuredK)),support,measuredK,constraintConfidenceWeights,0);%here the weight parameter is just 0 which means all Fourier points will be enforced for all iterations
    else
        [GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,initialObject,support,measuredK,constraintConfidenceWeights,0);

    end
end


reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('GENFIRE: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);


%display results
if runAutoResolutionCutoffFinder
    GENFIRE_rec_correlationExtensionSuppression = GENFIRE_rec(ncK1-n2X:ncK1+n2X-1,ncK2-n2Y:ncK2+n2Y-1,ncK3-n2X:ncK3+n2X-1);%extract the original sized result from the oversampled reconstruction
    nc = round((size(GENFIRE_rec_correlationExtensionSuppression,2)+1)/2);
    figure,
    subplot(2,3,4), imagesc(squeeze(sum(GENFIRE_rec_correlationExtensionSuppression,1))),title('GENFIRE projection 1')
    subplot(2,3,5), imagesc(squeeze(sum(GENFIRE_rec_correlationExtensionSuppression,2))),title('GENFIRE projection 2')
    subplot(2,3,6), imagesc(squeeze(sum(GENFIRE_rec_correlationExtensionSuppression,3))),title('GENFIRE projection 3')
    subplot(2,3,1), imagesc(squeeze(sum(recIFFTsmall,1))),title('before iteration projection 1')
    subplot(2,3,2), imagesc(squeeze(sum(recIFFTsmall,2))),title('before iteration projection 2')
    subplot(2,3,3), imagesc(squeeze(sum(recIFFTsmall,3))),title('before iteration projection 3')


    figure,
    subplot(2,3,4), imagesc(squeeze(GENFIRE_rec_correlationExtensionSuppression(size(GENFIRE_rec_correlationExtensionSuppression,1)/2+1,:,:))),title('GENFIRE slice 1')
    subplot(2,3,5), imagesc(squeeze(GENFIRE_rec_correlationExtensionSuppression(:,size(GENFIRE_rec_correlationExtensionSuppression,2)/2+1,:))),title('GENFIRE slice 2')
    subplot(2,3,6), imagesc(squeeze(GENFIRE_rec_correlationExtensionSuppression(:,:,size(GENFIRE_rec_correlationExtensionSuppression,3)/2+1))),title('GENFIRE slice 3')
    subplot(2,3,1), imagesc(squeeze(recIFFTsmall(size(GENFIRE_rec_correlationExtensionSuppression,1)/2+1,:,:))),title('before iteration slice 1')
    subplot(2,3,2), imagesc(squeeze(recIFFTsmall(:,size(GENFIRE_rec_correlationExtensionSuppression,2)/2+1,:))),title('before iteration slice 2')
    subplot(2,3,3), imagesc(squeeze(recIFFTsmall(:,:,size(GENFIRE_rec_correlationExtensionSuppression,3)/2+1))),title('before iteration slice 3')
else
    GENFIRE_rec_correlationExtensionSuppression = GENFIRE_rec(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction
    nc = round((size(GENFIRE_rec_correlationExtensionSuppression,2)+1)/2);
    %display results
    figure,
    subplot(2,3,4), imagesc(squeeze(sum(GENFIRE_rec_correlationExtensionSuppression,1))),title('GENFIRE projection 1')
    subplot(2,3,5), imagesc(squeeze(sum(GENFIRE_rec_correlationExtensionSuppression,2))),title('GENFIRE projection 2')
    subplot(2,3,6), imagesc(squeeze(sum(GENFIRE_rec_correlationExtensionSuppression,3))),title('GENFIRE projection 3')
    subplot(2,3,1), imagesc(squeeze(sum(recIFFTsmall,1))),title('before iteration projection 1')
    subplot(2,3,2), imagesc(squeeze(sum(recIFFTsmall,2))),title('before iteration projection 2')
    subplot(2,3,3), imagesc(squeeze(sum(recIFFTsmall,3))),title('before iteration projection 3')


    figure,
    subplot(2,3,4), imagesc(squeeze(GENFIRE_rec_correlationExtensionSuppression(size(GENFIRE_rec_correlationExtensionSuppression,1)/2+1,:,:))),title('GENFIRE slice 1')
    subplot(2,3,5), imagesc(squeeze(GENFIRE_rec_correlationExtensionSuppression(:,size(GENFIRE_rec_correlationExtensionSuppression,2)/2+1,:))),title('GENFIRE slice 2')
    subplot(2,3,6), imagesc(squeeze(GENFIRE_rec_correlationExtensionSuppression(:,:,size(GENFIRE_rec_correlationExtensionSuppression,3)/2+1))),title('GENFIRE slice 3')
    subplot(2,3,1), imagesc(squeeze(recIFFTsmall(size(GENFIRE_rec_correlationExtensionSuppression,1)/2+1,:,:))),title('before iteration slice 1')
    subplot(2,3,2), imagesc(squeeze(recIFFTsmall(:,size(GENFIRE_rec_correlationExtensionSuppression,2)/2+1,:))),title('before iteration slice 2')
    subplot(2,3,3), imagesc(squeeze(recIFFTsmall(:,:,size(GENFIRE_rec_correlationExtensionSuppression,3)/2+1))),title('before iteration slice 3')
    
    
end
%save results
if runAutoResolutionCutoffFinder
    if ComputeFourierShellCorrelation
        save(filename_results,'GENFIRE_rec_correlationExtensionSuppression','GENFIRE_rec_resolutionExtensionSuppression','errK','corrCoeffsInversion','corrCoeffsIteration','percentFourierGridFilled','spatialFrequency')
    else
        save(filename_results,'GENFIRE_rec_correlationExtensionSuppression','GENFIRE_rec_resolutionExtensionSuppression','errK','spatialFrequency','percentFourierGridFilled')
    end
else
    if ComputeFourierShellCorrelation
        save(filename_results,'GENFIRE_rec','errK','corrCoeffsInversion','corrCoeffsIteration','spatialFrequency')
    else
        save(filename_results,'GENFIRE_rec','errK')
    end
end














