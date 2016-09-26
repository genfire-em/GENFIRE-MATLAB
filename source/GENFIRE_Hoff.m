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
%% Filtered Back-Projection. In addition, GENFIRE is capable of dynamically identifying
%% the noise level of the input projections (out to what resolution does useful information exist
%% before noise dominates the signal) without any prior knowledge from the user. 

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   User Parameters   %%%

% filename_Projections = '.\models\vesicle\projections.mat';%%filename of projections, which should be size NxNxN_projections where N_projections is the number of projections
% filename_Projections = 'C:\Users\apryo_000\Documents\MATLAB\singleParticleIterativeReconstruction\simulation\relion\100Fresh\projections100.mat'  ;%%filename of projections, which should be size NxNxN_projections where N_projections is the number of projections
% filename_Projections = 'C:\Users\apryo_000\Documents\MATLAB\singleParticleIterativeReconstruction\simulation\relion\100Fresh\projections100.mat'  ;%%filename of projections, which should be size NxNxN_projections where N_projections is the number of projections
filename_Projections = './projections.mat';%%filename of projections, which should be size NxNxN_projections where N_projections is the number of projections


% filename_Angles = '.\models\betagal\angles.mat';%%angles can be either a 1xN_projections array containing a single tilt series, or
% filename_Angles = 'C:\Users\apryo_000\Documents\MATLAB\singleParticleIterativeReconstruction\simulation\relion\100Fresh\angles100.mat';%%angles can be either a 1xN_projections array containing a single tilt series, or
% filename_Angles = '.\models\betagal\angles.mat';%%angles can be either a 1xN_projections array containing a single tilt series, or
filename_Angles = './angles.mat';%%angles can be either a 1xN_projections array containing a single tilt series, or


%%a 3xN_projections array containing 3 Euler angles for each projections in the form [phi;theta;psi]

filename_Support = './support.mat'; %% NxNxN binary array specifying a region of 1's in which the reconstruction can exist 
% filename_Support = 'C:\Users\apryo_-1\Documents\MATLAB\singleParticleIterativeReconstruction\simulation\relion\100Fresh\support.mat'; %% NxNxN binary array specifying a region of 1's in which the reconstruction can exist 

% filename_InitialModel = 'C:\Users\apryo_000\Documents\MATLAB\singleParticleIterativeReconstruction\simulation\relion\100proj\RELION_rec.mat';
% 
filename_results = './GENFIRE_rec.mat';
% filename_results = 'C:\Users\apryo_000\Documents\MATLAB\singleParticleIterativeReconstruction\simulation\relion\100Fresh\tests\GENFIRE_rec_v1p3_weightedValueWeightedReplacement_linearInterp_50It.mat';

global numIterations 
numIterations = 100; 

global pixelSize
pixelSize = 2; 

oversamplingRatio =3; %%The code will zero-pad projections for you to the inputted oversampling ratio. If your projections are already oversampled
%%then set this to 1.

interpolationCutoffDistance =.7; %%radius of sphere (in pixels) within which to include measured datapoints 
%%when assembling the 3D Fourier grid

ComputeFourierShellCorrelation = 1; %%set to 1 to divide dataset in half, independently reconstruct, and compute Fourier Shell Correlation (FSC) between two halves.
numBins = 50; %number of bins for FRC averaging
%% If you do not need FRC, set ComputeFourierShellCorrelation to 0 for speed as the FRC calculation requires reconstructing everything twice

%%%%  choose whether or not to run auto resolution cutoff finder or to
%%%%  enter your own cutoffs
runAutoResolutionCutoffFinder = 1;%%set to 1 to automatically find a suitable resolution cutoff for enforceable data
useAverageProjectionCutoff = 0;%Recommended to set this parameter to 0, change with caution. Used by the auto resolution cutoff finder. Set to 1 to use a single resolution cutoff for all projections (fast), or set to 0 to use a separate cutoff for each projection (slower because it requires reassembling the grid multiple times, but more accurate)

%to override the auto-search and set your own cutoffs for resolution to enforce, set runAutoResolutionCutoffFinder = 0 then input
%your desired cutoff here
complexCutoff = 1; %value between 0 and 1 indicating percentage of maximum resolution to enforce as complex values
magnitudeCutoff = .0;%value between 0 and 1 indicating percentage of maximum resolution to enforce as magnitude only
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%     Parameters that are unlikely to need to changing        %%%
targetCorrelation = 0.75;
Rfree_complex_target = 0.55;%value between 0 and 1 indicating the maximum acceptable convergence value for Rfree_complex for complex constraint cutoff (I recommend values between 0.5 to 0.65)
% Rfree_magnitude_target = 0.40;%value between 0 and 1 indicating the maximum acceptable convergence value for Rfree_magnitude for magnitude constraint cutoff (I recommend values between 0.3 to 0.4)
%%To override and manually set your own cutoffs, set runAutoResolutionCutoffFinder = 0 and input your
%%choices for complexCutoff and magnitudeCutoff below
%%%%****IN CURRENT IMPLEMENTATION ONLY COMPLEX VALUES ARE ENFORCED
%%%%****FUTURE IMPLEMENTATIONS WILL CONTAIN PHASE RETRIEVAL FOR PARTIAL
%%%%****MAGNITUDES
percentValuesForRfree = 0.05;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Begin Reconstruction   %%%
clc
angles = single(importdata(filename_Angles));%casting to single saves memory
projections = single(importdata(filename_Projections));

% projections = projections(:,:,1:10);
% angles = angles(:,1:10);

if exist('filename_InitialModel','var')
    initialObject = single(importdata(filename_InitialModel));
else
    initialObject = [];
end


global support %make support variable globally accessable to avoid passing copies of large arrays around to different functions
support = single(importdata(filename_Support));
n2 = size(support,2)/2;%array radius


%zero pad projections to user-inputted oversampling ratio
numProj = size(projections,3);
dim = size(projections,2);%original array size
newDim = dim*oversamplingRatio;%size of oversampled array
padding = (newDim-dim)/2;%how many zeros to add
nc = newDim/2+1;%central pixel of larger array
n2 = dim/2;%radius of smaller array
projections = padarray(projections,[padding padding 0]);%zero pad projections
support = padarray(support,[padding padding padding]); %%zero pad support
if exist('filename_InitialModel','var')
    initialObject =  padarray(initialObject,[padding padding padding]); 
end
%If there is only one tilt angle (such as in tomography), fill in other
%Euler angles with zeros

if size(angles,1)>3
    error('The dimension of the angles is incorrect.')
end

if size(angles,1) ==1 
angles = [angles;zeros(1,length(angles));zeros(1,length(angles))];
end

%construct K-space indices
if mod(size(support,1),2)==0
ncK1 = size(support,1)/2+1;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1-1)./n2K1;
elseif size(support,1)==1
vec1 = 0;
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
else
ncK3 = (size(support,3)+1)/2;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3)./n2K3; 
end
[Kx Ky Kz] = meshgrid(vec2,vec1,vec3);%grid of Fourier indices
Kmags = sqrt(Kx.^2+Ky.^2+Kz.^2);%take magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if runAutoResolutionCutoffFinder %automatically find suitable resolution limit to enforce
    foundBestCutoff = 0;
    complexCutoff = 1;
else
    fprintf('GENFIRE: Automatic resolution finder overwritten, continuing with user-inputted complex enforcement cutoff at %.12g%%\n\n',100*complexCutoff)
    foundBestCutoff = 1;
end
        
resRange = -0.05;%thickness of resolution ring to use for removal of datapoints for Rfree test
count = 0;%counts number of times the resoluton cutoff finding loop has run
if ComputeFourierShellCorrelation
    fprintf('GENFIRE: Dividing dataset in half...\n\n')
    pjCutoff1 = round((numProj+1)/2);%divide dataset in half
    pjCutoff2 = numProj;
    fprintf('GENFIRE: Assembling Fourier grids...\n\n');
    
    %%INTERPOLATE TO GRID
    %%fillInFourierGrid_C contains the c++ code "weightVals.cpp" that must be
    %%compiled on your local machine to run. If you cannot get it to work with
    %%the command "mex weightVals.cpp" then use the version fillInFourierGrid
    %%which does the same thing but slower.
    projectionCutoffValues1 = ones(1,size(projections(:,:,1:pjCutoff1),3));%initially try full resolution of every projection
    projectionCutoffValues2 = ones(1,size(projections(:,:,(pjCutoff1+1):pjCutoff2),3));%initially try full resolution of every projection
    initWeights = ones(size(projections,3),numBins);
    [initObject1 measuredK1,constraintConfidenceWeights1] = fillInFourierGrid_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,initWeights);%interpolate projections to Fourier grid
    [initObject2 measuredK2,constraintConfidenceWeights2] = fillInFourierGrid_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,initWeights);
  % [initObject1 measuredK1,constraintConfidenceWeights1] = fillInFourierGrid(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,initWeights);%interpolate projections to Fourier grid
  % [initObject2 measuredK2,constraintConfidenceWeights2] = fillInFourierGrid(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,initWeights);
    rec1 = initObject1;
    rec2 = initObject2;

    if runAutoResolutionCutoffFinder %automatically find suitable resolution limit to enforce
        foundBestCutoff = 0;
        complexCutoff = 1;
    else
        foundBestCutoff = 1;
    end
    
    while ~foundBestCutoff
        if useAverageProjectionCutoff || ~runAutoResolutionCutoffFinder%if the average is used then only enforce points within resolution of cutoff
            constraintInd_complex1 = find(measuredK1~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
            constraintInd_complex2 = find(measuredK2~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
        else
            constraintInd_complex1 = find(measuredK1~=0);%indices of measured Fourier points to enforce with complex values
            constraintInd_complex2 = find(measuredK2~=0);
        end
        constraintInd_magnitude1 = [];%This is a null placeholder in current version.
        constraintInd_magnitude2 = [];
        
        measuredPointInd_complex = find(measuredK1~=0&Kmags>=(complexCutoff+resRange)&Kmags<complexCutoff); %candidate values for Rfree_complex
        P = randperm(numel(measuredPointInd_complex)); %shuffle values
        measuredPointInd_complex = measuredPointInd_complex(P);
        cutoffInd_complex = floor(numel(measuredPointInd_complex).*percentValuesForRfree); %take 5% of values
        R_freeInd_complex = measuredPointInd_complex(1:cutoffInd_complex);
        constraintInd_complex1 = setdiff(constraintInd_complex1,R_freeInd_complex); %withhold the Rfree values from constraints

        %perform reconstruction
        fprintf('GENFIRE: Reconstructing... \n\n');
        if isempty(initialObject) %if no initial object was provided, just use the IFFT of the constraints as initial object
            [recHalf1, errK1,Rfree_complex] = GENFIRE_iterate(numIterations,smooth3D(initObject1,2*complexCutoff),support,measuredK1,constraintConfidenceWeights1,constraintInd_complex1,constraintInd_complex1,[]);
            [recHalf2, errK2] = GENFIRE_iterate(numIterations,smooth3D(initObject2,2*complexCutoff),support,measuredK2,constraintConfidenceWeights2,constraintInd_complex2,[],[]);
        else
            [recHalf1, errK1,Rfree_complex] = GENFIRE_iterate(numIterations,initialObject,support,measuredK1,constraintConfidenceWeights1,constraintInd_complex1,constraintInd_complex1,[]);
            [recHalf2, errK2] = GENFIRE_iterate(numIterations,initialObject,support,measuredK2,constraintConfidenceWeights2,constraintInd_complex2,[],[]); 
        end
        
        GENFIRE_recSmall1 = recHalf1(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction
        GENFIRE_recSmall2 = recHalf2(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction

        padding = round(sqrt(2)*max(size(GENFIRE_recSmall1))/4);%pad array for calculating backprojections
        GENFIRE_recSmall1 = padarray(GENFIRE_recSmall1,[padding padding padding]);
        GENFIRE_recSmall2 = padarray(GENFIRE_recSmall2,[padding padding padding]);

        unpaddedProjections = single(importdata(filename_Projections));
        unpaddedProjections = unpaddedProjections(:,:,:);
        
        %backproject to same angles used for reconstruction and compare to
        %input projections
        fprintf('GENFIRE: Computing FRC between backprojections of reconstruction and input projections...\n\n')
        [meanFRCcoeff, invResInd,newCutoff,newProjectionCutoffValues1,FRCcoeffs1] = backprojectAndCompareFRC(padarray(unpaddedProjections(:,:,1:pjCutoff1),[padding padding 0]),angles(:,1:pjCutoff1),GENFIRE_recSmall1,numBins,pixelSize,complexCutoff,targetCorrelation);
        [meanFRCcoeff, invResInd,newCutoff,newProjectionCutoffValues2,FRCcoeffs2] = backprojectAndCompareFRC(padarray(unpaddedProjections(:,:,(pjCutoff1+1):pjCutoff2),[padding padding 0]),angles(:,(pjCutoff1+1):pjCutoff2),GENFIRE_recSmall2,numBins,pixelSize,complexCutoff,targetCorrelation);

        convergenceStatus = testConvergence(Rfree_complex);%check that the Rfree converged 
        finalRfree_complex = Rfree_complex(end); %record final value of Rfree_complex
        if useAverageProjectionCutoff
            %if the change between the current cutoff and the new one is
            %very small, and if Rfree_complex converged to a value below
            %the target, then convergence has been found and we will use
            %that resolution cutoff
            if abs(newCutoff-complexCutoff)<0.02&&convergenceStatus&&finalRfree_complex<Rfree_complex_target;
                foundBestCutoff = 1;%turn on flag to terminate loop
                complexCutoff = max([newCutoff complexCutoff]);%take maximum of the two very similar values
                fprintf('GENFIRE: Convergence reached! Successfully found complex enforceable-resolution cutoff at %4.4g Angstrom cutoff (%.12g%% of maximum resolution)\n\n',1/(round(100*complexCutoff)./2./pixelSize./100),round(10000*complexCutoff)./100)
                
            else
                complexCutoff = newCutoff;%update cutoff
                count = count+1;%update number of times loop has run

                if count>10%if the loop does not converge after 10 iterations, something is wrong
                    error('GENFIRE: Cutoff search has stagnated; exiting.\n\n')
                end
                fprintf('GENFIRE: Assembling newest Fourier grid...\n\n');

            end
        else
            %if individual cutoffs are being used for each projection, then
            %convergence is found when the average change between cutoffs
            %is very small. The condition for Rfree_complex is the same as
            %in the first case
            if mean(abs(projectionCutoffValues1-newProjectionCutoffValues1))<0.02&&convergenceStatus&&finalRfree_complex<Rfree_complex_target;%check for convergence
                foundBestCutoff = 1;
                complexCutoff = max([newCutoff complexCutoff]);%take highest cutoff
                fprintf('GENFIRE: Convergence reached! Successfully found complex enforceable-resolution cutoffs with an average cutoff of %4.4g Angstroms (%.12g%% of maximum resolution)\n\n',1/(round(100*(mean([newProjectionCutoffValues1 newProjectionCutoffValues2])))./2./pixelSize./100),round(10000*mean([newProjectionCutoffValues1 newProjectionCutoffValues2]))./100)
            else
                count = count+1;%track how many times the cutoff search has run

                if count>10
                    error('GENFIRE: Cutoff search has stagnated; exiting.\n\n')
                end
                fprintf('GENFIRE: Assembling newest Fourier grid...\n\n');
                projectionCutoffValues1 = newProjectionCutoffValues1;%update resolution cutoffs
                projectionCutoffValues2 = newProjectionCutoffValues2;
                [rec1, measuredK1,constraintConfidenceWeights1] = fillInFourierGrid_C(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,FRCcoeffs1);%reassemble grid with final cutoff values
                [rec2, measuredK2,constraintConfidenceWeights2] = fillInFourierGrid_C(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,FRCcoeffs2);
%                 [rec1 measuredK1,constraintConfidenceWeights1] = fillInFourierGrid(projections(:,:,1:pjCutoff1),angles(:,1:pjCutoff1),interpolationCutoffDistance,FRCcoeffs1);%reassemble grid with final cutoff values
%                 [rec2 measuredK2,constraintConfidenceWeights2] = fillInFourierGrid(projections(:,:,(pjCutoff1+1):pjCutoff2),angles(:,(pjCutoff1+1):pjCutoff2),interpolationCutoffDistance,FRCcoeffs2);
                complexCutoff = mean([projectionCutoffValues1 projectionCutoffValues2]);
            end
        end


    end
    if ~runAutoResolutionCutoffFinder%if auto resolution searching is not turned on, then do the half reconstructions now
        if useAverageProjectionCutoff || ~runAutoResolutionCutoffFinder
            constraintInd_complex1 = find(measuredK1~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
            constraintInd_complex2 = find(measuredK2~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
        else
            constraintInd_complex1 = find(measuredK1~=0);%if the average cutoff is not being used, then the individual resolution cutoffs have already been included in the gridding process
            constraintInd_complex2 = find(measuredK2~=0);
        end
%         constraintInd_complex1 = find(measuredK1~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
        constraintInd_magnitude1 = [];%This is a null placeholder in current version.
%         constraintInd_complex2 = find(measuredK2~=0&Kmags<=complexCutoff);
        constraintInd_magnitude2 = [];
        
        %reconstruct
        if isempty(initialObject) %if no initial object was provided, just use the IFFT of the constraints as initial object
%             [recHalf1, errK1] = GENFIRE_iterate(numIterations,rec1,support,measuredK1,constraintInd_complex1,[],[]);
%             [recHalf2, errK2] = GENFIRE_iterate(numIterations,rec2,support,measuredK2,constraintInd_complex2,[],[]);
            [recHalf1, errK1] = GENFIRE_iterate(numIterations,smooth3D(initObject1,2*complexCutoff),support,measuredK1,constraintConfidenceWeights1,constraintInd_complex1,[],[]);
            [recHalf2, errK2] = GENFIRE_iterate(numIterations,smooth3D(initObject2,2*complexCutoff),support,measuredK2,constraintConfidenceWeights2,constraintInd_complex2,[],[]);
        else
            [recHalf1, errK1] = GENFIRE_iterate(numIterations,initialObject,support,measuredK1,constraintConfidenceWeights1,constraintInd_complex1,[],[]);
            [recHalf2, errK2] = GENFIRE_iterate(numIterations,initialObject,support,measuredK2,constraintConfidenceWeights2,constraintInd_complex2,[],[]); 
        end
    else
        foundBestCutoff = 0;
    end
        softMask1 = makeLooseMask(recHalf1);%making a unique loose mask prevents artificial biasing between two halves 
        softMask2 = makeLooseMask(recHalf2);
        %perform Fourier Shell correlation
        [corrCoeffsIteration invResInd meanIntensity] = FourierShellCorrelate(recHalf1.*softMask1,recHalf2.*softMask2,numBins,pixelSize);%title('iterative FSC')
        figure, plot(invResInd,corrCoeffsIteration,'b'),%title('FSC between 2 half datasets, iterative method')
        softMask1 = makeLooseMask(rec1.*support);
        softMask2 = makeLooseMask(rec2.*support);
        [corrCoeffsInversion invResInd meanIntensity] = FourierShellCorrelate(rec1.*softMask1,rec2.*softMask2,numBins,pixelSize);%title('inversion FRC')
        hold on, plot(invResInd,corrCoeffsInversion,'r'),
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


end

fprintf('GENFIRE: Assembling Fourier grid...\n\n');

%%INTERPOLATE TO GRID
%%fillInFourierGrid_C contains the c++ code "weightVals.cpp" that must be
%%compiled on your local machine to run. If you cannot get it to work with
%%the command "mex weightVals.cpp" then use the version fillInFourierGrid
%%which does the same thing but slower.
% projectionCutoffValues = complexCutoff.*ones(1,size(projections,3));%initially try full resolution of every projection
projectionCutoffValues = ones(1,size(projections,3));%initially try full resolution of every projection
initWeights = ones(size(projections,3),numBins);
FRCcoeffs = initWeights;
[recIFFT measuredK,constraintConfidenceWeights] = fillInFourierGrid_C(projections,angles,interpolationCutoffDistance,initWeights);%interpolate projections to Fourier grid
% [recIFFT, measuredK,constraintConfidenceWeights] = fillInFourierGrid(projections,angles,interpolationCutoffDistance,initWeights);%interpolate projections to Fourier grid
initObject = recIFFT;


recIFFTsmall = recIFFT(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%take back original sized array of the initial IFFT to compare before/after iteration
tic
count = 0;

doneFlag = 0;%Once convergence of the resolution cutoffs is found, the final reconstruction is already completed, so this flag is to skip re-doing the final reconstruction
if ~foundBestCutoff
    while ~foundBestCutoff
        if useAverageProjectionCutoff || ~runAutoResolutionCutoffFinder%if the average is used then only enforce points within resolution of cutoff
            constraintInd_complex = find(measuredK~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
            fprintf('GENFIRE: Trying enforceable-resolution cutoffs with an average cutoff at %4.4g Angstrom cutoff (%.12g%% of maximum resolution)...\n\n',1/(round(100*complexCutoff)./2./pixelSize./100),round(10000*complexCutoff)./100)

        else
%             constraintInd_complex = find(measuredK~=0);%if the average cutoff is not being used, then the individual resolution cutoffs have already been included in the gridding process
            constraintInd_complex = find(measuredK~=0&Kmags<=1);%if the average cutoff is not being used, then the individual resolution cutoffs have already been included in the gridding process

            fprintf('GENFIRE: Trying enforceable-resolution cutoffs with an average cutoff at %4.4g Angstrom cutoff (%.12g%% of maximum resolution)...\n\n',1/(round(100*mean(projectionCutoffValues))./2./pixelSize./100),round(10000*mean(projectionCutoffValues))./100)

        end
        constraintInd_magnitude = [];%This is a null placeholder in current version.
        
        measuredPointInd_complex = find(measuredK~=0&Kmags>=(complexCutoff+resRange)&Kmags<complexCutoff); %candidate values for Rfree_complex
        P = randperm(numel(measuredPointInd_complex)); %shuffle values
        measuredPointInd_complex = measuredPointInd_complex(P);
        cutoffInd_complex = floor(numel(measuredPointInd_complex).*percentValuesForRfree); %take 5% of values
        R_freeInd_complex = measuredPointInd_complex(1:cutoffInd_complex);
        constraintInd_complex = setdiff(constraintInd_complex,R_freeInd_complex); %withhold the Rfree values from constraints
        
        tmpMeasuredK = measuredK;
        tmpMeasuredK(constraintInd_complex) = measuredK(constraintInd_complex);
        R_freeVals_complex = measuredK(R_freeInd_complex);
        
        fprintf('GENFIRE: Reconstructing... \n\n');

%         [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,smooth3D(recIFFT,2*complexCutoff),support,measuredK,constraintInd_complex,R_freeInd_complex,[]);   
        [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,tmpMeasuredK,constraintConfidenceWeights,constraintInd_complex,R_freeInd_complex,[],R_freeVals_complex);   
%         [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(measuredK)),support,measuredK,constraintConfidenceWeights,constraintInd_complex,R_freeInd_complex,[]);   

        GENFIRE_recSmall = GENFIRE_rec(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction
        padding = round(sqrt(2)*max(size(GENFIRE_recSmall))/4);%padding to prepare for backprojection calculation
        GENFIRE_recSmall = padarray(GENFIRE_recSmall,[padding padding padding]);
        unpaddedProjections = single(importdata(filename_Projections));
%         unpaddedProjections = unpaddedProjections(:,:,1:10);
        %backproject to same angles used for reconstruction and compare FRC to
        %input projections
        fprintf('GENFIRE: Computing FRC between backprojections of reconstruction and input projections...\n\n')
        [meanFRCcoeff, invResInd,newCutoff,newProjectionCutoffValues,newFRCcoeffs] = backprojectAndCompareFRC(padarray(unpaddedProjections,[padding padding 0]),angles,GENFIRE_recSmall,numBins,pixelSize,complexCutoff,targetCorrelation);

        convergenceStatus = testConvergence(Rfree_complex);%check that the error converged 
        finalRfree_complex = Rfree_complex(end); %record final value of Rfree_complex
        if useAverageProjectionCutoff
            if abs(newCutoff-complexCutoff)<0.02&&convergenceStatus&&finalRfree_complex<Rfree_complex_target;
                foundBestCutoff = 1;
                complexCutoff = max([newCutoff complexCutoff]);
                fprintf('GENFIRE: Convergence reached! Successfully found complex cutoff at %4.4g Angstrom cutoff (%.12g%% of maximum resolution)\n\n',1/(round(100*complexCutoff)./2./pixelSize./100),round(10000*complexCutoff)./100)
                doneFlag = 1;

            else
                complexCutoff = newCutoff;
                count = count+1;

                if count>10
                    error('GENFIRE: Cutoff search has stagnated; exiting.\n\n')
                end
                fprintf('GENFIRE: Assembling newest Fourier grid...\n\n');

            end
        else
%             if mean(abs(projectionCutoffValues-newProjectionCutoffValues))<0.02&&convergenceStatus&&finalRfree_complex<Rfree_complex_target;%check for convergence
%             if mean(mean(abs(FRCcoeffs-newFRCcoeffs)))<0.02&&convergenceStatus&&finalRfree_complex<Rfree_complex_target 
            if mean(mean(abs(FRCcoeffs-newFRCcoeffs)))<0.02&&convergenceStatus&&finalRfree_complex<Rfree_complex_target || count==1;%check for convergence
    
                foundBestCutoff = 1;
                complexCutoff = max([newCutoff complexCutoff]);%take highest cutoff
                fprintf('GENFIRE: Convergence reached! Successfully found complex cutoff at average of %4.4g Angstrom cutoff (%.12g%% of maximum resolution)\n\n',1/(round(100*(mean(newProjectionCutoffValues)))./2./pixelSize./100),round(10000*mean(projectionCutoffValues))./100)
                doneFlag = 1;
            else
                count = count+1;%track how many times the cutoff search has run
                complexCutoff =  newCutoff;
                FRCcoeffs = newFRCcoeffs;
                if count>10
                    error('GENFIRE: Cutoff search has stagnated; exiting.\n\n')
                end
                fprintf('GENFIRE: Assembling newest Fourier grid...\n\n');
                projectionCutoffValues = newProjectionCutoffValues;
                [recIFFT, measuredK,constraintConfidenceWeights] = fillInFourierGrid_C(projections,angles,interpolationCutoffDistance,FRCcoeffs);%reassemble grid
%                 [recIFFT measuredK,constraintConfidenceWeights] = fillInFourierGrid(projections,angles,interpolationCutoffDistance,FRCcoeffs);%reassemble grid
            end
        end
    end
end
if ~doneFlag
    if complexCutoff == 1;
        complexCutoff = 1.5; %% if all of the information is enforceable, I choose to also enforce super resolution as simulations 
        %% indicated doing so produced better reco nstructions
    end
    if useAverageProjectionCutoff || ~runAutoResolutionCutoffFinder
        constraintInd_complex = find(measuredK~=0&Kmags<=complexCutoff);%indices of measured Fourier points to enforce with complex values
    else
        constraintInd_complex = find(measuredK~=0);%indices of measured Fourier points to enforce with complex values

    end
    if magnitudeCutoff ==0
        constraintInd_magnitude = [];
    else

    constraintInd_magnitude = find(measuredK~=0&Kmags<=magnitudeCutoff); 
    end
    % constraintInd_magnitude = [];%This is a null placeholder for future versions

%     recIFFT=smooth3D(recIFFT,complexCutoff*2);%To low pass filter the IFFT reconstruction (which will be our initial
    %%object if none was inputted) to the resolution that we can enforce.
    fprintf('GENFIRE: Reconstructing... \n\n');

    if isempty(initialObject) 
%         [GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,recIFFT,support,measuredK,constraintInd_complex,[],[]);%run iterations
%         [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,smooth3D(initObject,2*complexCutoff),support,measuredK,constraintConfidenceWeights,constraintInd_complex,[],[]);   
        [GENFIRE_rec, errK,Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(measuredK)),support,measuredK,constraintConfidenceWeights,constraintInd_complex,[],[]);   

    else
        %[GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,initialObject,support,measuredK,constraintInd_complex,constraintInd_magnitude,[],[]);
        [GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,initialObject,support,measuredK,constraintConfidenceWeights,constraintInd_complex,[],[]);

    end
end


reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('GENFIRE: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);

GENFIRE_rec = GENFIRE_rec(ncK1-n2:ncK1+n2-1,ncK2-n2:ncK2+n2-1,ncK3-n2:ncK3+n2-1);%extract the original sized result from the oversampled reconstruction

%display results
figure,
subplot(2,3,4), imagesc(squeeze(sum(GENFIRE_rec,1))),title('GENFIRE projection 1')
subplot(2,3,5), imagesc(squeeze(sum(GENFIRE_rec,2))),title('GENFIRE projection 2')
subplot(2,3,6), imagesc(squeeze(sum(GENFIRE_rec,3))),title('GENFIRE projection 3')
subplot(2,3,1), imagesc(squeeze(sum(recIFFTsmall,1))),title('before iteration projection 1')
subplot(2,3,2), imagesc(squeeze(sum(recIFFTsmall,2))),title('before iteration projection 2')
subplot(2,3,3), imagesc(squeeze(sum(recIFFTsmall,3))),title('before iteration projection 3')

nc = round((size(GENFIRE_rec,2)+1)/2);
figure,
subplot(2,3,4), imagesc(squeeze(GENFIRE_rec(nc,:,:))),title('GENFIRE slice 1')
subplot(2,3,5), imagesc(squeeze(GENFIRE_rec(:,nc,:))),title('GENFIRE slice 2')
subplot(2,3,6), imagesc(squeeze(GENFIRE_rec(:,:,nc))),title('GENFIRE slice 3')
subplot(2,3,1), imagesc(squeeze(recIFFTsmall(nc,:,:))),title('before iteration slice 1')
subplot(2,3,2), imagesc(squeeze(recIFFTsmall(:,nc,:))),title('before iteration slice 2')
subplot(2,3,3), imagesc(squeeze(recIFFTsmall(:,:,nc))),title('before iteration slice 3')

%save results
if ComputeFourierShellCorrelation
    save(filename_results,'GENFIRE_rec','errK','corrCoeffsInversion','corrCoeffsIteration','invResInd','recIFFTsmall','complexCutoff','FRCcoeffs')
else
    save(filename_results,'GENFIRE_rec','errK','FRCcoeffs','complexCutoff','recIFFTsmall')
end














