%%  GENFIRE_reconstruct %%

%% Primary control function for reconstructions

%%inputs:
%%  GENFIRE_parameters - struct containing reconstruction parameters


%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function GENFIRE_reconstruct(GENFIRE_parameters)

%unpack reconstruction parameters
filename_Projections = GENFIRE_parameters.filename_Projections;
filename_Angles = GENFIRE_parameters.filename_Angles;
filename_Support = GENFIRE_parameters.filename_Support;
filename_Results = GENFIRE_parameters.filename_Results;
filename_InitialModel = GENFIRE_parameters.filename_InitialModel;
numIterations = GENFIRE_parameters.numIterations;
pixelSize = GENFIRE_parameters.pixelSize;
oversamplingRatio = GENFIRE_parameters.oversamplingRatio;
interpolationCutoffDistance = GENFIRE_parameters.interpolationCutoffDistance;
constraintPositivity = GENFIRE_parameters.constraintPositivity;
constraintSupport = GENFIRE_parameters.constraintSupport;
particleWindowSize = GENFIRE_parameters.particleWindowSize;
numBins = GENFIRE_parameters.numBins;
percentValuesForRfree = GENFIRE_parameters.percentValuesForRfree;
numBinsRfree = GENFIRE_parameters.numBinsRfree;
doCTFcorrection = GENFIRE_parameters.doCTFcorrection;
griddingMethod = GENFIRE_parameters.griddingMethod;
allowMultipleGridMatches = GENFIRE_parameters.allowMultipleGridMatches;
phaseErrorSigmaTolerance = GENFIRE_parameters.phaseErrorSigmaTolerance;
constraintEnforcementDelayIndicators = GENFIRE_parameters.constraintEnforcementDelayIndicators;
%%%   Begin Reconstruction   %%%

if griddingMethod>2
   error('GENFIRE: Unrecognized gridding method.') 
end

angles = single(importdata(filename_Angles));
projections = single(importdata(filename_Projections));

%Initialize the initial object
if filename_InitialModel
    initialObject = single(importdata(GENFIRE_parameters.filename_InitialModel));
else
    initialObject = [];
end



global support %make support variable globally accessable to avoid passing copies of large arrays around to different functions
support = single(importdata(filename_Support));

%get some values related to the size, center, etc of this array size
vecX = 1:size(support,1); ncX = round((size(support,1)+1)/2); vecX = vecX - ncX;
vecY = 1:size(support,2); ncY = round((size(support,2)+1)/2); vecY = vecY - ncY;
vecZ = 1:size(support,3); ncZ = round((size(support,3)+1)/2); vecZ = vecZ - ncZ;


if isempty(particleWindowSize)
    particleWindowSize = size(support,2);
else
    if particleWindowSize ~= size(support,2)
        error('GENFIRE: ERROR! The size of your projections (set by particleWindowSize) does not match the dimensions of your support!')
    end
end

n2 = particleWindowSize/2;%radius of smaller array
newDim = particleWindowSize*oversamplingRatio;%size of oversampled array
padding = round((newDim-particleWindowSize)/2);%how many zeros to add

%zero pad projections to user-inputted oversampling ratio
numProj = size(projections,3);
support = padarray(support,[padding padding padding]); %%zero pad support
if initialObject
    initialObject =  padarray(initialObject,[padding padding padding]); 
end

%If there is only one tilt angle (such as in tomography), fill in other
%Euler angles with zeros
if size(angles,2)>3
    error('The dimension of the angles is incorrect.')
end
if size(angles,2) ==1 
angles = [zeros(length(angles),1),  angles, zeros(length(angles),1)];%tomography tilt is the theta angle
end


Q = make_Kspace_indices(support);
resRange = -0.05;%thickness of resolution ring to use for removal of datapoints for Rfree test


%%INTERPOLATE TO GRID
%%fillInFourierGrid contains the C++ code "weightVals.cpp" that must be
%%compiled on your local machine to run. If you cannot get it to work with
%%the command "mex weightVals.cpp" then use the version fillInFourierGrid
%%which does the same thing but is slower.
fprintf('GENFIRE: Assembling Fourier grid...\n\n');

switch griddingMethod
    case 1
        [recIFFT, measuredK ] = fillInFourierGrid(projections,angles,particleWindowSize,oversamplingRatio,interpolationCutoffDistance,doCTFcorrection, [], allowMultipleGridMatches);%interpolate projections to Fourier grid
    case 2
        [recIFFT, measuredK] = fillInFourierGrid_DFT(projections, angles, interpolationCutoffDistance, size(support,1), size(support,2), ones(size(projections,3),numBins), 1, 0, 0, []);
end

if exist('sigmaPhases','var') && ~isempty(phaseErrorSigmaTolerance)
    measuredK(sigmaPhases>phaseErrorSigmaTolerance) = 0;
end

resolutionIndicators = zeros(size(Q)); %% first reconstruction uses resolution extension/suppression, where lower resolution information is enforced
%%initially and the maximum enforced resolution increases. This is followed by resolution where suppression, which is the same process run backwards, so at the final iteration 
%%only the lowest resolution is being enforced again.
resolutionIndicators(measuredK~=0) = 1-Q(measuredK~=0);%make lower resolution have higher confidence. Higher confidence means enforced earlier in the reconstruction and for longer overall than low confidence


%remove datapoints for Rfree calculation
tmpMeasuredK = measuredK;
spatialFrequencyForRfree = linspace(0,1,numBinsRfree+1);%compute spatial frequency bins
for shellNum = 1:numBinsRfree %loop over each frequency shell
    
    measuredPointInd_complex = find(measuredK~=0&Q>=(spatialFrequencyForRfree(shellNum)+resRange)&Q<spatialFrequencyForRfree(shellNum+1)); %candidate values for Rfree_complex
    P = randperm(numel(measuredPointInd_complex)); %shuffle values
    measuredPointInd_complex = measuredPointInd_complex(P); %apply shuffle
    cutoffInd_complex = floor(numel(measuredPointInd_complex).*percentValuesForRfree); %take indices for 5% of measured data
    if cutoffInd_complex == 0 %make sure to include at least one point
        cutoffInd_complex = 1;
    end
    R_freeInd_complex{shellNum} = measuredPointInd_complex(1:cutoffInd_complex);%take complex value for 5% of measured data

    %now create a temporary set of constraints that have this 5% of
    %datapoints removed
    tmpMeasuredK(R_freeInd_complex{shellNum}) = 0;
    R_freeVals_complex{shellNum} = measuredK(R_freeInd_complex{shellNum});
end

%run the actual reconstruction
fprintf('GENFIRE: Reconstructing... \n\n');

if isempty(initialObject)
    initialObject = zeros(size(support));
end
[GENFIRE_rec, errK, Rfree_complex] = GENFIRE_iterate(numIterations,zeros(size(support),'single'),support,tmpMeasuredK,resolutionIndicators,constraintEnforcementDelayIndicators,R_freeInd_complex,R_freeVals_complex, constraintPositivity, constraintSupport);   

reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('GENFIRE: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);

ncX_big = round((size(GENFIRE_rec,1)+1)/2);
ncY_big = round((size(GENFIRE_rec,2)+1)/2);
ncZ_big = round((size(GENFIRE_rec,3)+1)/2);
GENFIRE_rec = GENFIRE_rec(vecX + ncX_big, vecY + ncY_big, vecZ + ncZ_big);%extract the original sized result from the oversampled reconstruction
recIFFT = recIFFT(vecX + ncX_big, vecY + ncY_big, vecZ + ncZ_big);%take back original sized array of the initial IFFT to compare before/after iteration

GENFIRE_parameters.reconstruction = GENFIRE_rec;
GENFIRE_parameters.errK = errK;
GENFIRE_parameters.Rfree_complex_bybin = Rfree_complex;
% GENFIRE_parameters.Rfree_complex_bybin = Rfree_complex_total;
% %display the results
figure,
subplot(2,3,4), imagesc(squeeze(sum(GENFIRE_rec,1))),title('GENFIRE projection 1')
subplot(2,3,5), imagesc(squeeze(sum(GENFIRE_rec,2))),title('GENFIRE projection 2')
subplot(2,3,6), imagesc(squeeze(sum(GENFIRE_rec,3))),title('GENFIRE projection 3')
subplot(2,3,1), imagesc(squeeze(sum(recIFFT,1))),title('before iteration projection 1')
subplot(2,3,2), imagesc(squeeze(sum(recIFFT,2))),title('before iteration projection 2')
subplot(2,3,3), imagesc(squeeze(sum(recIFFT,3))),title('before iteration projection 3')


figure,
subplot(2,3,4), imagesc(squeeze(GENFIRE_rec(ncX,:,:))),title('GENFIRE slice 1')
subplot(2,3,5), imagesc(squeeze(GENFIRE_rec(:,ncY,:))),title('GENFIRE slice 2')
subplot(2,3,6), imagesc(squeeze(GENFIRE_rec(:,:,ncZ))),title('GENFIRE slice 3')
subplot(2,3,1), imagesc(squeeze(recIFFT(ncX,:,:))),title('before iteration slice 1')
subplot(2,3,2), imagesc(squeeze(recIFFT(:,ncY,:))),title('before iteration slice 2')
subplot(2,3,3), imagesc(squeeze(recIFFT(:,:,ncZ))),title('before iteration slice 3')

%save results
save(filename_Results, 'GENFIRE_parameters')

end

