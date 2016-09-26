%%  GENFIRE_OSS %%

%% An implementation of OverSampling Smoothness (OSS). The primary difference
%% is that this implementation allows for enforcing some of the Fourier data
%% as complex values, and others as magnitude. This allows for interesting
%% reconstruction techniques such as phase retrieval with a known low resolution
%% envelope.This code is heavily  modeled after the original source code. If you use it, 
%% please cite the following paper: http://arxiv.org/abs/1211.4519

%%inputs:
%%  numIterations - number of iterations to run
%%  initialObject - inital guess of object
%%  support - region of 1's and 0's defining where the reconstruction can exist (voxels that are 0 in support will be set to 0 in reconstruction)
%%  measuredK - measured Fourier points
%%  constraintInd_complex - indices of measuredK that can be enforced as complex (magnitude and phase) values
%%  constraintInd_magnitude - indices of measuredK that can be enforced as magnitude values
%%  R_freeInd_mag - indices of 5% of points in the highest resolution shell of measuredK that are being withheld from 
%%          reconstruction and compared to after iteration.

%%outputs:
%%  rec - reconstruction after iteration
%%  errK - reciprocal space error
%%  Rfree_magnitude - value of magnitude Rfree vs iteration
%%  bestErrors - minimum reciprocal error for each filter

%% Authors: AJ Pryor, Jose Rodriguez
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [rec, errK, Rfree_magnitude, bestErrors] = GENFIRE_OSS(numIterations,initialObject,support,measuredK,constraintInd_complex,constraintInd_magnitude,R_freeInd_mag)

beta = 0.9; %beta parameter for reconstruction update
numFilters = 10; %number of progressively stronger smoothing filters to use

bestErr = 1e30;%initialize best error
Rfree_magnitude = -1*ones(1,numIterations,'single');%% initialize Rfree_magnitude curve
errK = zeros(1,numIterations,'single'); % initialize error
combinedConstraintInd = union(constraintInd_complex,constraintInd_magnitude);%combined indices of voxels that are enforced in any way
bufferObject = initialObject; %object from last iteration
store = 0;

%setup indices for filter and calculate the sigma values to use
Rsize = size(initialObject,1);
Csize = size(initialObject,2);
Lsize = size(initialObject,3);
Rcenter = round((Rsize+1)/2);
Ccenter = round((Csize+1)/2);
Lcenter = round((Lsize+1)/2);
a=1:1:Rsize;
b=1:1:Csize;
c=1:1:Lsize;
[bb,aa,cc]=meshgrid(b,a,c);
X=1:numIterations;
sigmas=(numFilters+1-ceil(X*numFilters/(numIterations)))*ceil(numIterations/(1*numFilters)); 
sigmas=((sigmas-ceil(numIterations/numFilters))*(2*Rsize)/max(sigmas))+(Rsize/20);

figure(98), plot(X,sigmas), axis tight; title('OSS Filter Size v. Iterations');
totalIterationNum = 1;
for filterNum = 1:numFilters
    %compute filter
    kfilter=exp( -( ( ((sqrt((aa-Rcenter).^2+(bb-Ccenter).^2 + (cc-Lcenter).^2)).^2 ) ./ (2* sigmas((filterNum-1)*(numIterations/numFilters)+1).^2) )));
    kfilter=kfilter/max(kfilter(:)); %normalize
    for iterationNum = 1:numIterations/numFilters
        if totalIterationNum == 1; %for magnitude enforcement, provide random phases for first iteration
            k = my_fft(initialObject);
            rng('shuffle','twister')
            randPhases = single(rand(size(k))*2*pi);
            k(constraintInd_magnitude) = abs(measuredK(constraintInd_magnitude)).*exp(1*i*randPhases(constraintInd_magnitude));
            bestErr = 1e30;%initialize best error
            initialObject = real(my_ifft(k));
        end
        sample = initialObject.*support; %isolate region inside support
        initialObject = bufferObject-beta*initialObject; %beta update
        sample(sample<0)=initialObject(sample<0); %positivity constraint

        initialObject=single(real(my_ifft(my_fft(initialObject).*kfilter)));%smoothing constraint

        if mod(totalIterationNum,ceil(numIterations/numFilters))==0&&totalIterationNum>10
            %start each new filter with the current best reconstruction
            initialObject=bestRec; 
        else
            %otherwise put back the sample that was just isolated and
            %updated
            initialObject(support==1)=sample(support==1);
        end
        %update the buffer object
        bufferObject = initialObject;
        k = my_fft(initialObject);%take FFT of initial object

        %monitor error using the masked object
        ksample = my_fft(sample);
        errK(totalIterationNum) = sum(abs(abs(ksample(combinedConstraintInd))-abs(measuredK(combinedConstraintInd))))./sum(abs(measuredK(combinedConstraintInd)));%magnitudes only

        if nargin > 6 %if values have been withheld from measuredK for monitoring R_free, check them accordingly
            if ~isempty(R_freeInd_mag)
                Rfree_magnitude(iterationNum) = sum(abs(abs(k(R_freeInd_mag))-abs(measuredK(R_freeInd_mag))))./sum(abs(measuredK(R_freeInd_mag)));
            end
        end
        
        if errK(totalIterationNum)<bestErr&&totalIterationNum >store+2 %if current reconstruction has better error, update best error and best reconstruction
            bestErr = errK(totalIterationNum);
            rec = sample;
            bestRec = initialObject;
            store=totalIterationNum;
        end


        %Fourier constraint -- replace magnitudes with measured data
        k(constraintInd_magnitude) = abs(measuredK(constraintInd_magnitude)).*exp(1i*angle(k(constraintInd_magnitude)));
        totalIterationNum = totalIterationNum+1; %update total number of iterations
        k(constraintInd_complex) = measuredK(constraintInd_complex);%replace complex measured values
        initialObject = real(my_ifft(k));%obtain next object with IFFT
        figure(3), subplot(1,2,1),plot(errK),subplot(1,2,2), imagesc(squeeze(sum(sample,3)));axis image,title('current object')
    end
    bestErrors(filterNum) = bestErr;
end
