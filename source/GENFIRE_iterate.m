%%  GENFIRE_iterate %%

%%inputs:
%%  numIterations - number of iterations to run
%%  initialObject - inital guess of object
%%  support - region of 1's and 0's defining where the reconstruction can exist (voxels that are 0 in support will be set to 0 in reconstruction)
%%  measuredK - measured Fourier points
%%  constraintIndicators - a flag value that is used to determine when a particular Fourier point is enforced, i.e. by resolution.
%%      The values to enforce at any particular iteration are determined by constraintEnforcementDelayIndicators.
%%  constraintEnforcementDelayIndicators - vector of values indicating the Fourier enforcement cutoff. Fourier grid points with a constraintIndicators value greater than
%%      or equal to the current constraintEnforcementDelayIndicators value will be enforced. The values in constraintEnforcementDelayIndicators are spread evenly over the number of iterations.
%%  R_freeInd_complex - indices of 5% of points in the highest resolution shell of measuredK that are being withheld from 
%%          reconstruction and compared to after iteration. 
%%  R_freeVals_complex - corresponding complex values at the indices in R_freeInd_complex

%%outputs:
%%  rec - reconstruction after iteration
%%  errK - reciprocal space error
%%  Rfree_complex - value of complex Rfree in each resolution shell vs iteration

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [rec, errK, Rfree_complex] = GENFIRE_iterate(numIterations,initialObject,support,measuredK,constraintIndicators,constraintEnforcementDelayIndicators,R_freeInd_complex,R_freeVals_complex, enforce_support)
if nargin < 9
    enforce_support = 1;
end
bestErr = 1e30;%initialize best error
Rfree_complex = -1*ones(1,numIterations,'single');%% initialize Rfree_complex curve , -1 is a flag that means undefined
errK = zeros(1,numIterations,'single');

%prefetch indices to use for error metric to avoid having to lookup each
%iteration
errInd = find(measuredK~=0);

%determine how to spread the provided weighting cutoffs over the iterations
iterationNumsToChangeCutoff = round(linspace(1,numIterations,numel(constraintEnforcementDelayIndicators)));
constraintInd_complex = find(constraintIndicators~=0&measuredK~=0);
currentCutoffNum = 1;
for iterationNum = 1:numIterations
    if iterationNum == iterationNumsToChangeCutoff(currentCutoffNum)
        currentCutoffNum = find(iterationNumsToChangeCutoff==iterationNum,1,'last');
        constraintInd_complex = find(constraintIndicators>(constraintEnforcementDelayIndicators(currentCutoffNum))&measuredK~=0);
        currentCutoffNum = currentCutoffNum+1;
        bestErr = 1e30;%reset best error
    end
if mod(iterationNum,1)==0
    iterationNum
end
initialObject(initialObject<0) = 0; %enforce positivity
if enforce_support
initialObject = initialObject.*support;%enforce support
end

k = my_fft(initialObject);%take FFT of initial object
%monitor error
errK(iterationNum) = sum(abs(k(errInd)-measuredK(errInd)))./sum(abs(measuredK(errInd)));

if nargin > 6 %if values have been withheld from measuredK for monitoring R_free, check them accordingly
    if ~isempty(R_freeInd_complex)
        %calculate Rfree in each resolution shell
        for shellNum = 1:numel(R_freeInd_complex)
            tmpInd =R_freeInd_complex{shellNum};
            tmpVals = R_freeVals_complex{shellNum};
            Rfree_complex(shellNum,iterationNum) = sum(abs(k(tmpInd)-tmpVals))./sum(abs(tmpVals));
        end
    end
end
if errK(iterationNum)<bestErr %if current reconstruction has better error, update best error and best reconstruction
%     fprintf('GENFIRE: new best object, iteration %d\n',iterationNum)
    bestErr = errK(iterationNum);
    rec = initialObject;
end

%enforce Fourier constraint
k(constraintInd_complex) = measuredK(constraintInd_complex);

initialObject = real(my_ifft(k));%obtain next object with IFFT

end
