%%  weightAverageCombiner %%

%Weight average two matrices that have been individually composed of a
%weighted average. If a value is the weighted average of several initial values,
%the result is the same as if all of the initial values were replaced with
%the final weighted average at the weighted distance. For example consider
%the simple linear average of two values, 1 and 3. Their average is 2. If I
%want to know the average of three numbers 1,3, and 5 I can obtain that by (1+3+5)/3 = 3
%or equivalently if I know the average of the first two values is 2 then I
%could compute the average of the 3 values as (2+2+5)/3 = 3. For GENFIRE
%this process is used when handling large numbers of images, as it is
%computationally cumbersome, if not impossible, to simultaneously store all images in RAM at once.
%Instead, the images are dynamically read in small batches, and a "rolling
%weighted average" is kept. The final result is the same as if all
%projections were weighted at once.

%%inputs:
%%  in1 - first matrix
%%  in2 - second matrix
%%  weightedDistances1 - weighted average distance of matrix 1 at time of composition
%%  weightedDistances2 - weighted average distance of matrix 2 at time of composition
%%  confidenceWeights1 - number of values that were weighted averaged into each voxel in matrix 1
%%  confidenceWeights2 - number of values that were weighted averaged into each voxel in matrix 2

%%outputs:
%%  out - combined weight average

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function out = weightAverageCombiner(in1,in2,weightedDistances1,weightedDistances2,confidenceWeights1,confidenceWeights2)
weightedDistances1(confidenceWeights1==0) = 1e30;
weightedDistances2(confidenceWeights2==0) = 1e30;

normMatrix = zeros(size(in1),'single');
weights1 = zeros(size(in1),'single');
weights2 = zeros(size(in1),'single');

confidenceWeights1 = single(confidenceWeights1);
confidenceWeights2 = single(confidenceWeights2);
ind1 = find(in1~=0&confidenceWeights1~=0);
ind2 = find(in2~=0&confidenceWeights2~=0);
ind = union(ind1,ind2);

normMatrix(ind) = (confidenceWeights1(ind).*(1./weightedDistances1(ind)))+(confidenceWeights2(ind).*(1./weightedDistances2(ind)));
weights1(ind1) = (confidenceWeights1(ind1).*(1./weightedDistances1(ind1)))./(normMatrix(ind1));
weights2(ind2) = (confidenceWeights2(ind2).*(1./weightedDistances2(ind2)))./(normMatrix(ind2));


outCx = weights1.*in1+weights2.*in2;
outMags = weights1.*abs(in1)+weights2.*abs(in2);
outPhases = angle(outCx);
out = outMags.*exp(1i.*outPhases);
end