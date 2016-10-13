%%  FourierShellCorrelate %%

%computes Fourier shell correlation (in 3D) or Fourier ring correlation (in
%2D) depending upon the size of the objects to compare
%%inputs:
%%  obj1 - first object to compare
%%  obj2 - second object, should be same size as obj1
%%  numBins - number of spatial frequency bins in which to compute average correlation
%%  pixSize - pixel size, used to display FSC with spatial frequency values that match the data

%%outputs:
%%  corrCoeffs - correlation coefficients
%%  spatialFrequency - "inverse resolution indices", equivalent to normalized spatial frequency
%%  meanIntensity - average intensity at each shell

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [corrCoeffs, spatialFrequency, meanIntensity] = FourierShellCorrelate(obj1,obj2,numBins,pixSize)

%check array sizes match
try sum(size(obj1)==size(obj2));
catch
    error('GENFIRE: ERROR! Objects cannot have mismatching number of dimensions');
end

if sum(size(obj1)==size(obj2)) ~= ndims(obj1)
    error('GENFIRE: ERROR! Objects cannot be different sizes');
end

%calculate FFT
k1 = my_fft(obj1);
k2 = my_fft(obj2);

%get indices of Fourier points
Q = make_Kspace_indices(obj1);

spatialFrequency = linspace(0,1,numBins+1);%compute spatial frequency bins
% spatialFrequency = linspace(0,1,numBins);%compute spatial frequency bins

%initialize ouput vectors
corrCoeffs = zeros(1,length(spatialFrequency)-1);
meanIntensity = zeros(1,length(spatialFrequency)-1);

%loop over shells
for i = 1:length(spatialFrequency)-1
    ind = (Q>=spatialFrequency(i)&Q<spatialFrequency(i+1));%indices of reciprocal voxels within frequency range of interest
    normC = sqrt(sum((abs(k1(ind)).^2)).*sum((abs(k2(ind)).^2))); %denominator of FSC, see http://en.wikipedia.org/wiki/Fourier_shell_correlation
    corrVals = (sum(k1(ind).*conj(k2(ind))))./normC; %FSC value
    corrCoeffs(i) = real(corrVals); %take real part, the imaginary part should be essentially 0
    meanIntensity(i) = mean(abs(k1(ind))+abs(k2(ind))); %average intensity within this frequency
end
spatialFrequency = linspace(0,1,numBins);%compute spatial frequency bins
% spatialFrequency(1) = [];
% spacing = spatialFrequency(2)-spatialFrequency(1);halfSpacing = spacing./2;
% spatialFrequency = spatialFrequency + halfSpacing; %%center bins
if nargin>3 %if user provided a pixel size, compute actual spatial frequency values
 maxInverseResolution = (1./(2*pixSize));
 spatialFrequency = spatialFrequency.*maxInverseResolution;
end
end