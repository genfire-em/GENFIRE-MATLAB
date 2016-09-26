function [percentFilled, spatialFrequency] = percentageFourierGridFilledIn(measuredK,numBins,pixSize)
n2 = size(measuredK,2)/2;%array radius

if mod(size(measuredK,1),2)==0
ncK1 = size(measuredK,1)/2+1;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1-1)./n2K1;
elseif size(measuredK,1)==1
vec1 = 0;
else
ncK1 = (size(measuredK,1)+1)/2;%central pixel
n2K1 = ncK1-1;%max radius
vec1 = (-n2K1:n2K1)./n2K1; 
end

if  mod(size(measuredK,2),2)==0
ncK2 = size(measuredK,2)/2+1;%central pixel
n2K2 = ncK2-1;%max radius
vec2 = (-n2K2:n2K2-1)./n2K2;
elseif size(measuredK,2)==1
vec2 = 0;
else
ncK2 = (size(measuredK,2)+1)/2;%central pixel
n2K2 = ncK2-1;%max radius
vec2 = (-n2K2:n2K2)./n2K2; 
end

if  mod(size(measuredK,3),2)==0
ncK3 = size(measuredK,3)/2+1;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3-1)./n2K3;
elseif size(measuredK,3)==1
vec3 = 0;
else
ncK3 = (size(measuredK,3)+1)/2;%central pixel
n2K3 = ncK3-1;%max radius
vec3 = (-n2K3:n2K3)./n2K3; 
end
[Kx Ky Kz] = meshgrid(vec2,vec1,vec3);%grid of Fourier indices
Kmags = sqrt(Kx.^2+Ky.^2+Kz.^2);%take magnitude


spatialFrequency = linspace(0,1,numBins+1);%compute spatial frequency bins
for i = 1:length(spatialFrequency)-1
    ind = (Kmags>=spatialFrequency(i)&Kmags<spatialFrequency(i+1));%indices of reciprocal voxels within frequency range of interest
    percentFilled(i) = sum(sum(sum(measuredK(ind)~=0)))./sum(sum(sum(ind)));
end
spatialFrequency(end) = [];
spacing = spatialFrequency(2)-spatialFrequency(1);halfSpacing = spacing./2;
spatialFrequency = spatialFrequency + halfSpacing; %%center bins
if nargin>2 %if user provided a pixel size, compute actual spatial frequency values
 maxInverseResolution = (1./pixSize)/2;
 spatialFrequency = spatialFrequency.*maxInverseResolution;
end
end