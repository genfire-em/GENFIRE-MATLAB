%%  CrossCorrelate %%

%computes n-D cross correlation that is same size as input objects
%shiftX, shiftY, and shiftZ are the shifts to be applied to obj2 to match
%it to obj1
%%inputs:
%%  obj1,obj2 - matrices to correlate

%%outputs:
%%  xcorr - cross correlation result
%%  shiftX(Y)(Z) - shift amount in pixels that will move obj2 such that 
%%      a recalculation of the cross correlation would yield the maximum at the center

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.


function [xcorr, shiftX, shiftY, shiftZ] = CrossCorrelate(obj1,obj2)

%check array sizes match
try sum(size(obj1)==size(obj2));
catch
    error('GENFIRE: ERROR! Objects cannot have mismatching number of dimensions');
end

if sum(size(obj1)==size(obj2)) ~= ndims(obj1)
    error('GENFIRE: ERROR! Objects cannot be different sizes');
end

%get initial sizes
[startdimx,startdimy,startdimz] = size(obj1);

%copy sizes because they may be modified momentarily
dimx = startdimx;
dimy = startdimy;
dimz = startdimz;

%If any dimensions are even, remake that dimension to be 
%one pixel larger, and copy the input array.
if mod(startdimx,2)==0
    dimx = dimx+1;
end

if mod(startdimy,2)==0
    dimy = dimy+1;
end

if mod(startdimz,2)==0
    dimz = dimz+1;
end
%initialize input with all odd dimensions
input1 = zeros(dimx,dimy,dimz);

%copy input
input1(1:startdimx,1:startdimy,1:startdimz) = obj1;clear obj1

%initialize input with all odd dimensions
input2 = zeros(dimx,dimy,dimz);

%copy input
input2(1:startdimx,1:startdimy,1:startdimz) = obj2;clear obj2

xcorr = my_ifft( (my_fft(input1)) .* (my_fft(conj(input2(end:-1:1,end:-1:1,end:-1:1)))) );

if nargin>1 %compute shifts
    
    ncX = round((size(input1,1)+1)/2);
    ncY = round((size(input1,2)+1)/2);
    ncZ = round((size(input1,3)+1)/2);

    [Xmax, Ymax, Zmax] = ind2sub(size(input1),find(xcorr==max(xcorr(:))));
    shiftX = (Xmax - ncX);
    shiftY = (Ymax - ncY);
    shiftZ = (Zmax - ncZ);
end
