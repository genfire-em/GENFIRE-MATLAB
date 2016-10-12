%%  hermitianSymmetrize_v2 %%

%%Applies Hermitian symmetry so that voxels are equal to complex
%%conjugate of their Hermitian counterparts. If a value only exists at one
%%of the two symmetry mates, that value is used for both voxels. If values
%%exist at both positions, their average is used for the final value. This
%%implementation handles identification of symmetry mates by reading the
%%input matrix in reverse.

%%inputs:
%%  in  -  input object to which apply Hermitian symmetry

%%outputs:
%%  out  -  input object to which apply Hermitian symmetry


%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2016. All Rights Reserved.

function out = hermitianSymmetrize_v2(in)
%get initial sizes
[startdimx,startdimy,startdimz] = size(in);

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
input = zeros(dimx,dimy,dimz);

%copy input
input(1:startdimx,1:startdimy,1:startdimz) = in;

%initialize output and numberOfValues (to account for possible averaging)
out = zeros(dimx,dimy,dimz);
numberOfValues = out;

%combine with symmetry mates by reversing elements and complex conjugating
out(:,:,:) = input(:,:,:) + conj(input(end:-1:1,end:-1:1,end:-1:1));
numberOfValues(:,:,:) = single(input(:,:,:)~=0);
numberOfValues(:,:,:) = numberOfValues(:,:,:) + numberOfValues(end:-1:1,end:-1:1,end:-1:1);

out(numberOfValues~=0) = out(numberOfValues~=0)./ numberOfValues(numberOfValues~=0);

%retake original array size
out = out(1:startdimx,1:startdimy,1:startdimz);

end




