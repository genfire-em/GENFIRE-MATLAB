%%  alignByCOMx %%

%%align 2D input along x direction using center of mass
%%inputs:
%%  in - input 2D image

%%outputs:
%%  imgOut - aligned image
%%  shiftX  -  applied shift to align image

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [imgOut, shiftX] = alignByCOMx(in)

[dimx, ~] = size(in);
centerX = round((dimx+1)/2);

COMx = round(sum(sum(in,2).*(1:dimx)')./sum(sum(in)));
shiftX = centerX-COMx;
imgOut = circshift(in,[shiftX 0  0]);
