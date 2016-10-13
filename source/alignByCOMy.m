%%  alignByCOMy %%

%%align 2D input along y direction using center of mass
%%inputs:
%%  in - input 2D image

%%outputs:
%%  imgOut - aligned image
%%  shiftY  -  applied shift to align image

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [imgOut, shiftY] = alignByCOMy(in)

[~, dimy] = size(in);
centerY = round((dimy+1)/2);
 
COMy = round(sum(sum(in,1).*(1:dimy))./sum(sum(in)));
shiftY = centerY-COMy;

imgOut = circshift(in,[0 shiftY 0]);