%%  smooth3D %%

%% Smooth object with a low pass filter

%%inputs:
%%  img - object to smooth
%%  resolutionCutoff - resolution for 1 sigma falloff of Gaussian smoothing filter

%%outputs:
%%  smoothImg - smoothed object
%%  kfilter  -  smoothing filter


%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [smoothImg kfilter]= smooth3D(img,resolutionCutoff)
    %construct low pass filter and direct multiply in reciprocal space to
    %smooth
    [dimx, dimy, dimz] = size(img);
    ncx = round((dimx+1)/2);
    ncy = round((dimy+1)/2);
    ncz = round((dimz+1)/2);
    a=1:dimx;
    b=1:dimy;
    c=1:dimz;
    [bb,aa,cc]=meshgrid(b,a,c);
    sigma = dimx/2*resolutionCutoff;
    kfilter=exp( -( ( ((sqrt((aa-ncx).^2+(bb-ncy).^2 + (cc-ncz).^2).^2) ) ./ (2* sigma.^2) )));
    kfilter=kfilter/max(kfilter(:));
    kbinned = my_fft(img);      
    kbinned = kbinned.*kfilter;
    smoothImg = my_ifft(kbinned);
end