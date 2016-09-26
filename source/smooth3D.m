%%  smooth3D %%

%% Smooth object with a low pass filter

%%inputs:
%%  img - object to smooth
%%  resolutionCutoff - resolution for 1 sigma falloff of Gaussian smoothing filter

%%outputs:
%%  smoothImg - smoothed object
%%  kfilter  -  smoothing filter


%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function [smoothImg kfilter]= smooth3D(img,resolutionCutoff)
    %construct low pass filter and direct multiply in reciprocal space to
    %smooth
    Rsize = size(img,1);
    Csize = size(img,2);
    Lsize = size(img,3);
    Rcenter = round((Rsize+1)/2);
    Ccenter = round((Csize+1)/2);
    Lcenter = round((Lsize+1)/2);
    a=1:1:Rsize;
    b=1:1:Csize;
    c=1:1:Lsize;
    [bb,aa,cc]=meshgrid(b,a,c);
    sigma = Rsize/2*resolutionCutoff;
    kfilter=exp( -( ( ((sqrt((aa-Rcenter).^2+(bb-Ccenter).^2 + (cc-Lcenter).^2).^2) ) ./ (2* sigma.^2) )));
    kfilter=kfilter/max(kfilter(:));

    kbinned = my_fft(img);      
    kbinned = kbinned.*kfilter;
    smoothImg = my_ifft(kbinned);
end