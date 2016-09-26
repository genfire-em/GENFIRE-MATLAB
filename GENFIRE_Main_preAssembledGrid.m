%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                         %%
%%                        Welcome to GENFIRE!                              %%
%%           GENeralized Fourier Iterative REconstruction                  %%
%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENFIRE is a robust, Fourier-based reconstruction algorithm that is
%% capable of using a limited set of input projections to generate a 3D reconstruction
%% while also partially retrieving missing projection information. It does this by iterating 
%% between real and reciprocal space and applying simple constraints in each to find an optimal solution that  
%% correlates with the input projections while simultaneously obeying real space conditions
%% such as positivity (the assumption that there is no negative electron density)
%% and support (we require that the 3D object resulting from zero padded projections
%% exists entirely within some smaller region). The result is a more consistent 
%% and faithful reconstruction with superior contrast and, in many cases, resolution when
%% compared with more traditional 3D reconstruction algorithms such as
%% Filtered Back-Projection. In addition, GENFIRE is capable of dynamically identifying
%% the noise level of the input projections (out to what resolution does useful information exist
%% before noise dominates the signal) without any prior knowledge from the user. 

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a modified version of GENFIRE that is designed to take in a pre-assembled
%% grid of Fourier constraints, such as obtained through MRI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath ./source/
%%%   User Parameters   %%%

filename_FourierGrid = './data/preAssembledGrid.mat'; %%filename of pre-assembled Fourier constraint grid


filename_Support = './data/support180.mat'; %% NxNxN binary array specifying a region of 1's in which the reconstruction can exist 

%filename_InitialModel = '';

filename_results = './results/GENFIRE_rec_preAssembledGrid.mat';


global numIterations 
numIterations = 100; 

global pixelSize
pixelSize = 4; 

numBins = 50; %number of bins for FSC averaging

%%%   Begin Reconstruction   %%%
clc
measuredK = importdata(filename_FourierGrid);


if exist('filename_InitialModel','var')
    initialObject = single(importdata(filename_InitialModel));
else
    initialObject = [];
end


global support %make support variable globally accessable to avoid passing copies of large arrays around to different functions
support = single(importdata(filename_Support));

Q = make_Kspace_indices(support);

recIFFT = real(my_ifft(measuredK));%take back original sized array of the initial IFFT to compare before/after iteration
tic


fprintf('GENFIRE: Reconstructing... \n\n');
if isempty(initialObject) 
    [GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,recIFFT,support,measuredK,ones(size(measuredK)),1);%run iterations
else
    [GENFIRE_rec, errK] = GENFIRE_iterate(numIterations,recIFFT,support,measuredK,ones(size(measuredK)),1);

end

reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('GENFIRE: Reconstruction completed in %.12g seconds\n\n',reconstructionTime);

%display results
figure,
subplot(2,3,4), imagesc(squeeze(sum(GENFIRE_rec,1))),title('GENFIRE projection 1')
subplot(2,3,5), imagesc(squeeze(sum(GENFIRE_rec,2))),title('GENFIRE projection 2')
subplot(2,3,6), imagesc(squeeze(sum(GENFIRE_rec,3))),title('GENFIRE projection 3')
subplot(2,3,1), imagesc(squeeze(sum(recIFFT,1))),title('before iteration projection 1')
subplot(2,3,2), imagesc(squeeze(sum(recIFFT,2))),title('before iteration projection 2')
subplot(2,3,3), imagesc(squeeze(sum(recIFFT,3))),title('before iteration projection 3')

nc = round((size(GENFIRE_rec,2)+1)/2);
figure,
subplot(2,3,4), imagesc(squeeze(GENFIRE_rec(ncK1,:,:))),title('GENFIRE slice 1')
subplot(2,3,5), imagesc(squeeze(GENFIRE_rec(:,ncK2,:))),title('GENFIRE slice 2')
subplot(2,3,6), imagesc(squeeze(GENFIRE_rec(:,:,ncK3))),title('GENFIRE slice 3')
subplot(2,3,1), imagesc(squeeze(recIFFT(ncK1,:,:))),title('before iteration slice 1')
subplot(2,3,2), imagesc(squeeze(recIFFT(:,ncK2,:))),title('before iteration slice 2')
subplot(2,3,3), imagesc(squeeze(recIFFT(:,:,ncK3))),title('before iteration slice 3')

%save results
save(filename_results,'GENFIRE_rec','errK')














