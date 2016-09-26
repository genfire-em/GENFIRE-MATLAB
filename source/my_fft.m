%%  my_fft %%

%% Helper function to calculate forward fft

%% Author: AJ Pryor
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function kout = my_fft(img)
kout = fftshift(fftn((ifftshift(img))));
end
