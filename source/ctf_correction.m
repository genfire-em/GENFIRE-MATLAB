function [CTF, gamma] = ctf_correction(imin,DeltaFu,DeltaFv,azimuthal_ang,ignore_first_peak)

%% Initial parameters
if nargin < 5
    ignore_first_peak = 0;
end

Voltage = 200; % Accelarating voltage in KeV
% pix_size = 0.6375; % pixel size in Angstrom
pix_size = 1.4; % pixel size in Angstrom

temp = imin;
nc = size(temp,2)/2+1;
n2 = nc-1;

lambda=12.2643247 / sqrt(Voltage*1e3 * (1. + Voltage *1e3* 0.978466e-6));
Qn = 0.1; % Amplitude contrast stated to work best for CTF correction of vitrious ice
Qn = 0;


Cs= 2.0*10^(7); % Spherical abberation of the microscope from mm into Angstrom
% Cc = 2*10^(7); % Chromatic abberation of the microscope from mm into Angstrom
% deltaE = 0.02; % Energy width in eV
% alph = 0.02/1000; % illumination aperture in mrad (default 0.02 - reasonable for FEG)
% % azimuthal_ang = -2.04; % angle of astygmatism in degrees
rad_azimuth  = degtorad(azimuthal_ang);
defocav = -(DeltaFu + DeltaFv)/2;
defocdev = -(DeltaFu - DeltaFv)/2;
qq = (-n2:n2-1);
[Q1 Q2] = meshgrid(qq,qq);
QQ = sqrt(Q2.^2+Q1.^2).*1./(2*pix_size)./(n2); % Convert into a radially symmetric qspace
elipsoid_ang = atan2(Q1,Q2) - rad_azimuth;
cos_elipsoid_ang_2 = cos(2*elipsoid_ang);
Dz = (cos_elipsoid_ang_2.*defocdev + defocav); % 2D array containing the defocus as a function of radius

gamma = (2*pi).*(((Cs*lambda.^3.*QQ.^4)./4) + ((Dz.*lambda.*QQ.^2)./2 ));
CTF =  -(sqrt(1-Qn^2).*sin(gamma)  - Qn.*cos( gamma)); % calculate the contrast transfer function
if ignore_first_peak == 1;
    first_peak = find(abs(gamma) < (pi/2));
    CTF(first_peak) = 1;
end

%% Perform CTF correction

% dp = my_fft(temp);
% flip_ind = find(CTF < 0);
% dp(flip_ind) = -1*dp(flip_ind);
% 
% ctf_corrected = dp;
% 
% corrected_image = real(my_ifft(ctf_corrected));
% corrected_imageTimesCTF = real(my_ifft(ctf_corrected.*abs(CTF)));

end

function out = degtorad(in)
out = pi./180.*in;
end
