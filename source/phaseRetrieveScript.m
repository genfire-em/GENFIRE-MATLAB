%%script to demonstrate GENFIRE with OSS
% obj = padarray(phantom(64),[64*2 64*2]);
% obj = phantom(64)+rand(64,64);
obj = phantom(64);
% obj = rand(64,64);

obj = padarray(obj,[64*4 64*4]);

k = my_fft(obj);
measuredK = abs(k);
numIterations = 500;
initialObject = zeros(size(k));
support = ones(64,64);
support(1,1) = 0;
support = obj>0;
% support = convn(support,ones(3,3,3),'same')>0;
paddingX= size(obj,2)-size(support,2);
paddingY= size(obj,1)-size(support,1);

support =padarray(support,[paddingY/2 paddingX/2]);
constraintInd_complex = [];
constraintInd_magnitude = find(measuredK~=0);
R_freeInd_mag = [];

[rec, errK, Rfree_magnitude,bestErrors] = GENFIRE_OSS(numIterations,initialObject,...
support,measuredK,constraintInd_complex,constraintInd_magnitude,R_freeInd_mag);
