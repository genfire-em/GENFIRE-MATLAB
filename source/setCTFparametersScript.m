numStacks = 8;
ignore_first_peak = 0;
throwOutThreshhold = 0.09;
phaseFlip = 0;
multiplyByCTFabs = 1;

for stackNum = 1:numStacks
    stackNum
    filename = ['CTFparameters' num2str(stackNum) '.mat'];
    CTFparameters = importdata(filename);
    for partNum = 1:size(CTFparameters,2)
        CTFparameters(partNum).ignore_first_peak = ignore_first_peak;
        CTFparameters(partNum).throwOutThreshhold = throwOutThreshhold;
        CTFparameters(partNum).phaseFlip = phaseFlip;
        CTFparameters(partNum).multiplyByCTFabs = multiplyByCTFabs;
    end
    save(filename,'CTFparameters')
end
        
    
    