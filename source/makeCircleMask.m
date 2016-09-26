function out = makeCircleMask(radius,imgSize);

out = zeros(imgSize,imgSize);
nc = imgSize/2+1;
n2 = nc-1;
[xx yy] = meshgrid(-n2:n2-1,-n2:n2-1);
R = sqrt(xx.^2 + yy.^2);
out = R<=radius;