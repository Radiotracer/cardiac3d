function [mse] = diffFunc(p)
%  This function measures the difference between the measured image and the
%  modeled one. 
%  p--Model parameters
global imgMd;
img=createActImg3Dv3( p );
tmp=abs(img-imgMd);
mse=sum(tmp(:));
end

