clear;close all;
lv3d=readim('lv.sa.im');
figure;imshow3D(lv3d,[]);

p=[63, 64,40,10,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...  % end of apex
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,... % start of base
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8
        ,25,8];  
    
tic
img= createActImg3D(p);
toc
figure;imshow3D(img,[]);


p1=[63, 64,40,10,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...  % end of apex
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,... % start of base
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8];  
tic
img1= createActImg3Dv1(p1);
toc
figure;imshow3D(img,[]);


p=[63, 64,40,10,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...  % end of apex
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,... % start of base
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8];  
 p2=p/2;
tic
img2= createActImg3Dv2(p2);
toc
figure;imshow3D(img2,[]);


tp_w=[128.5 128.5 15 5 ...
    50 20 50 20 50 20 50 20 50 20 50 20 50 20 50 20 ...
    60 65 70 75 80 85 90 95 100 105 110 115 120 125  130 135 140 145];
global gaussFilter;
gaussFilter= fspecial('gaussian', [29 29], 6.37);
global nRad; nRad=8;
global rAng; rAng=2*pi/nRad;
global hrAng; hrAng=pi/nRad;
global nSeg;nSeg=floor((length(tp_w)-4-2*nRad));
global sAng;sAng=2*pi/nSeg;
global dimX;dimX=256;
global dimY;dimY=256;
tic
img_w=createActImg2D(tp_w);
toc
figure;imshow(img_w,[]);title('Truth(Weight)');






p=[63, 64,40,10,25,8,...
        25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,...  % end of apex
        25,8,25,8,25,8,25,8,... % start of base
        25,8,25,8,25,8,25,8];  
 p3=p/2;
tic
img3= createActImg3Dv3(p3);
toc
figure;imshow3D(img3,[]);

