clear; close all;
global dimX;dimX=128;
global dimY;dimY=128;
global dimZ;dimZ=100;    
global nApex;nApex=4;
global nBase;nBase=5;
global nSeg;nSeg=8;
global dsAng;dsAng=2*pi/nSeg;
global daAng;daAng=(0.5*pi)/nApex;
global imgMd;
tp=[63, 64,40,10,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...  % end of apex
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,... % start of base
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8];  
imgMd=createActImg3D(tp);
figure;imshow3D(imgMd,[]);

initP=[63, 64,40,10,22,13,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...  % end of apex
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,... % start of base
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
        25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8];  
imgInit=createActImg3D(initP);
figure;imshow3D(imgInit,[]);

options = optimoptions(@fmincon,...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'FinDiffRelStep',0.001,...
    'MaxFunEvals',10000 ...
    );

tic
[ncep,ncfval] = fmincon(@diffFunc,initP,[],[],[],[],[],[],@noconstraint,options);
toc