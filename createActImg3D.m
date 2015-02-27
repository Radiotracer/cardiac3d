function [ img ] = createActImg3D(p)
% This function creates a  3D image simulating the left
% ventricle (LV).
% Pre-defined Constants:
% 1. nApex=4: Number of sampling radials of the apex in long axis slice
% (22.5, 45, 67.5, 90)
% 2. nBase=5: Number of sampling slices in the base 
% 2. nSeg=8: Number of sampling radials in short axis slice. 
% 
%  Parameters: 
% a. (p(1),p(2)): Define the axis parallel to the z-axis
% b. p(3): z coordinate of the apex center, which is on the axis
% c. p(4): Sampling interval of the z-axis in the base
% % Endocardium radius and myocardium thickness come in pair
% Apex
% d. (p(5),p(6)): Apex radius and myocardium thickness on the axis
% e. (p(7),p(8)) -> ( p(7+2*(nSeg-1)), p(8+2*(nSeg-1)):     90-90/nApex
%    (p(7+2*nSeg),p(8+2*nSeg)) ->( p(7+2*nSeg+2*(nSeg-1)), p(8+2*nSeg+2*(nSeg-1)):    90 -90/nApex *2
%    ........
%    p(7+k*2*nSeg),p(8+k*2*nSeg)) ->( p(7+k*2*nSeg+2*(nSeg-1)), p(8+k*2*nSeg+2*(nSeg-1)):  0<= k <= (nApex-1)
% Base radii
% f: ( p(7+2*nApex*nSeg+k*2*nSeg),p(8+2*nApex*nSeg+k*2*nSeg))
%          --> ( p(7+2*nApex*nSeg+k*2*nSeg+2*(nSeg-1)),p(8+2*nApex*nSeg+k*2*nSeg+2*(nSeg-1))): 0<= k <=(nBase-1)
% 
% p=[63, 64,40,10,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...  % end of apex
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,... % start of base
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8,...
%         25,8,25,8,25,8,25,8,25,8,25,8,25,8,25,8];  
dimX=128;
dimY=128;
dimZ=100;    
nApex=4;
nBase=5;
nSeg=8;
dsAng=2*pi/nSeg;
daAng=(0.5*pi)/nApex;

% global dimX;
% global dimY;
% global dimZ;    
% global nApex;
% global nBase;
% global nSeg;
% global dsAng;
% global daAng;
invx=zeros((1+nApex+nBase),nSeg+1);
invy=zeros((1+nApex+nBase),nSeg+1);
invz=zeros((1+nApex+nBase),nSeg+1);
invx(1,:)=p(1);
invy(1,:)=p(2);
invz(1,:)=p(3)-p(5);

outvx=zeros((1+nApex+nBase),nSeg+1);
outvy=zeros((1+nApex+nBase),nSeg+1);
outvz=zeros((1+nApex+nBase),nSeg+1);
outvx(1,:)=p(1);
outvy(1,:)=p(2);
outvz(1,:)=p(3)-p(5)-p(6);

for m=1:nApex
    for n=1:nSeg
        idx=7+2*((m-1)*nSeg+ (n-1));
        angX=sin(daAng*m) * cos((n-1)*dsAng);
        angY=sin(daAng*m) * sin((n-1)*dsAng);
        invx(1+m,n)=p(1)+ p(idx) * angX;
        invy(1+m,n)=p(2)+ p(idx) * angY; 
        outvx(1+m,n)=p(1)+ (p(idx)+p(idx+1)) * angX;
        outvy(1+m,n)=p(2)+ (p(idx)+p(idx+1)) * angY;    
    end
    invz(1+m,:)=p(3)-p(5)*cos(daAng*m);
    outvz(1+m,:)=p(3)-(p(5)+p(6))*cos(daAng*m);
end

for m=1:nBase
    for n=1:nSeg
        idx=7+2*nApex*nSeg+2*((m-1)*nSeg+ (n-1));
        angX=cos((n-1)*dsAng);
        angY=sin((n-1)*dsAng);
        invx(1+nApex+m,n)=p(1)+p(idx)*angX;
        invy(1+nApex+m,n)=p(2)+p(idx)*angY;
        outvx(1+nApex+m,n)=p(1)+(p(idx)+p(idx+1))*angX;
        outvy(1+nApex+m,n)=p(2)+(p(idx)+p(idx+1))*angY;        
    end
    invz(1+nApex+m,:)=p(3)+m*p(4);
    outvz(1+nApex+m,:)=p(3)+m*p(4);
end
invx(:,end)=invx(:,1);
invy(:,end)=invy(:,1);
invz(:,end)=invz(:,1);

outvx(:,end)=outvx(:,1);
outvy(:,end)=outvy(:,1);
outvz(:,end)=outvz(:,1);

% figure;surf(invx,invy,invz);axis equal;
% xlabel('x');ylabel('y');zlabel('z');
% figure;surf(outvx,outvy,outvz);axis equal;
% xlabel('x');ylabel('y');zlabel('z');

inv=cat(3,invx,invy,invz);
inv=permute(inv,[3,1,2]);
outv=cat(3,outvx,outvy,outvz);
outv=permute(outv,[3,1,2]);


i=1:(1+nApex+nBase); 
j=1:(nSeg+1);        

incs = csape({i,j},inv,{'clamped','periodic'}); %spline interpolation with end conditions 
innz=floor(max(invz(:))-min(invz(:)))+1;
inPts=fnval(incs,{linspace(1,1+nApex+nBase,innz),linspace(1,nSeg+1,51)});

outcs = csape({i,j},outv,{'clamped','periodic'}); 
outnz=floor(max(outvz(:))-min(outvz(:)))+1;
outPts=fnval(outcs,{linspace(1,1+nApex+nBase,outnz),linspace(1,nSeg+1,51)});

inMask=[];
for k=1:innz
    mask=poly2mask(squeeze(inPts(1,k,:)),squeeze(inPts(2,k,:)),dimX,dimY);
    inMask=cat(3,inMask,mask);
end
inMask=cat(3,zeros(dimX,dimY,floor(p(3)-p(5))),inMask);
inMask=cat(3,inMask,zeros(dimX,dimY,dimZ-size(inMask,3)));

outMask=[];
for k=1:outnz
    mask=poly2mask(squeeze(outPts(1,k,:)),squeeze(outPts(2,k,:)),dimX,dimY);
    outMask=cat(3,outMask,mask);
end
outMask=cat(3,zeros(dimX,dimY,floor(p(3)-p(5)-p(6))),outMask);
outMask=cat(3,outMask,zeros(dimX,dimY,dimZ-size(outMask,3)));

btMask=outMask-inMask;
img=btMask;

end

