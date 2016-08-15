function Y=ea_detect_edges_3d(varargin)

%% function detect_edges_3d for Nii-Images
% Andreas Horn, 2014
% in: X (image matrix), alpha (threshold), useimgtoolbox (use image toolbox
% or replacement function).

X=varargin{1};
alpha=varargin{2};
if nargin==3
    useimgtbx=varargin{3};
else
    useimgtbx=0;
end

Y=X;
dimensionality=size(X,3);
ea_dispercent(0,'Sampling');

for slice=1:dimensionality
    if useimgtbx
    edgedslice=edge(squeeze(X(:,:,slice)),'canny',alpha);
    else
    edgedslice=ea_edgedetection(squeeze(X(:,:,slice)),alpha);
    end
    Y(:,:,slice)=edgedslice;
    ea_dispercent(slice/dimensionality);
end
ea_dispercent(1,'end');


function outimg=ea_edgedetection(inimg,alpha)

Nx1=10;Sigmax1=1;Nx2=10;Sigmax2=1;Theta1=pi/2;
Ny1=10;Sigmay1=1;Ny2=10;Sigmay2=1;Theta2=0;


% X-axis direction edge detection
filterx=d2dgauss(Nx1,Sigmax1,Nx2,Sigmax2,Theta1);
convimgx= conv2(inimg,filterx,'same');

% Y-axis direction edge detection
filtery=d2dgauss(Ny1,Sigmay1,Ny2,Sigmay2,Theta2);
convimgy=conv2(inimg,filtery,'same');

convimg=sqrt(convimgx.*convimgx+convimgy.*convimgy);

% Thresholding
maxintens=max(convimg(:));
minintens=min(convimg(:));
level=alpha*(maxintens-minintens)+minintens;
thinimg=max(convimg,level.*ones(size(convimg)));



% Thinning (Using interpolation to find the pixels where the norms of
% gradient are local maximum.)
[n,m]=size(thinimg);
outimg=zeros(n-1,m-1);



for i=2:n-1,
    for j=2:m-1,
        if thinimg(i,j) > level,
            X=[-1,0,+1;-1,0,+1;-1,0,+1];
            Y=[-1,-1,-1;0,0,0;+1,+1,+1];
            Z=[thinimg(i-1,j-1),thinimg(i-1,j),thinimg(i-1,j+1);
                thinimg(i,j-1),thinimg(i,j),thinimg(i,j+1);
                thinimg(i+1,j-1),thinimg(i+1,j),thinimg(i+1,j+1)];
            XI=[convimgx(i,j)/convimg(i,j), -convimgx(i,j)/convimg(i,j)];
            YI=[convimgy(i,j)/convimg(i,j), -convimgy(i,j)/convimg(i,j)];
            try
                ZI=interp2(X,Y,Z,XI,YI);
            end
            if thinimg(i,j) >= ZI(1) && thinimg(i,j) >= ZI(2)
                outimg(i,j)=maxintens;
            else
                outimg(i,j)=minintens;
            end
        else
            outimg(i,j)=minintens;
        end
    end
end


[xx,yy]=meshgrid(1:n-1,1:m-1);

[xxq,yyq]=meshgrid(1:n,1:m);
outimg = interp2(yy',xx',outimg,yyq,xxq,'cubic')';




function h = d2dgauss(n1,sigma1,n2,sigma2,theta)
r=[cos(theta) -sin(theta);
   sin(theta)  cos(theta)];
h=zeros(n1,n2);
for i = 1 : n2
    for j = 1 : n1
        u = r * [j-(n1+1)/2 i-(n2+1)/2]';
        h(i,j) = gauss(u(1),sigma1)*dgauss(u(2),sigma2);
    end
end
h = h / sqrt(sum(sum(abs(h).*abs(h))));

% Function "gauss.m":
function y = gauss(x,std)
y = exp(-x^2/(2*std^2)) / (std*sqrt(2*pi));

% Function "dgauss.m"(first order derivative of gauss function):
function y = dgauss(x,std)
y = -x * gauss(x,std) / std^2;
