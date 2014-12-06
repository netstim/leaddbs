function I_temp=ea_edgedetection(w,alfa)
%%%%%%%%%%%%% The main.m file  %%%%%%%%%%%%%%%

% The algorithm parameters:
% 1. Parameters of edge detecting filters:
%    X-axis direction filter:
    % Nx1=10;Sigmax1=1;Nx2=10;Sigmax2=1;Theta1=pi/2;

    Nx1=10;Sigmax1=1;Nx2=10;Sigmax2=1;Theta1=pi/2;
%    Y-axis direction filter:
%     Ny1=10;Sigmay1=1;Ny2=10;Sigmay2=1;Theta2=0;
Ny1=10;Sigmay1=1;Ny2=10;Sigmay2=1;Theta2=0;
% 2. The thresholding parameter alfa:
   %  alfa=0.5;
     
% Get the initial image lena.gif






%subplot(3,2,1);
%imagesc(w);
%title('Image: lena.gif');

% X-axis direction edge detection
%subplot(3,2,2);
filterx=d2dgauss(Nx1,Sigmax1,Nx2,Sigmax2,Theta1);
Ix= conv2(w,filterx,'same');
%imagesc(Ix);
%title('Ix');

% Y-axis direction edge detection
%subplot(3,2,3)
filtery=d2dgauss(Ny1,Sigmay1,Ny2,Sigmay2,Theta2);
Iy=conv2(w,filtery,'same'); 
%imagesc(Iy);
%title('Iy');

% Norm of the gradient (Combining the X and Y directional derivatives)
%subplot(3,2,4);
NVI=sqrt(Ix.*Ix+Iy.*Iy);
%imagesc(NVI);
%title('Norm of Gradient');

% Thresholding
I_max=max(max(NVI));
I_min=min(min(NVI));
level=alfa*(I_max-I_min)+I_min;
%subplot(3,2,5);
Ibw=max(NVI,level.*ones(size(NVI)));
%imagesc(Ibw);
%title('After Thresholding');

% Thinning (Using interpolation to find the pixels where the norms of 
% gradient are local maximum.)
%subplot(3,2,6);
[n,m]=size(Ibw);
for i=2:n-1,
for j=2:m-1,
	if Ibw(i,j) > level,
	X=[-1,0,+1;-1,0,+1;-1,0,+1];
	Y=[-1,-1,-1;0,0,0;+1,+1,+1];
	Z=[Ibw(i-1,j-1),Ibw(i-1,j),Ibw(i-1,j+1);
	   Ibw(i,j-1),Ibw(i,j),Ibw(i,j+1);
	   Ibw(i+1,j-1),Ibw(i+1,j),Ibw(i+1,j+1)];
	XI=[Ix(i,j)/NVI(i,j), -Ix(i,j)/NVI(i,j)];
	YI=[Iy(i,j)/NVI(i,j), -Iy(i,j)/NVI(i,j)];
	try
    ZI=interp2(X,Y,Z,XI,YI);
    end
		if Ibw(i,j) >= ZI(1) & Ibw(i,j) >= ZI(2)
		I_temp(i,j)=I_max;
		else
		I_temp(i,j)=I_min;
		end
	else
	I_temp(i,j)=I_min;
	end
end
end


%imagesc(I_temp);
%title('After Thinning');
%colormap(gray);
%%%%%%%%%%%%%% End of the main.m file %%%%%%%%%%%%%%%


%%%%%%% The functions used in the main.m file %%%%%%%
% Function "d2dgauss.m":
% This function returns a 2D edge detector (first order derivative
% of 2D Gaussian function) with size n1*n2; theta is the angle that
% the detector rotated counter clockwise; and sigma1 and sigma2 are the
% standard deviation of the gaussian functions.
function h = d2dgauss(n1,sigma1,n2,sigma2,theta)
r=[cos(theta) -sin(theta);
   sin(theta)  cos(theta)];
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
%%%%%%%%%%%%%% end of the functions %%%%%%%%%%%%%