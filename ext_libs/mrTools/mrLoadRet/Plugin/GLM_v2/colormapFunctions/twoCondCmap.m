function cmap = twoCondCmap(numGrays)
  minx = 0;maxx = 0.5;
  x = minx:(maxx-minx)/(numGrays/2):maxx;
  y = mygauss([1 0 0.25 0],x);
  
  red = [zeros(length(y),1) y' y'];
  y = fliplr(y);
  blue = [y' y' zeros(length(y),1)];
  
  cmap = 1-[red;blue];
return

% Gaussian
%
% usage: gauss(p,X,Y);
%   p is an array of parameters:
%     p(1) = height of Gaussian
%     p(2) = center x
%     p(3) = center y
%     p(4) = SD in x dimension
%     p(5) = SD in y dimension
%     p(6) = offset
%     p(7) = rotation in radians
%   X and Y are the position on which to evaluate the gaussian
%     to evaluate at a matrix of points use,
%     e.g. [X,Y] = meshgrid(-1:.1:1,-1:.1:1);
%
%  the function can also be called as follows for 1D
%  usage: gauss(p,X);
%
%     p(1) = height of Gaussian
%     p(2) = center
%     p(3) = SD
%     p(4) = offset
% 
%   by: justin gardner
% date: 6/6/97
function G=mygauss(p,X,Y)

% 2D Gaussian 
if nargin == 3

  % rotate coordinates
  % note that the negative sign is because
  % we are rotating the coordinates and not the function
  X1 = cos(-p(7)).*X - sin(-p(7)).*Y;
  Y1 = sin(-p(7)).*X + cos(-p(7)).*Y;

  % calculate the Gaussian
  G = p(1) * exp(-((((X1-p(2)).^2)/(2*p(4)^2))+(((Y1-p(3)).^2)/(2*p(5)^2))))+p(6);
  
% 1D Gaussian
elseif nargin == 2

  % calculate the Gaussian
  G = p(1) * exp(-(((X-p(2)).^2)/(2*p(3)^2)))+p(4);
  
else 
   % usage error
   disp('USAGE: gauss(parameters, X, Y)');
end
