function [xx,yy,zz] = ea_singlecylinder(varargin)
%CYLINDER Generate cylinder.
%   [X,Y,Z] = CYLINDER(R,N) forms the unit cylinder based on the generator
%   curve in the vector R. Vector R contains the radius at equally
%   spaced points along the unit height of the cylinder. The cylinder
%   has N points around the circumference. SURF(X,Y,Z) displays the
%   cylinder.
%
%   [X,Y,Z] = CYLINDER(R), and [X,Y,Z] = CYLINDER default to N = 20
%   and R = [1 1].
%
%   Omitting output arguments causes the cylinder to be displayed with
%   a SURF command and no outputs to be returned.
%
%   CYLINDER(AX,...) plots into AX instead of GCA.
%
%   See also SPHERE, ELLIPSOID.

%   Clay M. Thompson 4-24-91, CBM 8-21-92.
%   Copyright 1984-2002 The MathWorks, Inc. 

% Parse possible Axes input
narginchk(0,3);
[cax,args,nargs] = axescheck(varargin{:});

n = 20;
r = [1 1]';
if nargs > 0, r = args{1}; end
if nargs > 1, n = args{2}; end
r = r(:); % Make sure r is a vector.
m = length(r); if m==1, r = [r;r]; m = 2; end

theta=linspace(0,2*pi,n);

sintheta = sin(theta); sintheta(n) = 0;

x = r * cos(theta);
y = r * sintheta;

z = (0:m-1)'/(m-1) * ones(1,n);

if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax)
else
    xx = x; yy = y; zz = z;
end
