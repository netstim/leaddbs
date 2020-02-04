function nx=ea_vecnorm(x,p)
% Vector / matrix norm.
% This is a work-around for missing `vecnorm` function prior Matlab 2017b.

if nargin<2
    p=2;
end

if isvector(x)
    nx=norm(x,p);
else
    nx=x;
    for i=1:size(x,2)
        nx(:,i)=norm(x(:,i),p);
    end
end
