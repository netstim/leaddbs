function [res, errStr]= my_mod(x, y)

% overloded function:
%  MOD    Modulus (signed remainder after division).
%     MOD(x,y) is x - y.*floor(x./y) if y ~= 0.  By convention, MOD(x,0) is x.
%     The input x and y must be real arrays of the same size, or real scalars.
%  
%     The statement "x and y are congruent mod m" means mod(x,m) == mod(y,m).
%  
%     MOD(x,y) has the same sign as y while REM(x,y) has the same sign as x.
%     MOD(x,y) and REM(x,y) are equal if x and y have the same sign, but
%     differ by y if x and y have different signs.
%  
%     See also REM.
%
%  
% new feature: mod(i*y, y) == y instead of 0
%
%
% Bjoern W. Kreher
% 07/02
%
% UNIX

res= mod(x, y);
comp= res == 0;
res= (1 - comp).*res + comp.*y;