% Checks if the input argument is a vector.
%		truth = mn_isvector(inVc)
%
%		Returns 1 if inVc is a vector, else 0.
%		Single numbers (length(vec)==1) are understood as not beeing a vector.
%		Function accepts both, row and  column vectors.
%
%		Examples:
%		truth = mn_isvector([1;2;5]);		-> truth=1
%		truth = mn_isvector(1);			-> truth=0
%		truth = mn_isvector('string');	-> truth=0
%
%		Harald Fischer
%		1/00
%
%		PC

% ----------------------- development -----------------------------------
%
% Authors:
%   Björn Kreher (BK)
%   Michael von Mengershausen (MVM)
%
% Projects:
%   BK: compatibility of matlab_new 'isvector' function with builtin
%       'isvector' function
%
%   MVM: renaming of matlab_new 'isvector' --> 'mn_isvector' (060110)
    


function truth = mn_isvector(inVc)

if exist('isvector', 'builtin') 
    truth= builtin('isvector', inVc);
    return
end

%%%%% init and error check
truth = 0;
if nargin~=1
   warning('wrong number of input arguments');
   return;
end
if ~isnumeric(inVc)
   warning('input arguments is not numeric');
   return;
end
%%%%% End of: init and error check


if ndims(inVc)==2
   [size_y, size_x] = size(inVc);
   if (size_x==1 & size_y>1) | (size_y==1 & size_x>1)
      truth = 1;
   end
end
% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics
