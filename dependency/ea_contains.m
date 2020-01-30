function rv=ea_contains(str,pattern)
% Determine if pattern is in string(s).
%
% NOTE: This is a work-around for missing `contains()` function in
% Matlab version prior R2016b. Please note that it is not 100%
% compatible with the Matlab `contains()` function, it merely provides
% basic functionality. Specially, it does not support case-insensitive
% search and multiple-string patterns.
% It is based on the `strfind()` function.
% See https://www.mathworks.com/help/matlab/ref/strfind.html
% and https://www.mathworks.com/help/matlab/ref/contains.html
% for details.
%
% Arguments:
%   str -     character vector, cell array of character vectors, or string
%             array
%   pattern - character vector or string scalar
%
% Returns a logical value (or a vector of logical values) describing the
% presence of given pattern in given string(s).
%
% Examples:
%
%   ea_contains('blah foo','foo')
%       returns 1
%
%   ea_contains({'blah','foo'},'f')
%       returns [0, 1]
%
% Author: T. Sieger, tomas.sieger@seznam.cz
%

if nargin<2
    error('Not enough input arguments.');
end

x=strfind(str,pattern);
n=length(x);
if iscell(x) && n>0
    rv=repmat(logical(0),1,n);
    for i=1:n
        rv(i)=~isempty(x{i});
    end
else
    rv=~isempty(x);
end
