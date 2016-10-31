% indexMax.m: returns index of maximum value (and maximum value as second output) 
%
%        $Id$
%
%   example: if x=[4 5], y=[6 4] and z=[5 6], indexMax(x,y,z) returns 
%             index=[2 3] and max=[6 6] (because y(1)>x(1), y(1)>z(1) and z(2)>x(2), z(2)>y(2) )
%           

function [index, maximum] = indexMax(varargin)


nDims = length(size(varargin{1}));
array = varargin{1};
for i = 2:nargin
  %first check that all inputs have the same size
  if ~isequal(size(varargin{i}),size(varargin{1}))
    error('All inputs must have the same size')
  else
    %concatenate
    array=cat(nDims+1,array,varargin{i});
  end
end

[maximum,index] = max(array,[],nDims+1);

%put NaNs if all values are NaNs
index(all(isnan(array),nDims+1))=NaN;

