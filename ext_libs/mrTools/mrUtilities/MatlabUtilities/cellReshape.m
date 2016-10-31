function reshapedCellArray = cellReshape(inputCellArray,varargin)
%cellReshape equivalent of reshape for cell arrays
%
%     Syntax
%        B = cellReshape(A,m,n,p,...)
%        B = cellReshape(A,[m n p ...])
%     reshape cell array A into a new cell array B of dimension m * n * p ....
%   
%     Author: Julien Besle, 07/07/2010
%        $Id$


outputDims = reshape(cell2mat(varargin),1,[]);
 

if prod(outputDims) ~= numel(inputCellArray)
   error('To RESHAPE the number of elements must not change.')
end

inputCellArray = inputCellArray(:);

reshapedCellArray = cell(outputDims);
for iCell = 1:numel(reshapedCellArray)
   reshapedCellArray(iCell) = inputCellArray(iCell);
end

end

