function smoothedVals=connectionBasedSmooth(connectionMatrix,inputVals)
% smoothedVals=connectionBasedSmooth(connectionMatrix,inputVals)
%
% PURPOSE: Performs neighborhood averaging of the values in inputVals.
% connectionMatrix tells us which nodes in inputVals are connected. Those
% Each element in inputVals is replaced by the mean of itself and its neigbours.
% 
% AUTHOR: Wade
% Date written : 061603
% example : smoothedVals=connectionBasedSmooth(connectionMatrix,inputVals);

inputVals=inputVals(:);
[sy sx]=size(connectionMatrix);

if (sx~=length(inputVals))
    error('The width of the connection matrix and the length of inputVals should be the same');
end

connectionMatrix=(connectionMatrix~=0); % We do our own normalization
sumNeighbours=sum(connectionMatrix,2); % Although it should be symmetric, we specify row-summation
smoothedVals=connectionMatrix*inputVals;
smoothedVals=smoothedVals./sumNeighbours;
