%   function result= reshape_index(index, resolution)
%
%       change a 1-dim Index in a length(resolution)-dim index.
%       index is the 1-dim index resolution describes the size 
%       of the matrix und result is an Vektor of the indices of
%       the matrix
%
%       Example:
%           M= reshape(1:2000, [5, 10, 2, 5, 4]);
%           indexL= find(M == 345);
%           indexM= reshape_index(indexL, size(M));
%           M(indexM(1), indexM(2), indexM(3), indexM(4), indexM(5))
%               345
%
% Bjoern W. Kreher
% 3/02
%
% UNIX

function result= reshape_index(index, resolution)

result= [];%zeros(length(resolution), length(index));
index= reshape(index, [1 prod(size(index))]);
fak= prod(resolution);
if index > fak
    warning('index is out of range');
    result= [];
    return 
end
tmp= index - 1;
for i=length(resolution):-1:1
    fak= fak/resolution(i);
    result(i,:)= floor(tmp/fak) + 1;
    tmp= mod(tmp, fak);
end
result= result';