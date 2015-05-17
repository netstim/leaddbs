

function C = cmult(A,B,idxA,idxB)

szA = size(A);
szB = size(B);

szA = [szA ones(1,max(idxA)-length(szA))];
szB = [szB ones(1,max(idxB)-length(szB))];

if nargin == 2,
    idxA = length(szA);
    idxB = 1;
end;

permA = 1:length(szA);
remainA = setdiff(permA,idxA);
permA = [remainA idxA];

permB = 1:length(szB);
remainB = setdiff(permB,idxB);
permB = [idxB setdiff(permB,idxB)];


C = reshape( ...
    reshape(permute(A,permA),[prod(szA(remainA)) prod(szA(idxA))]) * ...
    reshape(permute(B,permB),[prod(szB(idxB)) prod(szB(remainB))]) ...
    ,[szA(remainA) szB(remainB) 1]);


