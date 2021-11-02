function [A, index] = ea_sortalike(A, B)
% Sort one vector (A) based on the order defined by another vector (B).

if ~isvector(A) || ~isvector(B)
    error('Inputs should be vector!');
elseif ~strcmp(class(A), class(B))
    error('Inputs should be the same class!')
end

[~, index] = ismember(A, B);

[~, NotFoundIndex] = sort(A(index==0));
index(index==0) = NotFoundIndex + max(index);

[~, index] = sort(index);

A = A(index);
