function [A, index] = ea_sortalike(A, B, mode)
% Sort one vector (A) based on the order defined by another vector (B).

if ~isvector(A) || ~isvector(B)
    error('Inputs should be vector!');
elseif ~strcmp(class(A), class(B))
    error('Inputs should be the same class!')
end

% Additional control flag for cell string. Can either use 'match' or
% 'fuzzy'. With 'match', elements in A should match exactly elements in B
% to be sorted. With 'fuzzy', elements in A starting with elements in B
% would be enough to be sorted.
if ~exist('mode', 'var')
    mode = 'match';
end

if strcmp(mode, 'match') || isnumeric(A)
    [~, index] = ismember(A, B);

    [~, NotFoundIndex] = sort(A(index==0));
    index(index==0) = NotFoundIndex + max(index);

    [~, index] = sort(index);

    A = A(index);
else
    index = zeros(1,length(A));
    for i=1:length(A)
        for j=1:length(B)
            if startsWith(A{i}, B{j})
                index(i) = j;
            end
        end
    end

    for i=1:max(index)
        duplicatedCheck = find(index==i);
        if length(duplicatedCheck)>1
            index(index>i) = index(index>i) + length(duplicatedCheck) - 1;
            [~, idx] = sort(A(index==i));
            for j=1:length(duplicatedCheck)
                index(duplicatedCheck(j)) = idx(j) + i - 1;
            end
        end
    end

    [~, NotFoundIndex] = sort(A(index==0));
    index(index==0) = NotFoundIndex + max(index);

    [~, index] = sort(index);

    A = A(index);
end
