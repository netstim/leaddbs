function s = size3(A)

s = size(A);
if length(s)<3
    s = [s ones(1,3-length(s))];
end
