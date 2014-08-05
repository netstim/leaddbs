%% imdilate replacement
% Author: Stanislaw Adaszewski, 2012
%%
function R=ls_imdilate(A, B)
    A_cls = class(A);
    if ~strcmp(A_cls, 'double')
        A=double(A);
    end
    s = size(B);
    if (size(s,2) < 3)
        s(3)=1;
    end
    s2 = zeros(1,3);
    for i=1:size(s,2)
        s2(i)=s(i)+1-mod(s(i),2);
    end
    z = zeros(s2);
    z((2-mod(s(1),2)):end,(2-mod(s(2),2)):end,(2-mod(s(3),2)):end) = B;
    B = z;
    if ~strcmp(class(B), 'double')
        B=double(B);
    end
    try
        n = feature('numCores');
    catch
        n = 1;
    end
    R = ls_improc(A, B, 0, n);
    if ~strcmp(A_cls, 'double')
        R = eval([A_cls, '(R)']);
    end
end