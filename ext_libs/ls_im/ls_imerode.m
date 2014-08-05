%% imerode replacement
% Author: Stanislaw Adaszewski, 2012
%%
function R=ls_imerode(A, B)
    A_cls = class(A);
    if ~strcmp(A_cls, 'double')
        A=double(A);
    end
    if ~strcmp(class(B), 'double')
        B=double(B);
    end
    try
        n = feature('numCores');
    catch
        n = 1;
    end
    R = ls_improc(A, B, 1, n);
    if ~strcmp(A_cls, 'double')
        R = eval([A_cls, '(R)']);
    end
end