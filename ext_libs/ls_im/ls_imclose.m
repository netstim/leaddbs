%% imclose replacement
% Author: Stanislaw Adaszewski, 2012
%%
function R=ls_imclose(A, B)
    R=ls_imerode(ls_imdilate(A, B), B);
end