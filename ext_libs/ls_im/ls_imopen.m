%% imopen replacement
% Author: Stanislaw Adaszewski, 2012
%%
function R=ls_imopen(A, B)
    R=ls_imdilate(ls_imerode(A, B), B);
end