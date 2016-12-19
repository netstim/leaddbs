function [res, errStr]= hm_trans(t_vect, M)
%
%   function [res, errStr]= hm_trans(t_vect, M)
%
%                       [1 0 0 tx]
%   M ist 4x4 Matrix    [0 1 0 ty] [* M]= t_vect [* M]
%                       [0 0 1 tz]
%                       [0 0 0  1]
%%   
% Bjoern W. Kreher
% 07/02
%
% UNIX

res= []; errStr= '';

if ~(isnumeric(t_vect) && (numel(t_vect, 2) == 3))
    errStr= 'Error in function [res, errStr]= hm_trans(t_vect, M): the param t_vect have to be numeric and size of [1 3]';
    return;
end

if ~exist('M') ||isempty(M)
    M= diag(ones(1, 4));
elseif (sum(size(M) == [4 4]) ~= 2)
    errStr= 'Error in function [res, errStr]= hm_trans(t_vect, M): the matrix M have to be 4x4';
    return;
end

%
%       [r11 r12 r13 tx]    [M(1, 1) M(1, 2) M(1, 3) M(1, 4)]
%   M = [r21 r22 r23 ty] =  [M(2, 1) M(2, 2) M(2, 3) M(2, 4)]
%       [r31 r32 r33 tz]    [M(3, 1) M(3, 2) M(3, 3) M(3, 4)]
%       [px  py  pz  s ]    [M(4, 1) M(4, 2) M(4, 3) M(4, 4)]
%
% M_tmp= zeros(4, 4);  M_tmp(1:5:16)= 1;
% 
% M_tmp(1:3, 4)= reshape(t_vect, [3 1]);
% 
% if exist('M')
%     res= M_tmp * M;
% else
%     res= M_tmp;
% end


tmp= zeros([size(t_vect, 1) 4]);
for i= 1:size(t_vect, 1)
    tmp(i, :)= M*[t_vect(i, :), 1];    
end
