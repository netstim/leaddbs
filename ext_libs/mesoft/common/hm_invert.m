function [res, errStr]= hm_invert(M)
%
%   function [res, errStr]= hm_invert(M)
%
%                       [r11 r12 r13 tx]
%   M ist 4x4 Matrix    [r21 r22 r23 ty]
%                       [r31 r32 r33 tz]
%                       [px  py  pz  s ]
%                       
%   
% Bjoern W. Kreher
% 07/02
%
% UNIX

res= []; errStr= '';

if sum(size(M) == [4 4]) ~= 2
    errStr= 'Error in function res= hm_invert(M): the matrix M have to be 4x4';
    return;
end

res= inv(M);
