function [res, errStr]= hm_mult(M, vec)
%
%   function res= hm_mult(M, vec)
%
%                       [r11 r12 r13 tx]
%   M ist 4x4 Matrix    [r21 r22 r23 ty]  * [vec(i,:) 0]' = erg(size(vec))
%                       [r31 r32 r33 tz]
%                       [px  py  pz  s ]
%                       
%   
% Bjoern W. Kreher
% 07/02
%
% UNIX


res= []; errStr= '';


if size(vec, 2) ~= 3
    errStr= 'Error in function res= hm_mult(M, vec): the matrix vector vec have to be the size [n, 3]';
    return;
end
    
res= zeros(size(vec));

if sum(size(M) == [4 4]) == 2
    tmpComp= vec(:, 1)*M(4, 1) + vec(:, 2)*M(4, 2) + vec(:, 3)*M(4, 3) + M(4, 4);

    for i= 1:3
        res(:, i)= (vec(:, 1)*M(i, 1) + vec(:, 2)*M(i, 2) + vec(:, 3)*M(i, 3) + M(i, 4))./tmpComp;
    end
elseif sum(size(M) == [3 3]) == 2
    for i= 1:3
        res(:, i)= vec(:, 1)*M(i, 1) + vec(:, 2)*M(i, 2) + vec(:, 3)*M(i, 3);
    end
else
    errStr= 'Error in function res= hm_mult(M, vec): the matrix M have to be 4x4 or 3x3';
    return;
end

