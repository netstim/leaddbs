function [res, errStr]= hm_rot(alpha, axis, M, flag)
%
%   function [res, errStr]= hm_rot(alpha[rad], axis, M, flag)
%
%                       [r11 r12 r13 0]
%   M ist 4x4 Matrix    [r21 r22 r23 0] = rot(alpha, axis) [* M]
%                       [r31 r32 r33 0]
%                       [ 0   0   0  1]
%
%   axis    either string ('x-axis' | 'y-axis' | 'z-axis')
%           or vector with size(axis) = [1 3]
%                       
%   flag    is set to '3x3', the function will create a 3x3 matrix, othewise a homogen 4x4
% Bjoern W. Kreher
% 07/02
%
% UNIX

res= []; errStr= '';

if ~isnumeric(alpha) && isvector(alpha)
    errStr= 'Error in function [res, errStr]= hm_rot(alpha, axis, M): the param alpha have to be numeric and a zize of [1 1]';
    return;
end

if ~ischar(axis) && ~(isnumeric(axis) && (sum(size(axis) == [1 3]) == 2))
    errStr= 'Error in function [res, errStr]= hm_rot(alpha, axis, M): the param alpha have to be numeric and a zize of [1 1]';
    return;
end


if (nargin >= 3) && ~(isequal(size(M), [4 4]) || isequal(size(M), [3 3]) || isempty(M))
    errStr= 'Error in function [res, errStr]= hm_rot(alpha, axis, M): the matrix M have to be 4x4 or 3x3';
    return;
elseif (nargin < 3)
    M= [];
end

%if ~exist('flag') | ~strcmp(flag, '3x3')
if (nargin < 4) || ~strcmp(flag, '3x3')
    flag= '4x4';
    if isequal(size(M), [3 3])
        tmp= diag([1 1 1 1]);
        tmp(1:3, 1:3)= M;
        M= tmp
    end
end 


%
%       [r11 r12 r13 tx]    [M(1, 1) M(1, 2) M(1, 3) M(1, 4)]
%   M = [r21 r22 r23 ty] =  [M(2, 1) M(2, 2) M(2, 3) M(2, 4)]
%       [r31 r32 r33 tz]    [M(3, 1) M(3, 2) M(3, 3) M(3, 4)]
%       [px  py  pz  s ]    [M(4, 1) M(4, 2) M(4, 3) M(4, 4)]
%
M_tmp= zeros(4, 4);  M_tmp(1:5:16)= 1;
axisStr= {'x-axis'; 'y-axis'; 'z-axis'};

if ischar(axis)
    if strcmp(axis, axisStr{1})
       M_tmp([6 11])= cos(alpha);
       M_tmp(2, 3)= -sin(alpha);  M_tmp(3, 2)= sin(alpha);       
    elseif strcmp(axis, axisStr{2})
       M_tmp([1 11])= cos(alpha);
       M_tmp(1, 3)= sin(alpha);  M_tmp(3, 1)= -sin(alpha);       
    elseif strcmp(axis, axisStr{3})
       M_tmp([1 6])= cos(alpha);
       M_tmp(1, 2)= -sin(alpha);  M_tmp(2, 1)= sin(alpha);       
    else    
       errStr= sprintf('Error in function [res, errStr]= hm_rot(alpha, axis, M): The axis ''%s'' is not supported', axis);    
       return
    end
else
    axis= axis/(norm(axis));
    rot_invM= zeros(4, 4);  rot_invM(1:5:16)= 1;
    
    [dummy, index]= max(axis);
    index= my_mod(index:1:index+3, 3);
    rot_invM(1:3, index(1))= reshape(axis, [3 1]);
    rot_invM(1:3, index(2))= cross(rot_invM(1:3, index(1)), rot_invM(1:3, index(2))); 
    rot_invM(1:3, index(3))= cross(rot_invM(1:3, index(1)), rot_invM(1:3, index(2)));

    rot_invM(1:3, index(2))= rot_invM(1:3, index(2))/norm(rot_invM(1:3, index(2)));
    rot_invM(1:3, index(3))= rot_invM(1:3, index(3))/norm(rot_invM(1:3, index(3)));

    M_tmp= rot_invM * hm_rot(alpha, axisStr{index(1)}) * rot_invM';
end

if ~isempty(M)
    if strcmp(flag, '3x3')
        res= M_tmp(1:3, 1:3) * M;
    else
        res= M_tmp * M;
    end
else
    res= M_tmp;
end

if strcmp(flag, '3x3')
    res= res(1:3, 1:3);
end