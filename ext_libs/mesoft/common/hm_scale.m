function [res, errStr]= hm_scale(scale, axis, M)
%
%   function [res, errStr]= hm_scale(scale, axis, M)
%
%                       [sx 0  0  0]
%   M ist 4x4 Matrix    [0  sy 0  0] 
%                       [0  0  sz 0]
%                       [0  0  0  1]
%
%   axis    either string ('x-axis' | 'y-axis' | 'z-axis')
%           or vector with size(axis) = [1 3]
%                       
%   
% Bjoern W. Kreher
% 07/02
%
% UNIX

res= []; errStr= '';

if ~isnumeric(scale) && isvector(scale)
    errStr= 'Error in function [res, errStr]= hm_scale(scale, axis, M): the param alpha have to be numeric and a size of [1 1] or [1 3]';
    return;
end

if ~exist('axis', 'var') && (numel(scale) ~= 3)
    errStr= 'Error in function [res, errStr]= hm_scale(scale, axis, M): the param alpha have to be numeric and a size of [1 1] or [1 3]';
    return;    
end
if exist('axis', 'var') && ~ischar(axis) && ~(isnumeric(axis) && (sum(size(axis) == [1 3]) == 2))
    errStr= 'Error in function [res, errStr]= hm_scale(scale, axis, M): the param alpha have to be numeric and a size of [1 1]';
    return;
end

if exist('M', 'var') && (sum(size(M) == [4 4]) ~= 2)
    errStr= 'Error in function [res, errStr]= hm_scale(scale, axis, M): the matrix M have to be 4x4';
    return;
end

    

%
%       [r11 r12 r13 tx]    [M(1, 1) M(1, 2) M(1, 3) M(1, 4)]
%   M = [r21 r22 r23 ty] =  [M(2, 1) M(2, 2) M(2, 3) M(2, 4)]
%       [r31 r32 r33 tz]    [M(3, 1) M(3, 2) M(3, 3) M(3, 4)]
%       [px  py  pz  s ]    [M(4, 1) M(4, 2) M(4, 3) M(4, 4)]
%
M_tmp= zeros(4, 4);  M_tmp(1:5:16)= 1;
axisStr= {'x-axis'; 'y-axis'; 'z-axis'};

if ~exist('axis', 'var')
    M_tmp([1 6 11])= scale;
elseif ischar(axis)
    if strcmp(axis, axisStr{1})
       M_tmp(1, 1)= scale;
    elseif strcmp(axis, axisStr{2})
       M_tmp(2, 2)= scale;
    elseif strcmp(axis, axisStr{3})
       M_tmp(3, 3)= scale;
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
    rot_invM(1:3, index(2))= rot_invM(1:3, index(2))/norm(rot_invM(1:3, index(2)))
    rot_invM(1:3, index(3))= rot_invM(1:3, index(3))/norm(rot_invM(1:3, index(3)))

    M_tmp= rot_invM * hm_scale(scale, axisStr{index(1)}) * rot_invM';

end

if exist('M')
    res= M_tmp * M;
else
    res= M_tmp;
end
