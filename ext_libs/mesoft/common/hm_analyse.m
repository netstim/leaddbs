function [P, errStr]= hm_analyse(M)
%
%   function [res, errStr]= hm_analyse(inAy)
%
%
%   res= [a b c x y z]
%
%   CODE WAS TAKEN FROM SPM5 (SPM_IMATRIX)
%
% Bjoern W. Kreher
% 09/02
%
% UNIX
%
% returns the parameters for creating an affine transformation
% FORMAT P = spm_imatrix(M)
% M      - Affine transformation matrix
% P      - Parameters (see spm_matrix for definitions)
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience


% res= zeros(1,6);
% res(3)= atan2((inAy(2, 1)), inAy(1, 1));
% res(2)= atan2(-inAy(3, 1), cos(res(3))*inAy(1, 1) + sin(res(3))*inAy(2, 1));
% res(1)= atan2(inAy(3, 2), inAy(3, 3));
%
% if prod(size(inAy)) == 16
%     res(4:6)= inAy(1:3, 4)';
% end

if nargin == 0
    errStr = 'Error in hm_analyse: no input matrix given';
else
    errStr = [];
end

if ~isequal(size(M), [4,4])
    errStr = 'Error in hm_analyse: matrix has wrong dimensions, has to be a 4 x 4 matrix!';
else
    errStr = [];
end

% Translations and zooms
%-----------------------------------------------------------------------
R         = M(1:3,1:3);
C         = chol(R'*R);
P         = [M(1:3,4)' 0 0 0  diag(C)'  0 0 0];
if det(R)<0, P(7)=-P(7);end % Fix for -ve determinants

% Shears
%-----------------------------------------------------------------------
C         = diag(diag(C))\C;
P(10:12)  = C([4 7 8]);
R0        = hm_create([0 0 0  0 0 0 P(7:12)]);
R0        = R0(1:3,1:3);
R1        = R/R0;

% This just leaves rotations in matrix R1
%-----------------------------------------------------------------------
%[          c5*c6,           c5*s6, s5]
%[-s4*s5*c6-c4*s6, -s4*s5*s6+c4*c6, s4*c5]
%[-c4*s5*c6+s4*s6, -c4*s5*s6-s4*c6, c4*c5]

P(5) = asin(rang(R1(1,3)));
if (abs(P(5))-pi/2)^2 < 1e-9,
    P(4) = 0;
    P(6) = atan2(-rang(R1(2,1)), rang(-R1(3,1)/R1(1,3)));
else
    c    = cos(P(5));
    P(4) = atan2(rang(R1(2,3)/c), rang(R1(3,3)/c));
    P(6) = atan2(rang(R1(1,2)/c), rang(R1(1,1)/c));
end
end
%return;

% There may be slight rounding errors making b>1 or b<-1.
function a = rang(b)
a = min(max(b, -1), 1);
%return;
end
