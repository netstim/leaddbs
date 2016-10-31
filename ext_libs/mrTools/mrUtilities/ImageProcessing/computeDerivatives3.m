function [fx,fy,fz,ft] = computeDerivatives3(in1,in2,shape)
% 
% function [fx,fy,fz,ft] = computeDerivatives3(in1,in2)
%
% in1 and in2 are volumes, 3d arrays
%
% [fx,fy,fz,ft] are volumes, derivatives of the volumes

if ~exist('shape','var')
    shape = 'valid';
end

filter = [0.03504 0.24878 0.43234 0.24878 0.03504];
dfilter = [0.10689 0.28461 0.0  -0.28461  -0.10689];

dz1 = convXYsep(convn(in1, permute(dfilter,[1 3 2]),shape),filter,filter,shape);
tmp1 = convn(in1, permute(filter,[1 3 2]),shape); %this does the same thing as tmp1 = convZ(in1,filter) if shape=='valid', but much faster, and you can choose the shape...
dx1 = convXYsep(tmp1,dfilter,filter,shape);
dy1 = convXYsep(tmp1,filter,dfilter,shape);
blur1 = convXYsep(tmp1,filter,filter,shape);

dz2 = convXYsep(convn(in2, permute(dfilter,[1 3 2]),shape),filter,filter,shape);
tmp2 = convn(in2, permute(filter,[1 3 2]),shape);
dx2 = convXYsep(tmp2,dfilter,filter,shape);
dy2 = convXYsep(tmp2,filter,dfilter,shape);
blur2 = convXYsep(tmp2,filter,filter,shape);

fx=(dx1+dx2)/2;
fy=(dy1+dy2)/2;
fz=(dz1+dz2)/2;
ft=(blur2-blur1);

