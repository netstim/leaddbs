% 
% corrects DE grads direction for the obligue slice
% 050613
% ika 
%
% inputs 
%           DeGradDirection         = nominal DE direction 
%           ImageOrientationPatient = orientation of the slice as written dicom header
% output    direction of the gradients in the coordinates of the ImageOrientationPatient-vestor and those normal 
% usage     NewDE_Gradients = rotateDeGradients(DeGradDirection, ImageOrientationPatient)
%
%   test up to now tested for slightly rotated axial slices
%
% history

function NewDE_Gradients = RotateDeGradients(DeGradDirection, ImageOrientationPatient)

DEBUG = 0;

if DEBUG 
    inPlaneVectors =[1     0     0     0     1     0]
    NormalToPlane =[    det([ inPlaneVectors(2), inPlaneVectors(3); inPlaneVectors(5) inPlaneVectors(6)]) ; ...
            det([ inPlaneVectors(3), inPlaneVectors(1); inPlaneVectors(6) inPlaneVectors(4)]) ; ...
            det([ inPlaneVectors(1), inPlaneVectors(2); inPlaneVectors(4) inPlaneVectors(5)] ) ]
    
    %DeGradDirection= [0 -0.8945 0.4472] % ctual unrotated grad dir =ImageComments 

    % ImageOrientationPatient =[0.9589   -0.0340   -0.2817    0.0000    0.9928  -0.1199]  double obligue 
    ImageOrientationPatient = [ 0.9397   -0.0000    0.3420    0.0000    1.0000       0]  % 20 grad 
    ImageOrientationPatient =  [1.0000    0.0000         0   -0.0000    0.9659    0.2588] % 15 grad 
    %ImageOrientationPatient = [ 1  -0.0000    0.    0.0000    1.0000       0]
end

sliceCoor = reshape(ImageOrientationPatient, 3,2);
normalToSlice =   [    det([ sliceCoor(2), sliceCoor(3); sliceCoor(5) sliceCoor(6)]) ; ...
        det([ sliceCoor(3), sliceCoor(1); sliceCoor(6) sliceCoor(4)]) ; ...
        det([ sliceCoor(1), sliceCoor(2); sliceCoor(4) sliceCoor(5)] ) ];
sliceCoor = [sliceCoor'; normalToSlice']';

IntialCoordinate = [1 0 0; 0 1 0; 0 0 1];

for i1=1:3
    for j1=1:3
        RotMatrix(i1,j1)= (sliceCoor(:,i1)'*IntialCoordinate(:,j1) );
    end
end
NewDE_Gradients = RotMatrix *DeGradDirection';

if DEBUG 
    RotMatrix 
    for i1=1:3
        for j1=i1:3
            angle(i1,j1)= 180*acos(sliceCoor(:,i1)'*sliceCoor(:,j1) )/pi
        end
    end
    % strange that I get complex value instad of just 0
    for i1=1:3
        for j1=1:3
            norm(i1,j1)= (sliceCoor(:,i1)'*sliceCoor(:,j1) );
        end
    end
    angle_to_reference =      180*acos(RotMatrix)/pi
end

