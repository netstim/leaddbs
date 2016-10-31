% rot2euler.m
%
%        $Id$ 
%      usage: rot2euler(rotationMatrix)
%         by: justin gardner
%       date: 11/23/08
%    purpose: get euler angles in radians from rotation matrix
%             see: http://en.wikipedia.org/wiki/Rotation_representation_(mathematics)#Rotation_matrix_.28or_direction_cosine_matrix.29
function [phi theta psi] = rot2euler(rotationMatrix)

% check arguments
if ~any(nargin == [1])
  help rot2euler
  return
end

phi = cart2pol(rotationMatrix(3,1),rotationMatrix(3,2));
theta = acos(rotationMatrix(3,3));
psi = cart2pol(rotationMatrix(1,3),rotationMatrix(2,3));
