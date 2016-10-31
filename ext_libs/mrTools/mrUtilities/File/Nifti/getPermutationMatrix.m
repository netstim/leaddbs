% getPermutationMatrix.m
%
%      usage: [permutationMatrix sliceIndex] = getPermutationMatrix(hdr)
%         by: david heeger
%       date: 10/24/07
%    purpose: extracts permutation matrix form nifti hdr
%             sliceIndex contains which index corresponds to
%             Sagittal, Coronal and Axial orientations.
%
function [permutationMatrix sliceIndex] = getPermutationMatrix(hdr)

% check arguments
if ~any(nargin == [1])
  help getPermutationMatrix
  return
end

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane. It is called by both 
% mrLoadRet and mrAlign to decide which way the cardinal axis
% point so that the sagittal/coronal/axial buttons work appropriately

% first get the correct xform
if hdr.sform_code >= 1
  % use sform if sform_code is set to 1
  xform = hdr.sform44(1:3,1:3);
else
  xform = hdr.qform44(1:3,1:3);
end

% now take the QR factorixation
[q,r] = qr(inv(xform));

% and return a permutation matrix
permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);
%permutationMatrix = [q(1,:); abs(q(2,:)); abs(q(3,:))];

% check the euler angles of the permutationMatrix
[phi theta psi] = rot2euler(permutationMatrix);

% some slice orientations that are really coronal
% (or atleast coronal-ish), come out looking axial
% from the permutation matrix. To fix this, we
% look at the theta euler angle and if it is between
% 30 and 45 degrees, we fix the permutaiton matrix
% so that the slices end up coronal.
% This seems to be working on a few test cases, but
% may need tweaking to work in general 11/23/08 -jg.
theta = 180*theta/pi;
%disp(sprintf('(getPermutationMatrix) phi: %0.2f theta: %0.2f psi: %0.2f',180*phi/pi,theta,180*psi/pi));
if (theta > 30) && (theta <= 45)
%  permutationMatrix = [1 0 0;0 0 1;0 1 0]*permutationMatrix;
end

% find out which index corresponds to what oreintation

% Sagittal
[m,sagIndex] = max(permutationMatrix * [1 0 0]');
% Coronal
[m,corIndex] = max(permutationMatrix * [0 1 0]');
% Axial
[m,axialIndex] = max(permutationMatrix * [0 0 1]');

sliceIndex = [sagIndex corIndex axialIndex];

