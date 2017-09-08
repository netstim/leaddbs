%% applyFSLTransformToPolyCoeffs - transform polynomial coefficient matrix in world coordinats of sourceRefImage to
% world coordinates in targetRefImage applying a FSL (FLIRT) transformation matrix
%
% Parameters: coeffs - Nx3 matrix, coeffs(end,:) is the bias term (x^0)
%            fslTransMat - FSL transformation matrix (4x4)
%            sourceRefImage - Source reference image (NiftiModality or filename)
%            targetRefImage - Target reference image (NiftiModality or filename)
%
% Returns: transCoeffs - Nx3 matrix
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2017
% mail@andreashusch.de, husch.andreas@chl.lu

function transCoeffs = applyFSLTransformToPolyCoeffs(coeffs, fslTransMat, sourceRefImage, targetRefImage)
LtoR = diag([1 1 1]);
transWorldToWorld = flirtmat2worldmatPaCER(fslTransMat,sourceRefImage,targetRefImage);

transfWoT = inv(transWorldToWorld); % take inv here
transfWoT(1:3,4) = 0;
transCoef = [coeffs(1:end-1,:)*LtoR ones(length(coeffs(1:end-1,:)), 1)] * (transfWoT)'; % note the transpose!!
transBias = [coeffs(end,:)*LtoR  1] * inv(transWorldToWorld)'; % note the transpose!!
transCoeffs = [transCoef(:,1:3); transBias(:,1:3)];
end