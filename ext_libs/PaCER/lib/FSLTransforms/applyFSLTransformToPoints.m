%% applyFSLTransformToPoints - transform points in world coordinats of sourceRefImage to
% world coordinates in targetRefImage applying a FSL (FLIRT) transformation
% matrix. world coordinates in PaCER style. (note L/R orientation in plots
% i.e. axis ij vs axis xy!)
%
% Parameters: points - Nx3 matrix 
%            fslTransMat - FSL transformation matrix (4x4)
%            sourceRefImage - Source reference image (NiftiModality or filename)
%            targetRefImage - Target reference image (NiftiModality or filename)
% 
% Returns: transPoints - Nx3 matrix 
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function transPoints = applyFSLTransformToPoints(points, fslTransMat, sourceRefImage, targetRefImage)
LtoR = diag([1 1 1]);
pointsR = points * LtoR;

transWorldToWorld = flirtmat2worldmatPaCER(fslTransMat,sourceRefImage,targetRefImage);
transPoints = (inv(transWorldToWorld) * [pointsR ones(length(pointsR), 1)]')' ; %#ok<MINV> %note the inv()!%
end