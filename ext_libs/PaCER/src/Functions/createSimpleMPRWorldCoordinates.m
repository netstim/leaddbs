%% Convenience Wrapper for the MPRWorldCoordinate Class by Florian Bernard
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dept. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function [mpr, oc] = createSimpleMPRWorldCoordinates(nii)
%mprPanel = uipanel('Parent', gcf);
axesObj = gca;%axesObj = axes('Parent', mprPanel);
mpr = MPRWorldCoordinates(axesObj, nii.img, [0 0 0], nii);
view(160,10);
%oc = OrientationCube(mprPanel, axesObj);
%axes(axesObj);
end