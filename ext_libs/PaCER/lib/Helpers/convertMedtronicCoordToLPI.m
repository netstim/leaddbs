%% convertMedtronicCoordToLPI - be careful with that
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dept. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function pointsLPI = convertMedtronicCoordToLPI(pointsInMm, refNii)
    refNii.load();
    % convert system of Medtronic to LPI [mm] in
    % reference image space (origin LPI image corner)
    origin = refNii.voxdim .* refNii.voxsize; %RAS corner in [mm]
    origin(1) = 0; %LAS corner
    
    pointsLPI = abs(origin - pointsInMm);
end