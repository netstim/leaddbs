%% getIntensityForLabel
%
% Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
function intensity = getIntensityForLabel(labelNames, label)
    intensity = find(strcmp(labelNames, label)) - 1;
end