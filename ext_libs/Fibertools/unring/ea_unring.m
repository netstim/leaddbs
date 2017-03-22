function v = ea_unring(v,params)
% unring - tool for removal of the Gibbs ringing artefact
% Usage: outvol = unring(invol,params)
% Options: invol - input volume 
%          params - 3x1 array with [minW maxW nsh]
%                     nsh discretization of subpixel spaceing (default 20)
%                     minW  left border of window used for TV computation (default 1)
%                     maxW  right border of window used for TV computation (default 3)

disp('Unring DTI...');
if nargin == 1,
    params = [1 3 20];
end;

v = double(v);
v = ringRm(v,params);



%% add methods dump:
cits={
    'Kellner, E., Dhital, B., Kiselev, V. G., & Reisert, M. (2015). Gibbs-ringing artifact removal based on local subvoxel-shifts. Magnetic Resonance in Medicine, n/a?n/a. http://doi.org/10.1002/mrm.26054'
    };
ea_methods(options,['Diffusion-weighted acquisitions were corrected for Gibbs'' ringing artefacts following the approach described in (Kellner 2015; https://bitbucket.org/reisert/unring).'],cits);
    