function ea_tonemapct(options, nativenorm)
if ~exist('nativenorm', 'var') || isempty(nativenorm)
    nativenorm = 'native';
end

switch nativenorm
    case 'native'
        if isfile(options.subj.postopAnat.CT.coreg)
            ea_dispt('Tonemapping coregistered CT...');
            ct = ea_load_nii(options.subj.postopAnat.CT.coreg);
            ct.fname = options.subj.postopAnat.CT.coregTonemap;
            ct.img = tonemap(ct.img);
            ct.dt(1) = 16;
            ea_write_nii(ct);
        else
            fprintf('Coregistered CT image not present. Skipping tonemapping.\n');
        end
    case 'norm'
        if isfile(options.subj.postopAnat.CT.norm)
            ea_dispt('Tonemapping normalized CT...');
            ct = ea_load_nii(options.subj.postopAnat.CT.norm);
            ct.fname = options.subj.postopAnat.CT.normTonemap;
            ct.img = tonemap(ct.img);
            ea_write_nii(ct);
        else
            fprintf('Normalized CT image not present. Skipping tonemapping.\n');
        end
end

ea_dispt('');


function tmct = tonemap(ct)
prefs = ea_prefs;
tmct = ct;
switch prefs.tonemap
    case 'heuristic'
        tmct(ct>0 & ct<80) = tmct(ct>0 & ct<80)/80;
        
        % bone window: center = 300, width = 2000
        tmct(ct>=80 & ct<1300) = (tmct(ct>=80 & ct<1300)-80)/1220;
        
        % saturate above and below levels:
        tmct(ct>=1300) = 1;
        tmct(ct<0)= 0;
    case 'albada'
        tmct(:) = ea_normal(ct(:), 1, 0, ' ', 0, 1, 'TRUE');
        tmct(:) = ea_contrast(tmct(:), 10);
end
