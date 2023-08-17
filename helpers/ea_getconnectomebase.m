function base=ea_getconnectomebase(cmd,prefs)

if ~exist('prefs','var') || isempty(prefs)
    prefs=ea_prefs('');
end

if isfield(prefs.lc,'datadir')
    base=prefs.lc.datadir;
else
    base=[ea_getearoot,'connectomes',filesep];
end

if exist('cmd','var')
    switch lower(cmd)
        case 'dmri'
            base=[base,'dMRI',filesep];
        case 'fmri'
            base=[base,'fMRI',filesep];
        case 'dmri_multitract'
            base=[base,'dMRI_MultiTract',filesep];
    
    end
end
