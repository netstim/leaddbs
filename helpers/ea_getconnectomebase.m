function base=ea_getconnectomebase(cmd)

prefs=ea_prefs('');
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
    end
end
