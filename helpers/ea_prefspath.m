function prefsPath = ea_prefspath(ext)
% Return prefs path

arguments
    ext {mustBeTextScalar} = '.m'
end

if ~startsWith(ext, '.')
    ext = ['.', ext];
end

prefsPath = fullfile(ea_prefsdir, ['ea_prefs_user', ext]);
