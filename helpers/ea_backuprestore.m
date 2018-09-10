function ea_backuprestore(file)
% Backup file to raw_file if raw_file doesn't exist.
% Restore file from raw_file if raw_file exists.
% Used mainly for coregistration pipeline.

[fpath, fname, ext] = ea_niifileparts(file);

file = [fpath, ext];
backup = [fileparts(fpath), filesep, 'raw_', fname, ext];

if ~exist(backup, 'file')
    fprintf('\nBacking up %s to raw_%s ...\n\n', [fname, ext], [fname, ext]);
    copyfile(file, backup);
else
    fprintf('\nRestoring %s from raw_%s ...\n\n', [fname, ext], [fname, ext]);
    copyfile(backup, file);
end
