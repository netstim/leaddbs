function ea_backuprestore(file)
% Restore file from raw_file if raw_file exists.
% Backup file to raw_file if raw_file doesn't exist.
% Used mainly for coregistration pipeline.

[fpath, fname, ext] = ea_niifileparts(file);

file = [fpath, ext];
backup = [fileparts(fpath), filesep, 'raw_', fname, ext];

if exist(backup, 'file')
    fprintf('\nRestoring %s from raw_%s ...\n\n', [fname, ext], [fname, ext]);
    copyfile(backup, file);
elseif exist(file, 'file')
    fprintf('\nBacking up %s to raw_%s ...\n\n', [fname, ext], [fname, ext]);
    copyfile(file, backup);
else
    fprintf('\nSkipping %s (file does not exist) ...\n\n', [fname, ext]);
end
