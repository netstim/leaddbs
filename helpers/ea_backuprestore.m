function ea_backuprestore(file)
% Restore file from preprocessing if preprocessing exists.
% Backup file to preprocessing if preprocessing doesn't exist.
% Used mainly for coregistration pipeline.

space_tag = regexp(file, '(?<=space-).+(?=desc)', 'match', 'once');

backup = strrep(file, [filesep 'coregistration' filesep], [filesep 'preprocessing' filesep]);
backup = strrep(backup, ['space-' space_tag], '');

[~,fname,~] = fileparts(file);
[~,bu_fname,~] = fileparts(backup);

if exist(backup, 'file')
    fprintf('\nRestoring %s from preprocessing %s ...\n\n', fname, bu_fname);
    copyfile(backup, file);
elseif exist(file, 'file')
    fprintf('\nBacking up %s to preprocessing %s ...\n\n', fname, bu_fname);
    copyfile(file, backup);
else
    fprintf('\nSkipping %s (file does not exist) ...\n\n', fname);
end
