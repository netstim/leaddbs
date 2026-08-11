function ea_unified_perm_readme_append(outdir, sectionTitle, body)
% Appends a titled, timestamped section to README.txt in a permutation-
% testing run folder (outdir). Every function that writes outputs into a
% run's folder (ea_unified_perm_sweetspot, ea_unified_perm_sweetspot_stats,
% ea_unified_perm_sweetspot_predict) calls this once, describing what it
% just added -- so the folder is self-documenting even outside MATLAB
% (e.g. browsing it in Finder months later, or handing it to a
% collaborator) without needing to find/read this codebase's source.
%
% sectionTitle - short header for this section, e.g. 'Voxelwise permutation test'.
% body         - the section's text (a plain string; embedded newlines are fine).

readmeFile = fullfile(outdir, 'README.txt');
fid = fopen(readmeFile, 'a');
if fid == -1
    warning('ea_unified_perm_readme_append:cannotWrite', 'Could not open %s for writing -- README not updated.', readmeFile);
    return
end
fprintf(fid, '%s\n', repmat('=', 1, numel(sectionTitle)));
fprintf(fid, '%s\n', sectionTitle);
fprintf(fid, '%s\n', repmat('=', 1, numel(sectionTitle)));
fprintf(fid, 'Written: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '%s\n\n', body);
fclose(fid);
end
