function ea_permtest_readme_append(rundir, title, lines)
% Append a timestamped section to <rundir>/README.md, the human-readable
% log of a permutation run (settings, shuffles, checks, timings, and later
% steps such as statistics and exports).
%
% rundir - run folder.
% title  - section heading.
% lines  - cellstr, one Markdown line each.

fid = fopen(fullfile(rundir, 'README.md'), 'a');
if fid < 0
    ea_error(['Cannot write README.md in ', strrep(strrep(rundir, '\', '\\'), '%', '%%')]);
end
closeFile = onCleanup(@() fclose(fid));
fprintf(fid, '## %s (%s)\n\n', title, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
fprintf(fid, '\n');
