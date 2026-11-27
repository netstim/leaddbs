function files = ea_permtest_export_fibers(stats, mapset, outdir)
% Export the significant fibers of a fiber permutation run for 3D
% visualisation: one file per analysis, family, direction and method
% (single-step / step-down), written only when at least one fiber
% survives. The files use the explorer's fiber export format -- fibcell
% {right, left} and vals {right, left} -- so they work with Lead-DBS
% tools that read that format (e.g. ea_discfibers2trk for TrackVis /
% Surf Ice, ea_discfibers2raw for patient space).
%
% File: <outdir>/<analysis>_<pos|neg>_<family>_<singlestep|stepdown>.mat
%   fibcell  {1 x Nsides} coordinates of the significant fibers per side
%   vals     {1 x Nsides} their real R
%   pFWE     {1 x Nsides} their FWE-corrected p-values
%   fibIdx   {1 x Nsides} their rows in the explorer's fiber list
%   info     analysis, family, direction, method, alpha, run and stats folders
%
% stats  - from ea_permtest_stats (on a fiber map set).
% mapset - from ea_permtest_read_fiberrun.
% outdir - export folder (created).
%
% files  - cellstr of written files.

if ~strcmp(mapset.modality, 'fibers')
    ea_error('ea_permtest_export_fibers expects a fiber map set.');
end
[ok, msg] = mkdir(outdir);
if ~ok
    ea_error(['Cannot create the export folder:\n', strrep(strrep(outdir, '\', '\\'), '%', '%%'), '\n', msg]);
end
F = load(mapset.meta.fibersFile, 'fibersets');
fibersets = F.fibersets;
Nsides = max([mapset.maps.side]);

files = {};
methods = {'singlestep', 'stepdown'};
for a = 1:numel(mapset.analyses)
    fs = fibersets(mapset.analyses(a).fiberset);
    for tail = unique({stats.maxstat.tail}, 'stable')
        for family = unique({stats.maxstat.family}, 'stable')
            % entries of this family/direction that contain maps of analysis a
            E = stats.maxstat(strcmp({stats.maxstat.tail}, tail{1}) & strcmp({stats.maxstat.family}, family{1}));
            for meth = methods
                fibcell = cell(1, Nsides);
                vals = cell(1, Nsides);
                pFWE = cell(1, Nsides);
                fibIdx = cell(1, Nsides);
                n = 0;
                for e = E
                    for j = 1:numel(e.members)
                        i = e.members(j);
                        if mapset.maps(i).analysis ~= a
                            continue
                        end
                        side = mapset.maps(i).side;
                        if strcmp(meth{1}, 'singlestep')
                            sig = e.sigSingle{j};
                            p = e.pSingle{j};
                        else
                            sig = e.sigStepdown{j};
                            p = e.pStepdown{j};
                        end
                        fibcell{side} = fs.fibcell{side}(sig);
                        vals{side} = mapset.maps(i).Remp(sig);
                        pFWE{side} = p(sig);
                        fibIdx{side} = mapset.maps(i).unitIdx(sig);
                        n = n + nnz(sig);
                    end
                end
                if n == 0
                    continue
                end
                info = struct('analysis', mapset.analyses(a).label, 'family', family{1}, 'direction', tail{1}, ...
                    'method', meth{1}, 'alpha', stats.settings.alphaPerTail, 'Nfibers', n, ...
                    'run', mapset.source, 'stats', stats.settings.outdir, 'created', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<NASGU>
                file = fullfile(outdir, sprintf('%s_%s_%s_%s.mat', regexprep(mapset.analyses(a).label, '[^\w-]', '_'), ...
                    tail{1}, family{1}, meth{1}));
                save(file, 'fibcell', 'vals', 'pFWE', 'fibIdx', 'info');
                files{end+1} = file; %#ok<AGROW>
            end
        end
    end
end
fprintf('%d fiber export file(s) written to %s\n', numel(files), outdir);
