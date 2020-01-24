function ea_export(options, clusterfunctionname)
% This function exports jobs created by the GUI of Lead-DBS. It is
% distributed within Lead-DBS toolbox (www.lead-dbs.org)
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

[fn, pth] = uiputfile({'*.m';'*.json'}, 'Specify location for new job file...','lead_job');
if ~fn % user pressed cancel
    return
end

try
    options = rmfield(options, 'root');
    options = rmfield(options, 'patientname');
end

fID = fopen([pth, fn], 'w');

% export json file if specified
[~, ~, ext] = fileparts(fn);
if strcmp(ext, '.json')
    fwrite(fID, jsonencode(options), 'char'); 
    fclose(fID);
    return
end

% export comments
fprintf(fID, '%s\n', ['function ', fn(1:end-2)]);

fprintf(fID, '%s\n',['% - Lead-DBS Job created on ', datestr(clock), ' -']);
fprintf(fID, '%s\n','% --------------------------------------');
fprintf(fID, '\n');

fprintf(fID, '%s\n',['lead path;']);
fprintf(fID, '\n');

if exist('clusterfunctionname','var') % submit to cluster instead of directly running
    fprintf(fID, '%s\n','% options.uipatdirs = ea_checknoerrorfolders(options.uipatdirs); % uncomment this line to only apply to folders that ran into errors the last time.');
    fprintf(fID, '%s\n','% options.uipatdirs = ea_checknofilefolders(options.uipatdirs, filename); % uncomment this line to only apply to folders that MISS a certain file.');
    fprintf(fID, '%s\n','% options.uipatdirs = ea_checkfilefolders(options.uipatdirs, filename); % uncomment this line to only apply to folders that CONTAIN a certain file.');
    fprintf(fID, '\n');
end

fprintf(fID, '%s\n','options = getoptslocal;');
if ~isempty(options.uipatdirs)
    if exist('clusterfunctionname', 'var') % submit to cluster instead of directly running
        fprintf(fID, '%s\n','patdirs = options.uipatdirs;');
        fprintf(fID, '%s\n', ['clusterfunctionname = ''', clusterfunctionname, ''';']);
        fprintf(fID, '%s\n','for pat = 1:length(patdirs)');
        fprintf(fID, '    %s\n','options.uipatdirs = patdirs{pat};');
        fprintf(fID, '    %s\n','options.root = [fileparts(patdirs{pat}), filesep];');
        fprintf(fID, '    %s\n','[~, options.patientname] = fileparts(patdirs{pat});');
        fprintf(fID, '    %s\n', 'feval(eval([''@'', clusterfunctionname]), options)');
        fprintf(fID, '%s\n', 'end');
    else
        fprintf(fID, '%s\n', 'ea_run(''run'', options);');
    end
else % no patient mode (also e.g. connectome mapper)
    if exist('clusterfunctionname', 'var') % submit to cluster instead of directly running
        fprintf(fID, '%s\n', ['clusterfunctionname = ''', clusterfunctionname, ''';']);
        fprintf(fID, '%s\n', 'feval(eval([''@'', clusterfunctionname]), options)');
    else
        fprintf(fID, '%s\n', 'ea_run(''run'', options);');
    end
end

fprintf(fID, '\n');
fprintf(fID, '\n');

fprintf(fID,'%s\n', 'function options = getoptslocal');

% remove 'modality' field in multiple patients case, will be determined
% automatically.
if length(options.uipatdirs) > 1
    options = rmfield(options, 'modality');
end

% generate code from 'options'
optionsCode = ea_gencode(options, 'options');
for e = 1:length(optionsCode)
    fprintf(fID, '%s\n', optionsCode{e});
end

if ~isdeployed
    edit(fullfile(pth, fn));
end
