function ea_methods(options, parsestr, refs)
% Show method text and dump it to patient directory if possible
%
% options can either be a struct or a string with the patient directory.

h = dbstack;
try
    callingfunction = h(2).name;
catch
    callingfunction = 'base';
end

% escape '\' in file paths for Windows
if ispc
    filePath = regexp(parsestr, '([a-zA-Z]:\\|\\\\).*?\.[a-z]+(?=,?\s+)', 'match');
    escapedFilePath = ea_path_escape(filePath);
    for i=1:numel(filePath)
        parsestr = strrep(parsestr, filePath{i}, escapedFilePath{i});
    end
end

expstr = '\n\n';
expstr = [expstr, char(datetime('now')), ': ', callingfunction, '\n', '--------------------------\n', parsestr];

if exist('refs','var') % add refs
    if isrow(refs)
        refs = refs';
    end
    expstr = [expstr, '\n\nReferences:\n', '--------------------------\n'];
    expstr = [expstr, strjoin(strcat(num2str((1:numel(refs))'), {') '}, refs), '\n')];
end

expstr=[expstr, '\n***\n\n'];

prefs = ea_prefs;
if prefs.machine.methods_show
    try
        ea_methodsdisp({expstr});
    end
else
    fprintf(expstr);
end

if ~isempty(options)
    % 'options' is a directory in BIDS dataset
    if ischar(options) && contains(options, ['derivatives', filesep, 'leaddbs'])
        options = ea_getptopts(options);
    end

    if isstruct(options) && isfield(options, 'subj') && isfield(options.subj, 'methodLog')
        % Get methodLog path in BIDS subj folder
        ea_mkdir(fileparts(options.subj.methodLog));
        methodfile = fopen(options.subj.methodLog, 'a');
    elseif ischar(options) && isfolder(options)
        % options is a directory outside of BIDS dataset
        methodfile = fopen(fullfile(options, 'methods.txt'), 'a');
    end

    if exist('methodfile', 'var')
        fprintf(methodfile, expstr);
        fclose(methodfile);
    end
end
