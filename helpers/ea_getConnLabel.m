function label = ea_getConnLabel(connectome, subset)
% Return a label for specific connectome

% fMRI connectome with subset in the 'connectome' parameter
if contains(connectome, '>')
    subset = regexprep(connectome, '.*> *', ''); % override 'subset' parameter
    connectome = regexprep(connectome, ' *>.*', '');
end

connBaseFolder = ea_getconnectomebase;
connFolder = fullfile(connBaseFolder, {'dMRI', 'fMRI', 'dmri_multitract'}', connectome);
connFolder = connFolder(isfolder(connFolder));

if isempty(connFolder)
    error('Specified connectome "%s" not found!', connectome);
elseif isfile(fullfile(connFolder{1}, 'dataset_info.json'))
    dataset = loadjson(fullfile(connFolder{1}, 'dataset_info.json'));
else
    error('No dataset_info found for the specified connectome "%s"!', connectome);
end

label = dataset.tag;

if exist('subset', 'var')
    label = [label, 'x', regexprep(subset, '[\W_]', '')];
end
