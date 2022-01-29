function label = ea_getConnLabel(connectome)
% Return a label for specific connectome

connBaseFolder = ea_getconnectomebase;
connFolder = fullfile(connBaseFolder, {'dMRI', 'fMRI', 'dmri_multitract'}', connectome);
connFolder = connFolder(isfolder(connFolder));

if isempty(connFolder)
    error('Specified connectome "%s" not found!', connectome);
elseif isfile(fullfile(connFolder{1}, 'dataset_info.json'))
    dataset = loadjson(fullfile(connFolder{1}, 'dataset_info.json'));
elseif isfile(fullfile(connFolder{1}, 'dataset_info.mat'))
    load(fullfile(connFolder{1}, 'dataset_info.mat'), 'dataset');
else
    error('No dataset_info found for the specified connectome "%s"!', connectome);
end

label = dataset.tag;
