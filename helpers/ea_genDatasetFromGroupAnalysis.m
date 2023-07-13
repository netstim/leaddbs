function [groupdir, analysisFile] = ea_genDatasetFromGroupAnalysis(analysisFile)
% Function to create dataset folder from group analysis file (single file
% outside of a dataset)

analysisFile = GetFullPath(analysisFile);

dataset = regexp(analysisFile, '(?<=dataset-)(.+)(?=_analysis-.+\.mat$)', 'match', 'once');
analysis = regexp(analysisFile, '(?<=dataset-.+_analysis-)(.+)(?=\.mat$)', 'match', 'once');
if isfolder(fullfile(fileparts(analysisFile), dataset))
    for i=1:20
        if ~isfolder(fullfile(fileparts(analysisFile), [dataset, num2str(i)]))
            dataset = inputdlg(sprintf('Folder ''%s'' already exists.\nPlease specify a new dataset name:', dataset), 'New Dataset Name', [1 35], {[dataset, num2str(i)]});
            break;
        end
    end
    if isempty(dataset)
        error('Please specify a new dataset name!');
    else
        dataset = dataset{1};
        movefile(analysisFile, regexprep(analysisFile, '(?<=dataset-)(.+)(?=_analysis-.+\.mat$)', dataset));
        analysisFile = regexprep(analysisFile, '(?<=dataset-)(.+)(?=_analysis-.+\.mat$)', dataset);
    end
end

ea_cprintf('CmdWinWarnings', 'Creating new dataset folder: %s\n', fullfile(fileparts(analysisFile), dataset));
groupdir = fullfile(fileparts(analysisFile), dataset, 'derivatives', 'leadgroup', analysis, filesep);
leaddbsFolder = fullfile(fileparts(analysisFile), dataset, 'derivatives', 'leaddbs');
datasetFolder = fullfile(fileparts(analysisFile), dataset);

ea_mkdir(groupdir);
ea_mkdir(leaddbsFolder);

load(analysisFile, 'M');
M.root = groupdir;

miniset.name = dataset;
miniset.numSubj = length(M.patient.list);
miniset.recon = 0;
miniset.vta = 0;
miniset.stats = 0;
miniset.groupAnalysis = analysis;
miniset.timeStamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

for p=1:length(M.patient.list)
    [oldPatientFolder, patientTag] = fileparts(M.patient.list{p});
    M.patient.list{p} = strrep(M.patient.list{p}, oldPatientFolder, leaddbsFolder);
    if isfield(M, 'stats')
        miniset.stats = 1;
        ea_mkdir(fullfile(leaddbsFolder, patientTag));
        ea_stats = M.stats(p).ea_stats;
        save(fullfile(leaddbsFolder, patientTag, [patientTag, '_desc-stats.mat']), 'ea_stats');
    end

    if isfield(M, 'elstruct')
        miniset.recon = 1;
        ea_mkdir(fullfile(leaddbsFolder, patientTag, 'reconstruction'));
        for e=1:length(M.elstruct(p).coords_mm)
            reco.props(e).elmodel = M.elstruct(p).elmodel;
            reco.props(e).manually_corrected = 1;
            reco.mni.coords_mm(e) = M.elstruct(p).coords_mm(e);
            reco.mni.markers(e) = M.elstruct(p).markers(e);
            reco.mni.trajectory(e) = M.elstruct(p).trajectory(e);
            reco.electrode(e).dbs.elmodel = M.elstruct(p).elmodel;
        end
        save(fullfile(leaddbsFolder, patientTag, 'reconstruction', [patientTag, '_desc-reconstruction.mat']), 'reco');
    end
end

movefile(analysisFile, groupdir);
analysisFile = ea_getGroupAnalysisFile(groupdir);
save(analysisFile, 'M');

savejson('', miniset, fullfile(datasetFolder, 'miniset.json'));
