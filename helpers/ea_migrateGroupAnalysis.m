function ea_migrateGroupAnalysis(source,dest)

LeadGroupFiles = ea_regexpdir(fileparts(source), '(?-i)LEAD_groupanalysis\.mat$', 0, 'f');

if isempty(LeadGroupFiles)
    return;
end

[~, op_filename] = fileparts(dest);
op_filename = regexprep(op_filename, '[\W_]', '');

evalin('base','WARNINGSILENT=1;');
ea_warning('Migrating lead group analysis in the patient folder. If you have different lead group files, please rename them manually.');

try
    load(LeadGroupFiles{1}, 'M');
    groupFolder = fullfile(dest, 'derivatives', 'leadgroup', M.guid);
    ea_mkdir(groupFolder);
    M.root = fullfile(dest, filesep);
    if isfield(M.ui,'groupdir')
        M.ui.groupdir = fullfile(dest, filesep);
    end
    
    bids_name = ['dataset-' op_filename  '_analysis-' M.guid '.mat'];
    save(fullfile(groupFolder, bids_name), 'M');
catch
    evalin('base','WARNINGSILENT=1;');
    ea_warning(['Could not load %s.', LeadGroupFiles]);
end
