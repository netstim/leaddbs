function ea_migrateGroupAnalysis(source,dest)

[root_dir,~,~] = fileparts(source);
[~,op_filename,~] = fileparts(dest); 
op_filename = regexprep(op_filename, '[\W_]', '');
LeadGroupFiles = ea_regexpdir(root_dir, 'lead_groupanalysis.*.mat$',0,'file');
evalin('base','WARNINGSILENT=1;');
ea_warning('Migrating lead group analysis in the patient folder. If you have different lead group files, please rename them manually.');

load(LeadGroupFiles{1})
if exist('M','var')
    lead_path = fullfile(dest,'derivatives','leadgroup',M.guid);
    if ~exist(lead_path,'dir')
        mkdir(lead_path)
    end
    M.root = [dest '/'];
    if isfield(M.ui,'groupdir')
        M.ui.groupdir = [dest '/'];
    end
    
    bids_name = ['dataset-' op_filename  '_analysis-' M.guid '.mat'];
    save(fullfile(lead_path,bids_name),'M');
else
    evalin('base','WARNINGSILENT=1;');
    ea_warning(['Could not load %s.',LeadGroupFiles]);
   
end
