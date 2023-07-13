function ea_migrateGroupAnalysis(source,dest,id_list)

% source :  full path of the specific lead group file to import (safest way 
%           in case several lead groups exist)
% dest :    BIDS root directory  
% id_list : optional - a N x 2 cell with the correspondance between old classic 
%           IDs and new bids IDs, in case patients were renamed

[~, op_filename] = fileparts(dest);
op_filename = regexprep(op_filename, '[\W_]', '');

try
    load(source, 'M');
    groupFolder = fullfile(dest, 'derivatives', 'leadgroup', M.guid);
    ea_mkdir(groupFolder);
    M.root = fullfile(dest, filesep);
    
    bids_name = ['dataset-' op_filename  '_analysis-' M.guid '.mat'];

    if ~isfield(M.ui, 'mirrorsides')
        M.ui.mirrorsides=0; 
    end 

    for ii = 1:length(M.patient.list)
        
        [~, old_pname] = fileparts(M.patient.list{ii}); 
 
        new_path = fullfile(dest, 'derivatives', 'leaddbs'); 

        if exist('id_list', 'var') % get ID from the list
            if sum(strcmp(old_pname, id_list(:,1))) == 1
                new_pname = id_list{strcmp(old_pname, id_list(:,1)), 2}; 
            else 
                ea_warning(['Lead group import failed - Could not rename ' old_pname '!'])
                return 
            end
        else % automatic renaming
            is_underscore = regexp(old_pname,'.*_|-\s.*','match');
            if ~isempty(is_underscore)
                new_pname = [ 'sub-' regexprep(old_pname, '[\W_]', '') ];
            else 
                new_pname = ['sub-' old_pname]; 
            end 
        end

        if ~isdir(fullfile(new_path, new_pname))
            ea_warning(['Could not complete the lead group import. Subject directory ' fullfile(new_path, new_pname) ' not found!'])
            return; 
        
        else % rename
            M.patient.list{ii} = fullfile(new_path, new_pname); 
            M.elstruct(ii).name = new_pname; 
            if isfield(M, 'stats') && ~isempty(M.stats)
                M.stats(ii).ea_stats.patname = { new_pname, new_pname};
                M.stats(ii).ea_stats.electrodes.name = new_pname; 
            end
        end

    end 

    save(fullfile(groupFolder, bids_name), 'M');

catch
    evalin('base','WARNINGSILENT=1;');
    ea_warning(['Could not complete the lead group import']);
end
