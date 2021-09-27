function [isBIDSRoot,uipatdir] = ea_checkbids(options,handles)

p='/'; % default use root
try
    p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
    load([ea_getearoot,'common',filesep,'ea_recentpatients.mat']);
    p=fileparts(fullrpts{1});
end

uipatdir=ea_uigetdir(p,'Please choose patient folder(s)...');

if isempty(uipatdir)
    return
end

if ~iscell(uipatdir)
    uipatdir = {uipatdir};
end

uipatdir = GetFullPath(uipatdir);

isSubjFolder = 0;
isBIDSRoot = 0;


if length(uipatdir) == 1 % Dragged single folder
    if contains(uipatdir{1}, ['derivatives', filesep, 'leaddbs']) % Is patient folder under derivatives
        isSubjFolder = 1;
    else % Check if it's BIDS root folder
        folders = dir(uipatdir{1});
        folders = {folders.name};
        if ismember('sourcedata', folders) || ismember('rawdata', folders)
            isBIDSRoot = 1;
        end
    end
    if ~isSubjFolder && ~isBIDSRoot
        warndlg('Neither BIDS dataset nor patient folder found, opening in lead migrate!');
        lead_migrate(uipatdir,options,handles)
     
    elseif isBIDSRoot % Is BIDS root folder but it is empty
        rawData = ea_regexpdir([uipatdir{1}, filesep, 'rawdata'], 'sub-', 0);
        rawData = regexprep(rawData, ['\', filesep, '$'], '');
        sourceData = ea_regexpdir([uipatdir{1}, filesep, 'sourcedata'], 'sub-', 0);
        sourceData = regexprep(sourceData, ['\', filesep, '$'], '');
        if isempty(rawData) && isempty(sourceData)
            isBIDSRoot = 0; %change value to zero
            error('Both sourcedata and rawdata folders are empty!');
        end
    end

else % Dragged multiple patient folders
    BIDSRoot = regexp(uipatdir{1}, ['^.*(?=\', filesep, 'derivatives\', filesep, 'leaddbs)'], 'match', 'once');
    patIndx = 1;
    for patdir=1:length(uipatdir)
       if ~contains(uipatdir{1}, ['derivatives', filesep, 'leaddbs'])
          is_not_bids{patIndx} = uipatdir{patdir}; 
       elseif isempty(BIDSRoot)
          error('Please select patient folders under DATASET/derivatives/leaddbs!');
       else
           isBIDSRoot = 1;
       end
       patIndx = patIndx+1;
    end
    warndlg('Neither BIDS dataset nor patient folder found for one or more patients. Opening the patients not in BIDS format inside lead migrate!');
    lead_migrate(uipatdir,options,handles)
    
end


