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
        choosedialog(uipatdir,options,handles);
        %answer = questdlg('This dataset appears to have an older lead-dbs data structure. Would you like to Migrate it to BIDS? Please note that support for older version will not provided following this version...','Migrate to BIDS','Yes','No','Read More','Read More');
        %switch answer
        %    case 'Yes'
        %        close(gcf);
        %        lead_migrate(uipatdir,options,handles);
        %    case 'No'
        %        return
        %    case 'Read More'
        %        disp(['Thank you for your interest in Lead-DBS. To better provide support to lead-DBS users and promote data sharing\n' ...
        %               'we have reorganized the structure in which data is stored. This means changes to your input and output data']);
        %end
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
    choosedialog(uipatdir,options,handles);
    
end

function choosedialog(uipatdir,options,handles)

     d = dialog('Position',[500 500 400 150],'Name','Migrate to BIDS');
     %,'Color',[1 1 1]
    
     txt = uicontrol('Parent',d,'Style','text','Position',[45 80 320 50],'String','This dataset appears to have an older lead-dbs data structure. Would you like to Migrate it to BIDS? Please note that support for older version will not provided following this version!');
        
     btn1 = uicontrol('Parent',d,...
            'Position',[50 20 70 25],...
            'String','Yes',...
            'Callback',{@callMigrate,d,uipatdir,options,handles});       
     btn2 = uicontrol('Parent',d,...
         'Position',[170 20 70 25],...
         'String','No',...
         'Callback','close(gcf)');     
     btn3 = uicontrol('Parent',d,...
         'Position', [290 20 70 25],...
         'String','Read More',...
         'Callback',{@callReadMore,d,txt});
  
        function callMigrate(btn1,event,d,uipatdir,options,handles)
           close(d);
           lead_migrate(uipatdir,options,handles);
        function callReadMore(btn3,event,d,txt)
            d.Position = [500 500 400 300];
            txt.Position = [45 230 320 50];
            version_message = ['Thank you for your interest in Lead-DBS! After version 2.5, we have re-organized the way Lead-DBS acesses and stores data.' ,...
                               'This implies changes to the organization of your input and output data. The main objective to set standards for data organization was to promote' ,...
                               'data sharing and open science initiatives. For more information and details on specific changes, please refer to our manual [insert url].' ,...
                               'Lead-migrate is a tool developed to automatically assist you in moving your dataset from the classic lead-dbs to the bidsy-fied version.',...
                               'If you wish to use it, please click on ''Yes''. You can exit this message by clicking on ''No'', but you will not be able to use Lead-DBS.'];
            txt2 = uicontrol('Parent',d,'Style','text','Position',[45 50 320 150],'String',version_message);
            



