function ea_update_data(varargin)
% update lead data

updurl = 'http://www.lead-dbs.org/release/download.php';

if nargin
    if strcmp(varargin{1},'full')
    	id = 'data';
    else
        id = 'updates_data';
    end
else
    id = 'updates_data';
end

earoot = ea_getearoot;

if strcmp(id, 'updates_data')
    status = getResponseCode(openConnection(java.net.URL([updurl,'?id=',id])));
    if status ~= 200
        info=sprintf(['No data update found!\n',...
                      'Alternatively, you can download the full data zip from:\n', ...
                      'http://www.lead-dbs.org/release/download.php?id=data']);
        disp(info);
        msgbox(info,'No Update','Help');
        return
    end
    disp('Downloading lead data update...');
elseif strcmp(id, 'data')
    disp('Downloading lead data...');
end

try
    webopts=weboptions('Timeout',5);
    websave([earoot,'tmp',filesep,'updates.zip'],updurl,'id',id,webopts);
catch
    try
        urlwrite([updurl,'?id=',id],[earoot,'tmp',filesep,'updates.zip'],'Timeout',5);
    catch
        info='Download error! Please retry later.';
        disp(info);
        msgbox(info,'Error','Error')
        return
    end
end

disp('Extracting code...');
try
    unzip([earoot,'tmp',filesep,'updates.zip'],[earoot,'tmp',filesep,id]);
catch
    system(['unzip -q ',earoot,'tmp',filesep,'updates.zip -d ', earoot,'tmp', filesep, id]);
end
delete([earoot,'tmp',filesep,'updates.zip']);

% delete unused data
% if strcmp(id, 'updates_data')
%     try
%         if exist([earoot,'tmp',filesep,id,filesep,'DELETE'], 'file')
%             disp('Deleting unused data...');
%             dfid = fopen([earoot,'tmp',filesep,id,filesep,'DELETE']);
%             dels = textscan(dfid,'%s');
%             fclose(dfid);
%             for f=1:length(dels{1})
%                 if isdir([earoot,dels{1}{f}])
%                     rmdir([earoot,dels{1}{f}],'s')
%                 else
%                     delete([earoot,dels{1}{f}])
%                 end
%             end
%             delete([earoot,'tmp',filesep,id,filesep,'DELETE'])
%         end
%     catch
%         disp('Error while deleting some data files. You may ignore this.');
%     end
% end

disp('Copying lead data...');
copyfile([earoot,'tmp',filesep,id,filesep,'*'],earoot,'f');

disp('Cleaning up...');
rmdir([earoot,'tmp'],'s')
disp('Done.');

disp('*** Successfully updated lead data. ***');
