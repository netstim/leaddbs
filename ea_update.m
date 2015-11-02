function update=ea_update(varargin)
% lead update

local=ea_getvsn('local');
web=ea_getvsn('web');
if strcmp(local,'Unknown') % full update
    update=1;
elseif strcmp(local,web) % no update
	update=0;
else % incremental update
    update=2;
end

if nargin    
    if strcmp(varargin{1},'force')
    	update=1; 
    else
        return
    end
end

earoot=[fileparts(which('lead')),filesep];

if update
    if update==1 % full update
        patch='lead_dbs';
        updurl=['http://www.lead-dbs.org/release/',patch,'.zip'];
    elseif update==2 % incremental update
        patch=['updates_',strrep(local,'.',''),'-',strrep(web,'.','')];
        updurl=['http://www.lead-dbs.org/release/',patch,'.zip'];
    end

    disp('*** Updating LEAD-DBS. Please do not quit MATLAB.');
    mkdir([earoot,'tmp'])
    disp('Downloading code...');
    urlwrite(updurl,[earoot,'tmp',filesep,'updates.zip']);
    disp('Extracting code...');
    try
        unzip([earoot,'tmp',filesep,'updates.zip'],[earoot,'tmp',filesep]);
    catch
        system(['unzip -q ',earoot,'tmp',filesep,'updates.zip -d ', earoot,'tmp']);
    end
    delete([earoot,'tmp',filesep,'updates.zip']);
    
    if update==2 % delete files during incremental updating
        dels=textscan(fopen([earoot,'tmp',filesep,patch,filesep,'DELETE']),'%s');
        for f=1:length(dels{1})
            if isdir([earoot,dels{1}{f}])
                rmdir([earoot,dels{1}{f}],'s')
            else
                delete([earoot,dels{1}{f}])
            end
        end
        delete([earoot,'tmp',filesep,patch,filesep,'DELETE'])
    end
    
    disp('copying code to the right place...');
    copyfile([earoot,'tmp',filesep,patch,filesep,'*'],earoot,'f');
    copyfile([earoot,'tmp',filesep,patch,filesep,'.version.txt'],earoot,'f');
    disp('Cleaning up...');
    rmdir([earoot,'tmp'],'s')
    disp('Done.');

    disp('Restarting LEAD-DBS.');
    lead;
    disp('*** Update finished.');

else
    disp('LEAD-DBS is already up-to-date!');
    msgbox({'LEAD-DBS aleady up-to-date!';...
            'To get the latest full version, please delete .version.txt and try again later.'},...
            'update','Help')
end






