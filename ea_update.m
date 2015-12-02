function ea_update(varargin)
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

updurl = 'http://www.lead-dbs.org/release/download.php';
if update
    try
        if update==1 % full update
            id='lead_dbs';
        elseif update==2 % incremental update
            id=['updates_',strrep(local,'.',''),'-',strrep(web,'.','')];
        end

        disp('*** Updating LEAD-DBS. Please do not quit MATLAB.');
        mkdir([earoot,'tmp'])
        disp('Downloading code...');
        try
            webopts=weboptions('Timeout',5);
            websave([earoot,'tmp',filesep,'updates.zip'],updurl,'id',id,webopts);
        catch
            try
                urlwrite([updurl,'?id=',id],[earoot,'tmp',filesep,'updates.zip'],'Timeout',5);
            catch
                info='Download error! Please retry later.';
                disp(info);
                msgbox(info,'Update','Error')
                return
            end
        end

        disp('Extracting code...');
        try
            unzip([earoot,'tmp',filesep,'updates.zip'],[earoot,'tmp',filesep]);
        catch
            system(['unzip -q ',earoot,'tmp',filesep,'updates.zip -d ', earoot,'tmp']);
        end
        delete([earoot,'tmp',filesep,'updates.zip']);

        disp('Deleting outdated code...');
        try
            if update==2 % delete files during incremental updating
                dels=textscan(fopen([earoot,'tmp',filesep,id,filesep,'DELETE']),'%s');
                for f=1:length(dels{1})
                    if isdir([earoot,dels{1}{f}])
                        rmdir([earoot,dels{1}{f}],'s')
                    else
                        delete([earoot,dels{1}{f}])
                    end
                end
                delete([earoot,'tmp',filesep,id,filesep,'DELETE'])
            end
        catch
            disp('Error while deleting some files. You may ignore this.');
        end

        disp('Copying new code...');
        copyfile([earoot,'tmp',filesep,id,filesep,'*'],earoot,'f');
        copyfile([earoot,'tmp',filesep,id,filesep,'.version.txt'],earoot,'f');
        disp('Cleaning up...');
        rmdir([earoot,'tmp'],'s')
        disp('Done.');

        disp('Restarting LEAD-DBS.');
        lead;
        disp('*** Update finished.');
    catch
        info=sprintf(['Failed to update!\n',...
                      'Alternatively, you can download the latest verion from: http://www.lead-dbs.org/release/download.php']);
        disp(info);
        msgbox(info,'Update','Error');
    end
else
    info=sprintf(['LEAD-DBS aleady up-to-date!\n',...
                  'Alternatively, you can re-download the latest verion from: http://www.lead-dbs.org/release/download.php']);
    disp(info);
    msgbox(info,'Update','Help');
end
