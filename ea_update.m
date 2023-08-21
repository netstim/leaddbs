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

earoot=ea_getearoot;

updurl = 'https://www.lead-dbs.org/release/download.php';
if update
    try
        if update==1 % full update
            id='leaddbs';
        elseif update==2 % incremental update
            id=['updates_',strrep(local,'.',''),'-',strrep(web,'.','')];
        end

        disp('*** Updating LEAD. Please do not quit MATLAB.');
        disp('Downloading code...');
        if ~exist([earoot,'tmp'],'dir')
            mkdir([earoot,'tmp']);
        end
        try
            webopts=weboptions('Timeout',Inf);
            websave([earoot,'tmp',filesep,'updates.zip'],updurl,'id',id,webopts);
        catch
            try
                urlwrite([updurl,'?id=',id],[earoot,'tmp',filesep,'updates.zip'],'Timeout',Inf);
            catch
                if update==1
                    ea_cprintf('CmdWinWarnings', ['\nDownload error! You may try to download Lead-DBS manually from:\n',...
                             '%s_pcloud or\n%s_onedrive\n\n'], [updurl,'?id=',id], [updurl,'?id=',id]);

                elseif update==2
                    ea_cprintf('CmdWinWarnings', ['\nDownload error! You may try to download the update package manually from:\n',...
                             '%s\nand then extract it into Lead-DBS installation folder.\n\n'], [updurl,'?id=',id]);
                end
                msgbox('Please check the command window for more information.','Download error!','Error')
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

        try
            % delete files during incremental updating
            if update==2 && exist([earoot,'tmp',filesep,id,filesep,'DELETE'], 'file')
                disp('Deleting outdated code...');
                dels = readlines([earoot,'tmp',filesep,id,filesep,'DELETE']);
                ea_delete(cellstr(dels));
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

        disp('Restarting LEAD.');
        close all force
        lead;
        disp('*** Update finished.');

        if getResponseCode(openConnection(java.net.URL([updurl,'?id=updates_data']))) == 200
            info='Start updating LEAD data...';
            disp(info);
            msgbox(info,'Data Update Found','Help');
            ea_update_data();
        end
    catch
        info=sprintf(['Failed to update!\n',...
                      'Alternatively, you can download the latest version from: https://www.lead-dbs.org/release/download.php']);
        disp(info);
        msgbox(info,'Update','Error');
    end
else
    info=sprintf(['LEAD aleady up-to-date!\n',...
                  'Alternatively, you can re-download the latest version from: https://www.lead-dbs.org/release/download.php']);
    disp(info);
    msgbox(info,'Update','Help');
end
