function success=ea_hotfix(varargin)
% lead hotfix
success=0;
local = ea_getvsn('local');
web = ea_getvsn('web');

% if local version unknow or outdated, skip applying hotfix
hotfix = 0;
if strcmp(local,web) % local version is the latest release, applicable for hotfix
    hotfix = 1;
end

earoot=ea_getearoot;

updurl = 'http://www.lead-dbs.org/release/download.php';
if hotfix
    try
        disp('*** Updating LEAD. Please do not quit MATLAB.');
        ea_delete([earoot,'tmp',filesep,'hotfix']);
        if ~exist([earoot,'tmp'] ,'dir')
            mkdir([earoot,'tmp']);
        end
        disp('Downloading hotfix...');
        try
            webopts=weboptions('Timeout',5);
            websave([earoot,'tmp',filesep,'hotfix.zip'],updurl,'id','hotfix',webopts);
        catch
            try
                urlwrite([updurl,'?id=hotfix'],[earoot,'tmp',filesep,'hotfix.zip'],'Timeout',5);
            catch
                info='Download error! Please retry later.';
                disp(info);
                msgbox(info,'Update','Error')
                return
            end
        end

        disp('Extracting code...');
        try
            unzip([earoot,'tmp',filesep,'hotfix.zip'],[earoot,'tmp',filesep,'hotfix']);
        catch
            system(['unzip -q ',earoot,'tmp',filesep,'hotfix.zip -d ', earoot,'tmp',filesep,'hotfix']);
        end
        delete([earoot,'tmp',filesep,'hotfix.zip']);

        disp('Deleting outdated code...');
        try
            if exist([earoot,'tmp',filesep,'hotfix',filesep,'DELETE'], 'file')
                dfid = fopen([earoot,'tmp',filesep,'hotfix',filesep,'DELETE']);
                dels=textscan(dfid,'%s');
                fclose(dfid);
                for f=1:length(dels{1})
                    if isdir([earoot,dels{1}{f}])
                        rmdir([earoot,dels{1}{f}],'s')
                    else
                        delete([earoot,dels{1}{f}])
                    end
                end
                delete([earoot,'tmp',filesep,'hotfix',filesep,'DELETE'])
            end
        catch
            disp('Error while deleting some files. You may ignore this.');
        end

        disp('Copying new code...');
        copyfile([earoot,'tmp',filesep,'hotfix',filesep,'*'],earoot,'f');
        disp('Cleaning up...');
        rmdir([earoot,'tmp'],'s')
        disp('Done.');

        disp('Restarting LEAD.');
        close all force
        lead;
        success=1;
        disp('*** Update finished.');
    catch
        info=sprintf(['Hotfix does not exist or failed to apply hotfix!\n',...
                      'Please wait for the next release.']);
        disp(info);
        msgbox(info,'Update','Error');
    end
else
    info=sprintf('Local version is not applicable for hotfix!');
    disp(info);
    msgbox(info,'Update','Help');
end
