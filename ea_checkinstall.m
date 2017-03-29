function [success,commands]=ea_checkinstall(cmd,checkonly,robot)
success=1;
if ~exist('checkonly','var')
    checkonly=0;
end
if ~exist('robot','var')
    robot=0;
end

switch cmd
    case 'list' % simply return list of installable datasets
        success={'Redownload Data files','Apply Hotfix','Big Brain 100um subcortical (Amunts 2013)','Structural group connectome (Horn 2013)'};
        commands={'leaddata','hotfix','bigbrain','groupconnectome2013'};
    case 'leaddata'
        checkf=[ea_space,'bb.nii'];
        force=ea_alreadyinstalled(checkf,checkonly,robot);
        if checkonly
            success=~force;
            return;
        end
        if force==-1
            success=-1;
            return;
        end

        if ~exist(checkf,'file') || force
            success=ea_downloadasset('Lead Datafiles',...
                [],...
                '');
        else
            disp('Lead datafiles is installed.')
        end

    case 'bigbrain'
        checkf=[ea_space,'bigbrain_2015_100um_bb.nii'];
        force=ea_alreadyinstalled(checkf,checkonly,robot);
        if checkonly
            success=~force;
            return;
        end
        if force==-1
            success=-1;
            return;
        end

        if ~exist(checkf,'file') || force
            success=ea_downloadasset('Bigbrain 100um subcortical',...
                [ea_space,'bigbrain_2015_100um_bb.nii.gz'],...
                'bigbrain');
        else
            %disp('BigBrain is installed.')
        end
    case 'groupconnectome2013'
        checkf=[ea_getconnectomebase('dmri'),'Groupconnectome (Horn 2013) full',filesep,'data.mat'];
        force=ea_alreadyinstalled(checkf,checkonly,robot);
        if checkonly
            success=~force;
            return;
        end
        if force==-1
            success=-1;
            return;
        end

        if ~exist(checkf,'file') || force
            success=ea_downloadasset('structural group connectome (Horn 2013)',...
                [ea_getconnectomebase('dmri'),'groupconnectome2013.zip'],...
                'group');
        else
            disp('Group Connectome (Horn 2013) is installed.')
        end
    case 'hotfix'
        if ~checkonly
            success=ea_hotfix;
        end
    otherwise
        success=0;
end


function success=ea_downloadasset(assetname,destination,id)

if strcmp(assetname,'Lead Datafiles')
    ea_update_data('full');

else
    downloadurl = 'http://www.lead-dbs.org/release/download.php';
    success=1;
    disp(['Downloading ',assetname,'...'])
    if ~exist(fileparts(destination), 'dir')
        mkdir(fileparts(destination));
    end
    try
        webopts=weboptions('Timeout',5);
        websave(destination,downloadurl,'id',id,webopts);
    catch
        try
            urlwrite([downloadurl,'?id=',id],destination,'Timeout',5);
        catch
            success=0;
        end
    end

    if success
        [loc,~,ext] = fileparts(destination);
        if strcmp(ext,'.gz')
            gunzip(destination, loc);
        elseif strcmp(ext,'.zip')
            unzip(destination, loc);
        end
    end

    ea_delete(destination);
end


function force=ea_alreadyinstalled(checkf,checkonly,robot)
if ~exist(checkf,'file') % then file not there, should install anyways.
    force=1;
    return
end
if checkonly % return here.
    force=0;
    return
end
if ~robot
    choice = questdlg('This dataset seems already to be installed. Do you wish to re-download it?', ...
        'Redownload Dataset?', ...
        'Yes','No','No');
    if strcmp(choice,'Yes')
        force=1;
    else
        force=-1;
    end
else
    force=0;
end
