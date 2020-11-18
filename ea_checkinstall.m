function [success, commands] = ea_checkinstall(cmd,checkonly,robot)
success=1;
if ~exist('checkonly','var')
    checkonly=0;
end
if ~exist('robot','var')
    robot=0;
end

switch cmd
    case 'list' % simply return list of installable datasets
        success={'Redownload data files'
                 'Install development version of Lead'
                 '2009b Nonlinear Flip Transform'
                 '7T Cardiac Gated FLASH MRI (Backdrop visualization)'
                 '7T Ex Vivo 100um Brain Atlas (Backdrop visualization)'
                 'Macroscale Human Connectome Atlas (Yeh 2018)'
                 'Structural group connectome 20 subjects Gibbs-tracker (Horn 2013)'
                 'Structural group connectome 169 NKI subjects Gibbs-tracker (Horn 2016)'
                 'Structural group connectome 32 Adult Diffusion HCP subjects GQI v1.1 (Horn 2017)'
                 'Structural group connectome 85 PPMI PD-patients GQI v1.1 (Ewert 2017)'
                 'Functional group connectome 74 PPMI PD-patients, 15 controls (Horn 2017)'};

        commands={'leaddata'
                  'hotfix'
                  'nlinflip'
                  '7tcgflash'
                  '7tev100um'
                  'macroscalehc'
                  'groupconnectome2013'
                  'groupconnectome2016'
                  'groupconnectome2017'
                  'groupconnectome_ppmi2017'
                  'fgroupconnectome_ppmi2017'};
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
    case 'nlinflip'
        checkf=[ea_space,'fliplr'];
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
            success=ea_downloadasset('FlipLR',...
                [ea_space,'fliplr.zip'],...
                'fliplr');
        else
            disp('2009b asym LR flip transform is installed.')
        end
    case '7tcgflash'
        if exist([ea_space,'backdrops',filesep,'7T_Flash_Horn_2018.mat'], 'file')
            movefile([ea_space,'backdrops',filesep,'7T_Flash_Horn_2018.mat'], ...
                     [ea_space,'backdrops',filesep,'7T_Flash_Horn_2019.mat']);

            fid = fopen([ea_space,'backdrops',filesep,'backdrops.txt'], 'r');
            text = fread(fid, inf, '*char')';
            fclose(fid);
            text = strrep(text, '7T_Flash_Horn_2018.mat', '7T_Flash_Horn_2019.mat');
            fid = fopen([ea_space,'backdrops',filesep,'backdrops.txt'], 'w');
            fwrite(fid, text);
            fclose(fid);
        end

        checkf=[ea_space,'backdrops',filesep,'7T_Flash_Horn_2019.mat'];
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
            ea_mkdir([ea_space,'backdrops']);
            success=ea_downloadasset('7T Cardiac Gated Flash MRI',...
                [ea_space,'backdrops',filesep,'7T_Flash_Horn_2019.mat'],...
                '7tcgflash');
            m = matfile([ea_space,'backdrops',filesep,'7T_Flash_Horn_2019.mat'],'Writable',true);
            m.fname = [ea_space,'backdrops',filesep,'7T_Flash_Horn_2019.nii'];
            fid=fopen([ea_space,'backdrops',filesep,'backdrops.txt'],'a');
            fprintf(fid,'%s %s\n','7T_Flash_Horn_2019.mat','7T_Cardiac_Gated_Flash_MRI_(Horn_2019)');
            fclose(fid);
        else
            disp('7T Cardiac Gated FLASH MRI (Backdrop visualization) is installed.')
        end
    case '7tev100um'
        checkf=[ea_space,'backdrops',filesep,'7T_100um_Edlow_2019.mat'];
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
            ea_mkdir([ea_space,'backdrops']);
            success=ea_downloadasset('7T Ex Vivo 100um Brain Atlas',...
                [ea_space,'backdrops',filesep,'7T_100um_Edlow_2019.mat'],...
                '7tev100um');
            m = matfile([ea_space,'backdrops',filesep,'7T_100um_Edlow_2019.mat'],'Writable',true);
            m.fname = [ea_space,'backdrops',filesep,'7T_100um_Edlow_2019.nii'];
            fid=fopen([ea_space,'backdrops',filesep,'backdrops.txt'],'a');
            fprintf(fid,'%s %s\n','7T_100um_Edlow_2019.mat','7T_Ex_Vivo_100um_Brain_Atlas_(Edlow_2019)');
            fclose(fid);
        else
            disp('7T Ex Vivo 100um Brain Atlas (Backdrop visualization) is installed.')
        end
    case 'macroscalehc'
        checkf=[ea_space,'atlases',filesep,'Macroscale Human Connectome Atlas (Yeh 2018)',filesep,'atlas_index.mat'];
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
            success=ea_downloadasset('Macroscale Human Connectome Atlas',...
                [ea_space,'atlases',filesep,'Macroscale_Human_Connectome_Atlas_Yeh_2018.zip'],...
                'macroscalehc');
        else
            disp('Macroscale Human Connectome Atlas is installed.')
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
            disp('BigBrain is installed.')
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
            success=ea_downloadasset('Structural Group Connectome (Horn 2013)',...
                [ea_getconnectomebase('dmri'),'groupconnectome2013.zip'],...
                'group2013');
        else
            disp('Structural Group Connectome (Horn 2013) is installed.')
        end
    case 'groupconnectome2016'
        checkf=[ea_getconnectomebase('dmri'),'Gibbsconnectome_169 (Horn 2016)',filesep,'data.mat'];
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
            success=ea_downloadasset('Structural Group Connectome (Horn 2016)',...
                [ea_getconnectomebase('dmri'),'groupconnectome2016.zip'],...
                'group2016');
        else
            disp('Structural Group Connectome (Horn 2016) is installed.')
        end
    case 'groupconnectome2017'
        checkf=[ea_getconnectomebase('dmri'),'HCP_MGH_32fold_groupconnectome (Horn 2017)',filesep,'data.mat'];
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
            success=ea_downloadasset('Structural Group Connectome (Horn 2017)',...
                [ea_getconnectomebase('dmri'),'groupconnectome2017.zip'],...
                'group2017');
        else
            disp('Structural Group Connectome (Horn 2017) is installed.')
        end
    case 'groupconnectome_ppmi2017'
        checkf=[ea_getconnectomebase('dmri'),'PPMI_85 (Ewert 2017)',filesep,'data.mat'];
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
            success=ea_downloadasset('Structural Group Connectome (Ewert 2017)',...
                [ea_getconnectomebase('dmri'),'groupconnectome_ppmi2017.zip'],...
                'group2017_ppmi');
        else
            disp('Structural Group Connectome (Ewert 2017) is installed.')
        end
    case 'fgroupconnectome_ppmi2017'
        checkf=[ea_getconnectomebase('fmri'),'PPMI 74_15 (Horn 2017)',filesep,'dataset_info.mat'];
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
            success=ea_downloadasset('Functional Group Connectome (Horn 2017)',...
                [ea_getconnectomebase('fmri'),'fgroupconnectome_ppmi2017.zip'],...
                'fgroup2017_ppmi');
        else
            disp('Functional Group Connectome (Horn 2017) is installed.')
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
    success=1;
else
    downloadurl = 'https://www.lead-dbs.org/release/download.php';
    success=1;
    
    if ~exist(fileparts(destination), 'dir')
        mkdir(fileparts(destination));
    end
    
    % get the file size and display it to inform user of size
    try
        fsize = ea_getassetfilesize(id);    % get the filesize in bytes from the server
        fprintf('Downloading %s with a size of %.2f GB\nFilename: %s\n', assetname, fsize*1e-9, destination);
    catch
        fsize = 0;
        fprintf('Downloading %s\nFilename: %s\n', assetname, destination);
    end
    
    % first see if parallel toolbox is installed and can be utilized to
    % download in the backgroung
    if license('test', 'Distrib_Computing_Toolbox') % check if parallel toolbox is installed
        disp('Parallel toolbox detected, downloading via background worker')
        try
            ea_downloadasset_parallel(downloadurl, assetname, destination, id, fsize);
        catch
            try
                delete(findall(0,'type','figure','tag','TMWWaitbar'));      % delete waitbar if error occured
                disp('Download using parallel toolbox failed, trying websave...')
                webopts=weboptions('Timeout',Inf);
                websave(destination,downloadurl,'id',id,webopts);
            catch
                try
                    disp('''websave'' failed, trying ''urlwrite''.');
                    urlwrite([downloadurl,'?id=',id],destination,'Timeout',Inf);
                catch
                    disp('''urlwrite'' failed.');
                    success=0;
                end
            end
        end
    else
        try
            disp('Trying ''websave''.');
            webopts=weboptions('Timeout',Inf);
            % uncomment this if you encounter problems with certificate validation
            % (see https://www.mathworks.com/matlabcentral/answers/400086-how-to-read-data-form-website-url)
            %webopts.CertificateFilename=('');
            websave(destination,downloadurl,'id',id,webopts);
        catch
            try
                disp('''websave'' failed, trying ''urlwrite''.');
                urlwrite([downloadurl,'?id=',id],destination,'Timeout',Inf);
            catch
                disp('''urlwrite'' failed.');
                success=0;
            end
        end
    end

    [loc,~,ext] = fileparts(destination);
    if success
        disp(['Installing ',assetname,'...'])
        try
            if strcmp(ext,'.gz')
                gunzip(destination, loc);
                ea_delete(destination);
            elseif strcmp(ext,'.zip')
                unzip(destination, loc);
                ea_delete(destination);
            end
        catch
            disp('Installation failed.');
            success=0;
        end
    end
    if ~success
        fprintf(['\nDownload / install error! You may try to download the file manually from:\n',...
                 '%s\nand then extract it into %s.\n\n'], [downloadurl,'?id=',id], loc);
        msgbox('Please check the command window for more information.','Download error!','Error')
        return
    end
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
