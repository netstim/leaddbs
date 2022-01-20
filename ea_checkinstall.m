function [success, commands] = ea_checkinstall(cmd,checkonly,robot,prefs)
success=1;
if ~exist('checkonly','var')
    checkonly=0;
end
if ~exist('robot','var')
    robot=0;
end
if ~exist('prefs','var')
    prefs = [];
end

python_envs = dir(fullfile(ea_getearoot, 'classes', 'conda_utils', 'environments', '*.yml'));
python_envs = cellfun(@(x) x(1:end-4), {python_envs.name}', 'UniformOutput', false);

switch cmd
    case 'list' % simply return list of installable datasets
        success={'Redownload data files'
                 'Install development version of Lead'
                 'Python Environment (beta)'
                 python_envs
                 '2009b Nonlinear Flip Transform'
                 'Allan Institute Genetics Database'
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
                  'pyenv'
                  python_envs
                  'nlinflip'
                  'allengenetics'
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
    case 'allengenetics'
        checkf=[ea_space,'genetics'];
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
            success=ea_downloadasset('Allan Institute Genetics Database',...
                [ea_space,'genetics.zip'],...
                'allengenetics');
        else
            disp('Allan Institute Genetics Database is installed.')
        end
    case '7tcgflash'
        space = ea_space;
        if exist([space,'backdrops',filesep,'7T_Flash_Horn_2018.mat'], 'file')
            movefile([space,'backdrops',filesep,'7T_Flash_Horn_2018.mat'], ...
                     [space,'backdrops',filesep,'7T_Flash_Horn_2019.mat']);

            fid = fopen([space,'backdrops',filesep,'backdrops.txt'], 'r');
            text = fread(fid, inf, '*char')';
            fclose(fid);
            text = strrep(text, '7T_Flash_Horn_2018.mat', '7T_Flash_Horn_2019.mat');
            fid = fopen([space,'backdrops',filesep,'backdrops.txt'], 'w');
            fwrite(fid, text);
            fclose(fid);
        end

        checkf=[space,'backdrops',filesep,'7T_Flash_Horn_2019.mat'];
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
            ea_mkdir([space,'backdrops']);
            success=ea_downloadasset('7T Cardiac Gated Flash MRI',...
                [ea_space,'backdrops',filesep,'7T_Flash_Horn_2019.mat'],...
                '7tcgflash');
            m = matfile([space,'backdrops',filesep,'7T_Flash_Horn_2019.mat'],'Writable',true);
            m.fname = [space,'backdrops',filesep,'7T_Flash_Horn_2019.nii'];
            fid=fopen([space,'backdrops',filesep,'backdrops.txt'],'a');
            fprintf(fid,'%s %s\n','7T_Flash_Horn_2019.mat','7T_Cardiac_Gated_Flash_MRI_(Horn_2019)');
            fclose(fid);
        else
            disp('7T Cardiac Gated FLASH MRI (Backdrop visualization) is installed.')
        end
    case '7tev100um'
        space = ea_space;
        checkf=[space,'backdrops',filesep,'7T_100um_Edlow_2019.mat'];
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
            ea_mkdir([space,'backdrops']);
            success=ea_downloadasset('7T Ex Vivo 100um Brain Atlas',...
                [ea_space,'backdrops',filesep,'7T_100um_Edlow_2019.mat'],...
                '7tev100um');
            m = matfile([space,'backdrops',filesep,'7T_100um_Edlow_2019.mat'],'Writable',true);
            m.fname = [space,'backdrops',filesep,'7T_100um_Edlow_2019.nii'];
            fid=fopen([space,'backdrops',filesep,'backdrops.txt'],'a');
            fprintf(fid,'%s %s\n','7T_100um_Edlow_2019.mat','7T_Ex_Vivo_100um_Brain_Atlas_(Edlow_2019)');
            fclose(fid);
        else
            disp('7T Ex Vivo 100um Brain Atlas (Backdrop visualization) is installed.')
        end
    case 'macroscalehc'
        space=ea_space;
        checkf=[space,'atlases',filesep,'Macroscale Human Connectome Atlas (Yeh 2018)',filesep,'atlas_index.mat'];
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
                [space,'atlases',filesep,'Macroscale_Human_Connectome_Atlas_Yeh_2018.zip'],...
                'macroscalehc');
        else
            disp('Macroscale Human Connectome Atlas is installed.')
        end
    case 'bigbrain'
        space=ea_space;
        checkf=[space,'bigbrain_2015_100um_bb.nii'];
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
                [space,'bigbrain_2015_100um_bb.nii.gz'],...
                'bigbrain');
        else
            disp('BigBrain is installed.')
        end
    case 'groupconnectome2013'
        checkf=[ea_getconnectomebase('dmri',prefs),'Groupconnectome (Horn 2013) full',filesep,'data.mat'];
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
                [ea_getconnectomebase('dmri',prefs),'groupconnectome2013.zip'],...
                'group2013');
        else
            disp('Structural Group Connectome (Horn 2013) is installed.')
        end
    case 'groupconnectome2016'
        checkf=[ea_getconnectomebase('dmri',prefs),'Gibbsconnectome_169 (Horn 2016)',filesep,'data.mat'];
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
                [ea_getconnectomebase('dmri',prefs),'groupconnectome2016.zip'],...
                'group2016');
        else
            disp('Structural Group Connectome (Horn 2016) is installed.')
        end
    case 'groupconnectome2017'
        checkf=[ea_getconnectomebase('dmri',prefs),'HCP_MGH_32fold_groupconnectome (Horn 2017)',filesep,'data.mat'];
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
                [ea_getconnectomebase('dmri',prefs),'groupconnectome2017.zip'],...
                'group2017');
        else
            disp('Structural Group Connectome (Horn 2017) is installed.')
        end
    case 'groupconnectome_ppmi2017'
        checkf=[ea_getconnectomebase('dmri',prefs),'PPMI_85 (Ewert 2017)',filesep,'data.mat'];
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
                [ea_getconnectomebase('dmri',prefs),'groupconnectome_ppmi2017.zip'],...
                'group2017_ppmi');
        else
            disp('Structural Group Connectome (Ewert 2017) is installed.')
        end
    case 'fgroupconnectome_ppmi2017'
        checkf=[ea_getconnectomebase('fmri',prefs),'PPMI 74_15 (Horn 2017)',filesep,'dataset_info.mat'];
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
                [ea_getconnectomebase('fmri',prefs),'fgroupconnectome_ppmi2017.zip'],...
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
        
        if any(strcmp(cmd,python_envs))
            py_env = ea_conda_env(cmd);
            if ~checkonly
                if ~ea_conda.is_installed
                    ea_conda.install;
                end
                py_env.force_create;
            end
            success = py_env.is_created;
        end
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
    catch
        fsize = 0;
    end
    
    if fsize ~= 0
        fprintf('Downloading %s with a size of %.2f GB\nFilename: %s\n', assetname, fsize*1e-9, destination);
    else
        fprintf('Downloading %s\nFilename: %s\n', ea_getspace, destination);
        disp('This could take a while...');
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
                disp('Parallel toolbox not detected, downloading via websave...')
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
