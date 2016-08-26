function success=ea_checkinstall(cmd,force)
success=1;
earoot=ea_getearoot;
if ~exist('force','var')
    force=0;
end
switch cmd
    case 'bigbrain'
        if ~exist([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii'],'file') || force
            disp('BigBrain not installed. Downloading...')
            try
            websave([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz'],'http://www.lead-dbs.org/release/data/bigbrain_2015_100um_bb.nii.gz');
            gunzip([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz']);
            catch
                success=0;
            end
            try delete([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz']); end
            try delete([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz.html']); end
        else
            disp('BigBrain is installed.')
        end
    case 'macaque'
        if ~exist([earoot,'toolbox',filesep,'macaque'],'file') || force
            disp('Macaque toolbox not installed. Downloading...')
            try
            websave([earoot,'toolbox',filesep,'macaque.zip'],'http://www.lead-dbs.org/release/download.php?id=macaque');
            unzip([earoot,'toolbox',filesep,'macaque.zip']);
            catch
                success=0;
            end
            try delete([earoot,'toolbox',filesep,'macaque.zip']); end
            try delete([earoot,'toolbox',filesep,'macaque.zip.html']); end
        else
            disp('Macaque toolbox is installed.')
        end
    case 'groupconnectome2013'
        if ~exist([earoot,ea_getconnectomebase('dmri'),filesep,'Groupconnectome (Horn 2013) full.mat'],'file') || force
            try
            websave([earoot,ea_getconnectomebase('dmri'),filesep,'groupconnectome2013.zip'],'http://www.lead-dbs.org/release/download.php?id=group');
            unzip([earoot,ea_getconnectomebase('dmri'),filesep,'groupconnectome2013.zip']);
            catch
                success=0;
            end
            try delete([earoot,ea_getconnectomebase('dmri'),filesep,'groupconnectome2013.zip']); end
            try delete([earoot,ea_getconnectomebase('dmri'),filesep,'groupconnectome2013.zip.html']); end    
        else
            disp('Group Connectome (Horn 2013) is installed.')
        end
    otherwise
        success=0;
end