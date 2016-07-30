function ea_checkinstall(cmd)

earoot=ea_getearoot;
switch cmd
    case 'bigbrain'
        if ~exist([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii'],'file')
            disp('BigBrain not installed. Downloading...')
            websave([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz'],'http://www.lead-dbs.org/release/data/bigbrain_2015_100um_bb.nii.gz');
            gunzip([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz']);
            try delete([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz']); end
            try delete([earoot,'templates',filesep,'bigbrain_2015_100um_bb.nii.gz.html']); end
        end
end