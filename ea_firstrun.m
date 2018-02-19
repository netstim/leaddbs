function ea_firstrun(handles,options)

% check if a newer version is available..
local=ea_getvsn('local',1);
web=ea_getvsn('web',1);
vcheck=0;
if strcmp(local,'Unknown')
    vcheck=1;
elseif ~strcmp(web,'Unknown')
    vcheck=(local<web);
end

try
    if vcheck
        set(handles.updatebutn,'BackgroundColor',[0.2,0.8,0.2]);
    else
        set(handles.updatebutn,'BackgroundColor',[0.94,0.94,0.94]);
    end
end

% if ~strcmp(handles.prod,'dbs_connectome')
%     try
%         webopts=weboptions('Timeout',5);
%         webread('http://www.lead-dbs.org/release/stats.php','id',handles.prod,'ver',['R',version('-release')],webopts);
%     catch
%         try
%             urlread(['http://www.lead-dbs.org/release/stats.php?id=',handles.prod,'&ver=R', version('-release')],'Timeout',5);
%         catch
%         end
%     end
% end


if ~isfield(options.prefs,'firstrun') % first run.
    fprintf(['Welcome to LEAD-DBS.\n \n',...
        'This seems to be your first run of the toolbox.\n',...
        'To start with your first reconstruction, create a folder for a patient and put either MR files \n',...
        'or CT files inside of it. \n',...
        'Naming of the files should be as the following: \n \n',...
        'MR images: \n',...
        '+ patient_directory \n',...
        '    + postop_tra.nii (for an unnormalized version of a transversal postoperative image, i.e. in native space). \n',...
        '    + postop_cor.nii (for an unnormalized version of a coronal postoperative image). \n',...
        '    + postop_sag.nii (for an unnormalized version of a coronal postoperative image). \n',...
        '    + anat_t1.nii (for an unnormalized version of a preoperative image). \n \n',...
        'If you only have one postoperative image, you should name that one postop_tra.nii and omit the rest. \n',...
        'CT images: \n',...
        '+ patient_directory \n',...
        '    + postop_ct.nii (for postoperative CT image). \n',...
        '    + anat_t1.nii (for preoperative MRI image). \n \n',...
        'Most naming-conventions can be modified by editing the file ea_prefs.m \n',...
        'However sticking to the given naming-conventions might make it easier e.g. to work on collaborations. \n \n',...
        'We hope that you enjoy using the Lead-DBS toolbox. \n \n',...
        'Any suggestions are more than welcome (andreas.horn@charite.de). \n'
        ]);

    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.m'],[ea_gethome,'.ea_prefs.m'], 'f');
    fid = fopen([ea_gethome,'.ea_prefs.m'],'a');
    fprintf(fid, '\n\nprefs.firstrun=''off'';\n');
    fclose(fid);

    % check dataset isntallation
    if ~exist([ea_space,'bb.nii'], 'file')
        fprintf(['\nIt seems that you don''t have LEAD dataset installed.\n' ...
                 'You can either install it via ''Install'' --> ''Redownload Data Files'' menu,\n' ...
                 'or download it from http://www.lead-dbs.org/release/download.php?id=data and extract it into LEAD folder.\n\n']);

        msg = sprintf(['It seems that you don''t have LEAD dataset installed.\n' ...
                      'Do you wish to download it now? Alternatively, you can also download it from http://goo.gl/yQz8bz and extract it into LEAD folder.']);
        choice = questdlg(msg, 'Download Dataset?', 'Yes', 'No', 'Yes');
        if strcmp(choice, 'Yes')
            ea_update_data('full');
        end
    end
end
