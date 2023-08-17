
space = ea_space;
root = ea_getearoot;

if strcmp(ea_getspace,'MNI152NLin2009bAsym')
    if exist([root,'templates',filesep,'mni_hires_t1.nii'],'file');
        try movefile([root,'templates'], [root,'templates_temp']); end
        try mkdir([root,'templates',filesep,'space',filesep]); end

        try movefile([root,'templates_temp'],[root,'templates',filesep,'space',filesep,ea_getspace]); end % do not use space in this line since wont be defined and result in error.
        try movefile([root,'atlases'],[space,'atlases']); end

        try movefile([space,'mni_hires_t1.nii'],[space,'t1.nii']); end
        try movefile([space,'mni_hires_t2.nii'],[space,'t2.nii']); end
        try movefile([space,'mni_hires_pd.nii'],[space,'pd.nii']); end
        try movefile([space,'mni_hires_fa.nii'],[space,'fa.nii']); end
        try movefile([space,'mni_hires_bb.nii'],[space,'bb.nii']); end
        try movefile([space,'mni_hires_c1mask.nii'],[space,'c1mask.nii']); end
        try movefile([space,'mni_hires_c2mask.nii'],[space,'c2mask.nii']); end
        try movefile([space,'TPM_2009b.nii'],[space,'TPM.nii']); end
        try movefile([space,'mni_hires_distal.nii'],[space,'atlas.nii']); end
        try movefile([space,'mni_wires.mat'],[space,'wires.mat']); end
    end

    if ~exist([space,'spacedef.mat'],'file')
        try spacedef=ea_gendefspacedef; end
        try save([space,'spacedef.mat'],'spacedef'); end

    end

    if exist([space,'distal.nii'],'file')
        try movefile( [space,'distal.nii'],[space,'atlas.nii']); end
    end

    if exist([space,'TPM_Lorio_Draganski.nii'],'file')
        try movefile([space,'TPM_Lorio_Draganski.nii'],[root,'templates',filesep,'TPM_Lorio_Draganski.nii']); end
    end

    if exist([space,'electrode_contacts'],'dir')
        try movefile( [space,'electrode_contacts'],[root,'templates']); end
    end

    if exist([space,'electrode_models'],'dir')
        try movefile( [space,'electrode_models'],[root,'templates']); end
    end

    % remove square bracket from folder name
    folders = ea_regexpdir([space, 'atlases'], '^(?!\.).*', 0, 'dir');
    for i=1:length(folders)
        oldfolder = folders{i};
        newfolder = regexprep(oldfolder, '(.*) \[(.*)\]', '$1 - $2');
        if ~strcmp(oldfolder, newfolder)
            movefile(oldfolder, newfolder)
        end
    end

    % remove square bracket from file name
    files = ea_regexpdir([space, 'labeling'], '^(?!\.).*', 0, 'file');
    for i=1:length(files)
        oldfile = files{i};
        newfile = regexprep(oldfile,'\[(.*)\]', '- $1');
         if ~strcmp(oldfolder, newfolder)
            movefile(oldfile, newfile)
         end
    end
end

if ~exist([root,'templates'],'dir')
    prefs=ea_prefs;
    if exist([prefs.lc.datadir,'spacedefinitions'], 'dir')
        try movefile([prefs.lc.datadir,'spacedefinitions'], [root,'templates']); end
    elseif exist([root,'spacedefinitions'], 'dir')
        try movefile([root,'spacedefinitions'], [root,'templates']); end
    end
end

