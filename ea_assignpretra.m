function [options, presentfiles] = ea_assignpretra(options,allowgz)

if ~exist('allowgz','var')
    allowgz=0;
end

directory = fullfile(options.root, options.patientname);

if allowgz
    sstr = [options.prefs.prenii_searchstring(1:end-1),'*'];
else
    sstr = [options.prefs.prenii_searchstring];
end

presfiles = dir(fullfile(directory, sstr));
pfcell = {presfiles.name}';

% order the anatomical images in accordance with 'prefs.prenii_order'
prenii_order = cellfun(@(x) strrep(options.prefs.prenii_searchstring,'*',x), options.prefs.prenii_order, 'UniformOutput', 0);
[~,idx] = ismember(ea_stripext(prenii_order), ea_stripext(pfcell));
presentfiles = pfcell([nonzeros(idx)',setdiff(1:numel(pfcell),nonzeros(idx))]);

load([ea_space,'ea_space_def.mat']);

options.primarytemplate = spacedef.templates{1}; % default T1.

if isempty(presentfiles)
    warning(['No anatomy information found! Please put either ', ...
        prenii_order{1},', ',prenii_order{2},' or ',prenii_order{3},' into subject folder.']);
    return
end

ism = ismember(presentfiles,{'anat_STN.nii','anat_GPi.nii','anat_GPe.nii','anat_RN.nii'});
addtoend = presentfiles(ism);
presentfiles(ism) = [];
presentfiles = [presentfiles;addtoend];

% set prenii_unnormalized
options.prefs.prenii_unnormalized = presentfiles{1};

% determine primary template
if any(idx)
    anchor = options.prefs.prenii_order(find(idx,1));
    templateFound = 0;
    for t=1:size(spacedef.norm_mapping,1)
       if ismember(anchor, spacedef.norm_mapping{t,1})
           templateFound = 1;
           options.primarytemplate = spacedef.norm_mapping{t,2};
           break
       end
    end

    if ~templateFound
        options.primarytemplate = spacedef.misfit_template;
    end
end
