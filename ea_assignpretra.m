function [options,presentfiles]=ea_assignpretra(options)

directory=[options.root,options.patientname,filesep];

presfiles=dir([directory,options.prefs.prenii_searchstring]);
pfcell = {presfiles.name}';

% order the anatomical images in accordance with 'prefs.prenii_order'
prenii_order = cellfun(@(x) strrep(options.prefs.prenii_searchstring,'*',x), options.prefs.prenii_order, 'UniformOutput', 0);
[~,idx] = ismember(prenii_order, pfcell);
pfcell_priority = pfcell([nonzeros(idx)',setdiff(1:numel(pfcell),nonzeros(idx))]);

% determine primary template
if any(idx)
    options.primarytemplate = options.prefs.prenii_order{find(idx,1)};
else % could happen if neither T2, T1 or PD is present but only custom sequences are being used
    options.primarytemplate = 't2'; % default T2.
end

% define correct order:
cnt=1;
for priority=1:length(options.prefs.prenii_order)
    [check,ix]=ismember(strrep(options.prefs.prenii_searchstring,'*',options.prefs.prenii_order{priority}),pfcell);
    if check;
        pfcell_priority{cnt}=pfcell{ix};
        pfcell(ix)=[];
        if cnt==1 % also assign template suffix
            ptemps=options.prefs.prenii_order{priority};
        end
        cnt=cnt+1;
    end
end


if ~exist('ptempts','var') % could happen if neither T2, T1 or PD is present but only custom sequences are being used
    ptemps='_t2'; % default T2.
end
try
pfcell_priority=[pfcell_priority,pfcell];
catch
    pfcell_priority=pfcell; % if no single anatomy file inside folder matching the search string.
end

if isempty(pfcell_priority)
    warning(['No anatomy information found. Please put either ',options.prefs.rawpreniis{1},', ',  options.prefs.rawpreniis{2},' or ', options.prefs.rawpreniis{3}, ' into subject folder.']);
end

% anat preprocess, only do once.
% a small hidden file '.pp' inside patient folder will show this has been done before.
if ~exist([directory,'.pp'],'file') && ~exist([directory,'ea_normmethod_applied.mat'],'file')
    for fi=1:length(pfcell_priority)
        % apply reorient/crop and biasfieldcorrection
        ea_anatpreprocess([directory,pfcell_priority{fi}]);
    end
    try
        fs=fopen([directory,'.pp'],'w');
        fprintf(fs,'%s','anat preprocess done');
        fclose(fs);
    end
end

options.prefs.prenii_unnormalized=pfcell_priority{1};

presentfiles=pfcell_priority;
