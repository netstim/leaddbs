function [options,presentfiles]=ea_assignpretra(options)

directory=[options.root,options.patientname,filesep];

presfiles=dir([directory,options.prefs.prenii_searchstring]);
for p=1:length(presfiles)
   pfcell{p}=presfiles(p).name; 
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
pfcell_priority=[pfcell_priority,pfcell];


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
options.primarytemplate=ptemps;


presentfiles=pfcell_priority;


