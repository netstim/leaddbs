function ea_MergeGroupAnalysis(inputLGfiles,datasetDir,guid)
    %This function merges multiple BIDS lead group files.
    %output path is a BIDS directory. All lead group files, including the
    %merged file will be saved here.
    %GUID is an optional argument.
    %if no GUID is provided, then the GUID of the first LG file will be
    %assigned to the merged GUID
    
    %initialize M
    M.patient.list = {};
    M.clinical.labels = {};
    M.groupdir = datasetDir;
    M.patient.group = [];
    %figure out if the lead group files are in classic or in BIDS
    for LGfiles=1:length(inputLGfiles)
        cohortM = load(inputLGfiles{LGfiles});
        if ~contains(cohortM.M.root,'derivatives/leadgroup')
            ea_warning("It looks like some of the lead groups are not in BIDS format. Please convert them to BIDS before you use this function");
        end
        %now merge the patients list
        M.patient.list = [M.patient.list;cohortM.M.patient.list];
        M.patient.group = [M.patient.group;repmat(LGfiles,length(cohortM.M.patient.list),1)];
        M.clinical.labels = [M.clinical.labels,cohortM.M.clinical.labels];
        %merging structs is a bit different
        if LGfiles == 1
            M.stats = cohortM.M.stats;
            M.S = cohortM.M.S;
            M.elstruct = cohortM.M.elstruct;
            %take the following from the first cohort
            M.vilist = cohortM.M.vilist;
            M.fclist = cohortM.M.fclist;
            M.ui = cohortM.M.ui;
            if ~exist('guid','var')
                M.guid = cohortM.M.guid;
            else
                M.guid = guid;
            end
            M.vatmodel = cohortM.M.vatmodel; %this has to be the same across all cohorts!!! the vatmodel of the first LG file will be used!!
            M.root = datasetDir;
        else
            M.stats = [M.stats,cohortM.M.stats];
            M.S = [M.S,cohortM.M.S];
            M.elstruct = [M.elstruct,cohortM.M.elstruct];
           
        end
       
    end
    %groups
    M.groups.group = 1:length(inputLGfiles);
    M.groups.color = ea_color_wes('lifeaquatic',length(inputLGfiles));
    M.groups.colorschosen = length(inputLGfiles);
    M.clinical.vars = cell(1,length(M.clinical.labels));
    for i=1:length(M.clinical.vars)
        %initialize with nan
        M.clinical.vars{i} = nan(length(M.patient.list),1);
    end
    varidx = 1;
    %open the loop for filling the vars it again
    for LGfiles_iter = 1:length(inputLGfiles)
        cohortM = load(inputLGfiles{LGfiles_iter});
        cohortLength = length(cohortM.M.patient.list);
        if LGfiles_iter == 1
            from = 1;
            to = cohortLength;
        else
            from = to+1;
            to = to+cohortLength;
        end
        for j = 1:length(cohortM.M.clinical.vars)
            try
                M.clinical.vars{varidx}(from:to) = cohortM.M.clinical.vars{j};
                varidx = varidx + 1;
            catch ME
                if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
                    disp("Looks like the lead group files you input does not have correct number of variables. Please review and try again")
                end
               
            end

        end
    end
    %save the new M file
    [~,datasetName,~] = fileparts(datasetDir);
    output_filename = ['dataset-',datasetName,'_','analysis-',M.guid,'.mat'];
    fullpath = fullfile(datasetDir,'derivatives','leadgroup',M.guid);
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    save(fullfile(fullpath,output_filename),'M');
    disp("***Process done***");

end