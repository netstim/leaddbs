function list=ea_bidsfiles2list(options)
list=cell(0);

transformfiles=ea_gettransformfiles(options);
if exist(transformfiles.forward,'file')
    list{end+1}=transformfiles.forward;
end
if exist(transformfiles.inverse,'file')
    list{end+1}=transformfiles.inverse;
end

procSteps={'preproc','coreg','norm','brainshift','recon','preopAnat','postopAnat'};

for procStep=1:length(procSteps)
    while 1
        thisfield=options.subj.(procSteps{procStep});
        path=cell(0);
        while ~ischar(thisfield)
            fn=fieldnames(thisfield);
            if isempty(fn)
                break
            end
            thisfield=thisfield.(fn{1});
            path{end+1}=fn{1};
        end
        % remove entry from struct
        switch length(path)
            case 0
                break
            case 1
                if ischar(thisfield) && exist(thisfield,'file')
                    list{end+1}=thisfield;
                end
                options.subj.(procSteps{procStep})=rmfield(options.subj.(procSteps{procStep}),path{1});
            case 2
                if ischar(thisfield) && exist(thisfield,'file')
                    list{end+1}=thisfield;
                end
                options.subj.(procSteps{procStep}).(path{1})=rmfield(options.subj.(procSteps{procStep}).(path{1}),path{2});
            case 3
                if ischar(thisfield) && exist(thisfield,'file')
                    list{end+1}=thisfield;
                end
                options.subj.(procSteps{procStep}).(path{1}).(path{2})=rmfield(options.subj.(procSteps{procStep}).(path{1}).(path{2}),path{3});
            case 4
                if ischar(thisfield) && exist(thisfield,'file')
                    list{end+1}=thisfield;
                end
                options.subj.(procSteps{procStep}).(path{1}).(path{2}).(path{3})=rmfield(options.subj.(procSteps{procStep}).(path{1}).(path{2}).(path{3}),path{4});
        end
    end
end
