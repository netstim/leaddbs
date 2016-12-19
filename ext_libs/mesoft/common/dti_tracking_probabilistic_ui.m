% interface program for batch program
%
% Author: Susanne Schnell and Volkmar Glauche
% PC 25.07.2008

function out = dti_tracking_probabilistic_ui(P)
%load DTD, start and stop masks
dtdStruct = dtdstruct_read(P.filename{1});

switch char(fieldnames(P.seedchoice))
    % P.seedchoice always has exactly one field, therefore this
    % switch is correct
    case 'seedroi'
        if isfield(P.seedchoice.seedroi.seedmask,'seedname')
            check = P.seedchoice.seedroi.seedmask.seedname;
            if iscellstr(check)
                check = size(check,1); % not tested!
            else
                check = 1;
            end
        else
            check = P.seedchoice.seedroi.seedmask.seednumber;
        end
        [SeedMask errStr] = maskstruct_read(P.seedchoice.seedroi.seedfile{1});
        
        [test, errStr, orInfo]= maskstruct_query(SeedMask, 'testMRDat', dtdStruct);
        if test > 0
            if ~orInfo  % no valid orientation info available in maskStruct
                SeedMask= maskstruct_modify(SeedMask, 'mrStructProb', ...
                    dtdstruct_query(dtdStruct, 'mrStructProb'));
            end
        elseif test == 0
            [SeedMask, errStr, ok]= maskstruct_modify(SeedMask, 'coregister', dtdStruct);
            if ~ok
                error(['could not open maskStruct:' errStr]);
            end
        else
            error('could not open maskStruct');
        end                                                           
        
        if numel(check) == 1 && isfinite(check)
            if isfield(P.seedchoice.seedroi.seedmask,'seedname')
                [SeedMask, errStr]= maskstruct_query(SeedMask, 'getMaskVc', P.seedchoice.seedroi.seedmask.seedname);
            elseif isfield(P.seedchoice.seedroi.seedmask,'seednumber')
                [SeedMask, errStr]= maskstruct_query(SeedMask, 'getMaskVc', P.seedchoice.seedroi.seedmask.seednumber);
            end
            if ~isempty(errStr)
                error(errStr);
            end
            if isempty(SeedMask)
                error('Selected ROI as seed region is empty');
            end
        else
            % prepare jobs for multiple seeds ROIs                                   
            
            maxMaskNo = maskstruct_query(SeedMask,'maskNo');
            if isfield(P.seedchoice.seedroi.seedmask,'seednumber')
                if any(~isfinite(P.seedchoice.seedroi.seedmask.seednumber))
                    newseedmask = 1:maxMaskNo;
                else
                    newseedmask = P.seedchoice.seedroi.seedmask.seednumber;
                end
            elseif isfielf(P.seedchoice.seedroi.seedmask,'seedname')
                %% build up newseedmask?
                newseedmask = P.seedchoice.seedroi.seedmask.seedname;
            end
            % set format of output filenames
            [path,nam,ext] = fileparts(P.filename{1});
            if isfield(P.newprobfile,'auto')
                if length(nam) > 4 && strcmp(nam(end-3:end), '_DTD')
                    nam = nam(1:end-4);
                end
                nam = [nam '_MAP'];
                P.newprobfile = struct('out',struct('dir',{{path}},'fname',''));
            else
                [p1,nam,ext] = fileparts(P.newprobfile.out.fname);
            end
            namfmt = sprintf('%s%%0%dd.mat', nam, floor(log10(maxMaskNo)+1));
            % build job list
            nj = cell(size(newseedmask));
            for k = 1:numel(newseedmask)
                nj{k}.dtijobs.tracking.probabilistic = P;
                nj{k}.dtijobs.tracking.probabilistic.seedchoice.seedroi.seedmask.seednumber = newseedmask(k);
                nj{k}.dtijobs.tracking.probabilistic.newprobfile.out.fname = sprintf(namfmt, newseedmask(k));
            end
            cj = cfg_util('initjob', nj);
            cfg_util('run', cj);
            jout = cfg_util('getAllOutputs', cj);
            cfg_util('deljob', cj);
            % construct cell array of filenames from the cell of
            % individual job outputs.
            jout = [jout{:}];
            out.files = [jout.files]';
            % finished after all jobs have completed
            return;
        end
    case 'seedposition'
        SeedMask = P.seedchoice.seedposition;
end
if isfield(P.defroi,'defthresh')
    [faData errStr] = dtdstruct_query(dtdStruct, 'getFA');
    if ~isempty(errStr)
        error(errStr);
    end
    [trdData errStr] = dtdstruct_query(dtdStruct, 'getTrace');
    if ~isempty(errStr)
        error(errStr);
    end
    DefArea = (trdData.dataAy < P.defroi.defthresh.deftrace) & (faData.dataAy > P.defroi.defthresh.deffa);
else
    DefArea = maskstruct_read(P.defroi.defarea.deffile{1});
    
    
    [test, errStr, orInfo]= maskstruct_query(DefArea, 'testMRDat', dtdStruct);
    if test > 0
        if ~orInfo  % no valid orientation info available in maskStruct
            DefArea= maskstruct_modify(DefArea, 'mrStructProb', ...
                dtdstruct_query(dtdStruct, 'mrStructProb'));
        end
    elseif test == 0
        [DefArea, errStr, ok]= maskstruct_modify(DefArea, 'coregister', dtdStruct);
        if ~ok
            error(['could not open maskStruct:' errStr]);
        end
    else
        error('could not open maskStruct');
    end                                                           
              
    
    if isfield(P.defroi.defarea.defmask,'defnumber')
        [DefArea errStr] = maskstruct_query(DefArea, 'getMask',P.defroi.defarea.defmask.defnumber);
        if ~isempty(errStr)
            error(errStr);
        end
    else
        [DefArea errStr] = maskstruct_query(DefArea,'getMask',P.defroi.defarea.defmask.defname);
        if ~isempty(errStr)
            error(errStr);
        end
    end
end
if P.probrand == 2
    algName = 'ProbRandExt';
    [probStruct errStr] = probstruct_op(strcat('apply_', algName),dtdStruct,double(SeedMask),DefArea,P.fibrelength,P.nowalks,P.exponent,P.revisits,400);
else
    algName = 'ProbRand';
    [probStruct errStr] = probstruct_op(strcat('apply_', algName),dtdStruct,double(SeedMask),DefArea,P.fibrelength,P.nowalks,P.exponent,[],[]);
end
if ~isempty(errStr)
    error(errStr);
end
if isempty(probStruct)
    error(lasterror);
end
% save file
[path,nam,ext] = fileparts(P.filename{1});
if isfield(P.newprobfile,'auto')
    if length(nam) > 4 && strcmp(nam(end-3:end), '_DTD')
        nam = nam(1:end-4);
    end
    out.files{1} = fullfile(path, [nam,'_MAP.mat']);
else
    out.files{1} = fullfile(P.newprobfile.out.dir{1}, P.newprobfile.out.fname);
end
mrstruct_write(probStruct,out.files{1});
