function ea_plot_fiber_state(varargin)

    % example 
    % ea_plot_fiber_state('/home/forel/Documents/data/TWEEDMiniset/derivatives/leaddbs/sub-TWEED02/stimulations/MNI152NLin2009bAsym/tst/subthalamopallidal/PAM/sub-TWEED02_sim-fiberActivation_model-ossdbs_hemi-R_tract-gpe2stn_ass_right.mat')

    fiberActivationProb = varargin{1};
    load(fiberActivationProb);

    if nargin >=2
        numfibers = varargin{2};
    else
        numfibers = size(idx,1);
    end

    if numfibers > 500
        downsamplefactor = 5;
    else
        downsamplefactor = 1;
    end

    col = [1,0,0];

    status = zeros(size(idx,1),1);
    jumper = 1;
    % get status at one compartment from each fiber 
    for fiber_i = 1:size(status,1)
        status(fiber_i) = fibers(jumper,5);
        jumper = jumper + idx(fiber_i);
    end

    %% old and slow
%     fibs = unique(fibers(:,4));
%     for k = 1:length(fibs)
%        fibersnew{k,1} =  fibers(fibers(:,4) == fibs(k),1:3);
%     end
    %% new and fast
    fibersnew=mat2cell(fibers(:,1:3),idx);
    %% downsampling
    fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/downsamplefactor))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);

    [maxvals,myfibs] = maxk(status,numfibers);
    
%     %% reduce number of visualized fibers
%     if ~isempty(numfibers) && length(fibersnew) > numfibers
%         myfibs = ceil(linspace(1, length(fibersnew), numfibers));
%     else
%         myfibs = 1:length(fibersnew);
%     end
    %%

    %exp_norm_status = (exp(status)-1.0)/max(exp(status));

    for fiber_i = 1:size(myfibs,1)
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),status(myfibs(fiber_i)));
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),status(myfibs(fiber_i))*0.5);
        if status(myfibs(fiber_i)) == 0
            % not activated
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.05);
            set(mytract,'FaceColor',[1,1,1],'FaceAlpha',0.25,'EdgeColor','none')
        elseif status(myfibs(fiber_i)) == -1 || status(myfibs(fiber_i)) == -3
            % damaged
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.1);
            set(mytract,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',1.0,'EdgeColor','none')
        elseif status(myfibs(fiber_i)) == -2
            % CSF
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.1);
            set(mytract,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',1.0,'EdgeColor','none')
        else
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.1);
            %mytract = streamtube(fibersnew(myfibs(fiber_i)),status(myfibs(fiber_i))*1.0);
            %set(mytract,'FaceColor',[status(myfibs(fiber_i))*0.8+0.2,(1-status(myfibs(fiber_i)))*0.33,0],'FaceAlpha',exp_norm_status(myfibs(fiber_i)),'EdgeColor','none')
            %set(mytract,'FaceColor',[status(myfibs(fiber_i))*0.75+0.25,(1-status(myfibs(fiber_i)))*0.75,0],'FaceAlpha',exp_norm_status(myfibs(fiber_i)),'EdgeColor','none')
            set(mytract,'FaceColor',[1 0 0],'FaceAlpha',1.0,'EdgeColor','none')
        end
    end

    
