function ea_plot_prob_fiber_state(varargin)

    % example 
    % ea_plot_prob_fiber_state('/media/konstantin/Konstantin/StimFit_Cohort/StimFitBIDS/derivatives/leaddbs/sub-SBN7F4E4/stimulations/native/gs_20230506011354/sub-SBN7F4E4_sim-fiberActivation_model-ossdbs_hemi-L_tract-cerebellothalamic_left.mat', 125)

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

    probability = zeros(size(idx,1),1);
    jumper = 1;
    % get status at one compartment from each fiber 
    for fiber_i = 1:size(probability,1)
        probability(fiber_i) = fibers(jumper,5);
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

    [maxvals,myfibs] = maxk(probability,numfibers);
    
%     %% reduce number of visualized fibers
%     if ~isempty(numfibers) && length(fibersnew) > numfibers
%         myfibs = ceil(linspace(1, length(fibersnew), numfibers));
%     else
%         myfibs = 1:length(fibersnew);
%     end
    %%

    exp_norm_probability = (exp(probability)-1.0)/max(exp(probability));

    for fiber_i = 1:size(myfibs,1)
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),probability(myfibs(fiber_i)));
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),probability(myfibs(fiber_i))*0.5);
        if probability(myfibs(fiber_i)) == 0
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.05);
            set(mytract,'FaceColor',[1,1,1],'FaceAlpha',0.25,'EdgeColor','none')
        else
            mytract = streamtube(fibersnew(myfibs(fiber_i)),exp_norm_probability(myfibs(fiber_i))*0.3+0.1);
            %mytract = streamtube(fibersnew(myfibs(fiber_i)),probability(myfibs(fiber_i))*1.0);
            %set(mytract,'FaceColor',[probability(myfibs(fiber_i))*0.8+0.2,(1-probability(myfibs(fiber_i)))*0.33,0],'FaceAlpha',exp_norm_probability(myfibs(fiber_i)),'EdgeColor','none')
            %set(mytract,'FaceColor',[probability(myfibs(fiber_i))*0.75+0.25,(1-probability(myfibs(fiber_i)))*0.75,0],'FaceAlpha',exp_norm_probability(myfibs(fiber_i)),'EdgeColor','none')
            set(mytract,'FaceColor',[1.0,(1-probability(myfibs(fiber_i))),0],'FaceAlpha',exp_norm_probability(myfibs(fiber_i)),'EdgeColor','none')
        end
    end

    
