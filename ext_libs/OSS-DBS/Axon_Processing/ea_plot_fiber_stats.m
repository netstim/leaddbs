function ea_plot_fiber_stats(varargin)

    % example 
    % ea_plot_fiber_stats('/home/forel/Documents/data/JS/python_ff/Fiber_stats.mat')

    fiberStats = varargin{1};
    load(fiberStats);

    idx = idx';  % can't export as column from python

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

    stat = zeros(size(idx,1),1);
    jumper = 1;
    % get status at one compartment from each fiber 
    for fiber_i = 1:size(stat,1)
        stat(fiber_i) = fibers(jumper,5);
        jumper = jumper + idx(fiber_i);
    end

    %% new and fast
    fibersnew=mat2cell(fibers(:,1:3),idx);
    %% downsampling
    fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/downsamplefactor))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);

    [maxvals,myfibs] = maxk(stat,numfibers);
    
%     %% reduce number of visualized fibers
%     if ~isempty(numfibers) && length(fibersnew) > numfibers
%         myfibs = ceil(linspace(1, length(fibersnew), numfibers));
%     else
%         myfibs = 1:length(fibersnew);
%     end
    %%


    stat_abs = abs(stat);
    exp_norm_stat = (exp(stat_abs)-1.0)/max(exp(stat_abs));

    for fiber_i = 1:size(myfibs,1)
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),stat(myfibs(fiber_i)));
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),stat(myfibs(fiber_i))*0.5);
        if stat(myfibs(fiber_i)) == 0
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.05);
            set(mytract,'FaceColor',[1,1,1],'FaceAlpha',0.25,'EdgeColor','none')
        elseif stat(myfibs(fiber_i)) > 0.0
            mytract = streamtube(fibersnew(myfibs(fiber_i)),exp_norm_stat(myfibs(fiber_i))*0.3+0.1);
            set(mytract,'FaceColor',[1.0,(1-stat_abs(myfibs(fiber_i))),0],'FaceAlpha',exp_norm_stat(myfibs(fiber_i)),'EdgeColor','none')
        else
            mytract = streamtube(fibersnew(myfibs(fiber_i)),exp_norm_stat(myfibs(fiber_i))*0.3+0.1);
            set(mytract,'FaceColor',[0.0,(1-stat_abs(myfibs(fiber_i))),1.0],'FaceAlpha',exp_norm_stat(myfibs(fiber_i)),'EdgeColor','none')
        end
    end

    
