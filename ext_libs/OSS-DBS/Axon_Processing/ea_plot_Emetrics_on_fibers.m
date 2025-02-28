function ea_plot_Emetrics_on_fibers(varargin)
    % example 
    % ea_plot_Emetrics_on_fibers('/media/konstantin/Konstantin/StimFit_Cohort/StimFitBIDS/derivatives/leaddbs/sub-SBN7F4E4/miscellaneous/LentassCerb/merged_pathways.mat','/media/konstantin/Konstantin/StimFit_Cohort/StimFitBIDS/derivatives/leaddbs/sub-SBN7F4E4/miscellaneous/LentassCerb/gs_20230506011354_lh/E_metrics.mat',500)

    connectome = varargin{1};
    load(connectome);

    Epeak_file = varargin{2};
    load(Epeak_file);

    if nargin >=3
        numfibers = varargin{3};
    else
        numfibers = size(idx,1);
    end

    if numfibers > 500
        downsamplefactor = 5;
    else
        downsamplefactor = 1;
    end

    col = [1,0,0];

    % probability = zeros(size(idx,1),1);
    % jumper = 1;
    % % get status at one compartment from each fiber 
    % for fiber_i = 1:size(probability,1)
    %     probability(fiber_i) = fibers(jumper,5);
    %     jumper = jumper + idx(fiber_i);
    % end

    %% old and slow
%     fibs = unique(fibers(:,4));
%     for k = 1:length(fibs)
%        fibersnew{k,1} =  fibers(fibers(:,4) == fibs(k),1:3);
%     end
    %% new and fast
    fibersnew=mat2cell(fibers(:,1:3),idx);
    %% downsampling
    fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/downsamplefactor))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);

    [maxvals,myfibs] = maxk(E_metrics.magn_peak,numfibers);
    
%     %% reduce number of visualized fibers
%     if ~isempty(numfibers) && length(fibersnew) > numfibers
%         myfibs = ceil(linspace(1, length(fibersnew), numfibers));
%     else
%         myfibs = 1:length(fibersnew);
%     end
    %%

    exp_norm_Epeak = (exp(E_metrics.magn_peak)-1.0)/max(exp(E_metrics.magn_peak));

    norm_Epeak = (E_metrics.magn_peak)/max(E_metrics.magn_peak);

    total_vis = 0;
    for fiber_i = 1:size(myfibs,1)
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),probability(myfibs(fiber_i)));
        %mytract = streamtube(fibersnew(myfibs(fiber_i)),probability(myfibs(fiber_i))*0.5);
        if E_metrics.magn_peak(myfibs(fiber_i)) == 0
            %continue
            mytract = streamtube(fibersnew(myfibs(fiber_i)),0.05);
            set(mytract,'FaceColor',[1,1,1],'FaceAlpha',0.25,'EdgeColor','none')
        else
            mytract = streamtube(fibersnew(myfibs(fiber_i)),norm_Epeak(myfibs(fiber_i))*0.3+0.1);
            %mytract = streamtube(fibersnew(myfibs(fiber_i)),probability(myfibs(fiber_i))*1.0);
            %set(mytract,'FaceColor',[probability(myfibs(fiber_i))*0.8+0.2,(1-probability(myfibs(fiber_i)))*0.33,0],'FaceAlpha',exp_norm_probability(myfibs(fiber_i)),'EdgeColor','none')


            %set(mytract,'FaceColor',[1-norm_Epeak(myfibs(fiber_i)),1,0],'FaceAlpha',norm_Epeak(myfibs(fiber_i)),'EdgeColor','none')   % yellow-green
            set(mytract,'FaceColor',[0,1-norm_Epeak(myfibs(fiber_i)),1],'FaceAlpha',norm_Epeak(myfibs(fiber_i)),'EdgeColor','none')   
            
            total_vis = total_vis + 1;
        end
    end

    disp(total_vis)

    
