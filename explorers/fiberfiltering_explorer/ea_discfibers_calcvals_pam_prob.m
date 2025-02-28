function [fibsvalBin, fibsvalProb, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals_pam_prob(pamlist, obj, cfile)
% Extract fiber connection values from OSS-DBS results (for a particular connectome)
% use "probabilistic" PAM results

disp('Load Connectome...');
load(cfile);

% allow mirroring only for mirrored connectomes
if exist('mirrored','var')
    if mirrored && obj.multi_pathways
        numPatient = length(obj.allpatients)*2; 
    else
        numPatient = length(obj.allpatients);  % no mirroring
    end
else
    numPatient = length(obj.allpatients);  % no mirroring
end

numSide = 2; % hardcoded for now (as in ...getvats.m)

fibsvalBin = cell(1, numSide);
fibsvalProb = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

for side = 1:numSide
    fibsvalBin{side} = zeros(length(idx), numPatient);
    fibsvalProb{side} = zeros(length(idx), numPatient);
    fibsvalMean{side} = zeros(length(idx), numPatient);
    fibsvalPeak{side} = zeros(length(idx), numPatient);
    fibsval5Peak{side} = zeros(length(idx), numPatient);

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient
 
        if obj.multi_pathways == 1 % fiberActivation_side.mat already contains all fibers (incl. filtered out by Kuncel-VTA)

            % for mirrored patients, we will find a counterpart fibers in
            % another hemisphere

            if strcmp(char(pamlist(pt,side)), 'skip')
                % no stimulation for this hemisphere
                continue
            else
                try
                    fib_state_raw = load(char(pamlist(pt,side)));
                catch
                    warning("No merged fiberActivation file for this patient")
                    continue
                end
            end

            total_fibers = length(fib_state_raw.idx);
            fib_state = zeros(total_fibers,1);
            last_loc_i = 1;  

            if pt <= length(obj.allpatients)

                % original activations across all pathways and fibers               
                for fib_i = 1:total_fibers
                    fib_state(fib_i) = fib_state_raw.fibers(last_loc_i,5);
                    last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
                end
            else

                % excessive step, but simplifies logic
                fib_state_non_mirror = zeros(total_fibers,1);
                for fib_i = 1:total_fibers
                     fib_state_non_mirror(fib_i) = fib_state_raw.fibers(last_loc_i,5);
                    last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
                end

                % for mirrored we find indices of pathway counterparts as defined in
                % obj.map_list (order is path1_rh,path1_lh,path2_rh...)
                for pathway_i = 1:length(obj.map_list)
                    path_start = obj.map_list(pathway_i);

                    if pathway_i ~= length(obj.map_list)
                        path_end = obj.map_list(pathway_i+1) - 1;
                    end

                    if rem(pathway_i,2)
                        path_start_counter = obj.map_list(pathway_i+1);
                        if pathway_i == length(obj.map_list)-1
                            disp("prelast pathway")
                        else
                            path_end_counter = obj.map_list(pathway_i+2) - 1;
                        end
                    else
                        path_start_counter = obj.map_list(pathway_i-1);
                        path_end_counter = obj.map_list(pathway_i) - 1;                        
                    end

                    % copy fiber state to the counterpart
                    if pathway_i == length(obj.map_list)-1

                        fib_state(path_start:path_end) = fib_state_non_mirror(path_start_counter:end);
                    elseif pathway_i == length(obj.map_list)
                        fib_state(path_start:end) = fib_state_non_mirror(path_start_counter:path_end_counter);
                    else
                        fib_state(path_start:path_end) = fib_state_non_mirror(path_start_counter:path_end_counter);
                    end
                    %last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
                end
            end
        else

            load(cfile, 'fibers', 'idx');
            total_fibers = fibers(end,4); % load the actual .mat
            fib_state = zeros(total_fibers,1);

            try
                fib_state_raw = load(char(pamlist(pt,side)));
            catch
                disp("=================== WARNING ==================")
                disp("fiberActivation was not found for this patient")
                disp("perhaps Kuncel-VTA removed all fibers")
                disp("assigning zero activation")
                disp("==============================================")
                continue
            end
                
            %if ~strcmp(obj.connectome, fib_state_raw.connectome_name)
            %    error("=== Fiber activation was computed for another connectome!!! ===") 
            %end    

            last_loc_i = 1;  
            sub_i = 1;
            for fib_i = 1:total_fibers
                if fib_i > fib_state_raw.fibers(end,4)
                    fib_state(fib_i) = 0;  % the fiber was pre-filtered out with Kuncel-VTA
                else
                    if fib_state_raw.fibers(last_loc_i,4) == fib_i
                        fib_state(fib_i) = fib_state_raw.fibers(last_loc_i,5);
                        last_loc_i = fib_state_raw.idx(sub_i)+last_loc_i;
                        sub_i = sub_i + 1;    
                    else
                        fib_state(fib_i) = 0;  % the fiber was pre-filtered out with Kuncel-VTA
                    end
                end
            end
        end

        
        %maybe this is a wrong way
        % Skip further calculation in case no fibers were activated
        if ~any(fib_state)
            continue;
        end
       
        % alternatively, you could also add fib_state == -1
        % probabilistic_PA
        activated = find(fib_state >= 0.05);     % use low threshold when doing pPAM
        %activated = find(fib_state >= 0.5);    % maybe use higher threshold when doing binary tests

        % needed
        % Generate binary fibsval for the T-test method
        fibsvalBin{side}(activated, pt)=1;
        fibsvalProb{side}(activated, pt)=fib_state(activated);

    end

    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected = any(fibsvalBin{side}, 2);
    
    fibsvalBin{side} = sparse(fibsvalBin{side}(fibIsConnected, :));
    fibsvalProb{side} = sparse(fibsvalProb{side}(fibIsConnected, :));
    %fibsvalMean{side} = sparse(fibsvalMean{side}(fibIsConnected, :));
    %fibsvalPeak{side} = sparse(fibsvalPeak{side}(fibIsConnected, :));
    %fibsval5Peak{side} = sparse(fibsval5Peak{side}(fibIsConnected, :));

    % Extract connected fiber cell
    connFiberInd{side} = find(fibIsConnected);
    connFiber = fibers(ismember(fibers(:,4), connFiberInd{side}), 1:3);
    fibcell{side} = mat2cell(connFiber, idx(connFiberInd{side}));
end
