function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals_pam(pamlist, obj, cfile)
% Extract fiber connection values from OSS-DBS results (for a particular connectome)

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');


numPatient = length(obj.allpatients);  % no mirroring
numSide = 2; % hardcoded for now (as in ...getvats.m)

fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

for side = 1:numSide
    fibsvalBin{side} = zeros(length(idx), numPatient);
    fibsvalSum{side} = zeros(length(idx), numPatient);
    fibsvalMean{side} = zeros(length(idx), numPatient);
    fibsvalPeak{side} = zeros(length(idx), numPatient);
    fibsval5Peak{side} = zeros(length(idx), numPatient);

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient
 
        if obj.multi_pathways == 1 % fiberActivation_side.mat already contains all fibers (incl. filtered out by Kuncel-VTA)
            fib_state_raw = load(char(pamlist(pt,side)));
            total_fibers = length(fib_state_raw.idx);
            fib_state = zeros(total_fibers,1);
            last_loc_i = 1;  
            
            for fib_i = 1:total_fibers
                fib_state(fib_i) = fib_state_raw.fibers(last_loc_i,5);
                last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
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
        activated = find(fib_state == 1);
        %activated = find(fib_state == 1  | fib_state == -1);

        % needed
        % Generate binary fibsval for the T-test method
        fibsvalBin{side}(activated, pt)=1;

    end

    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected = any(fibsvalBin{side}, 2);
    
    fibsvalBin{side} = sparse(fibsvalBin{side}(fibIsConnected, :));
    %fibsvalSum{side} = sparse(fibsvalSum{side}(fibIsConnected, :));
    %fibsvalMean{side} = sparse(fibsvalMean{side}(fibIsConnected, :));
    %fibsvalPeak{side} = sparse(fibsvalPeak{side}(fibIsConnected, :));
    %fibsval5Peak{side} = sparse(fibsval5Peak{side}(fibIsConnected, :));

    % Extract connected fiber cell
    connFiberInd{side} = find(fibIsConnected);
    connFiber = fibers(ismember(fibers(:,4), connFiberInd{side}), 1:3);
    fibcell{side} = mat2cell(connFiber, idx(connFiberInd{side}));
end
