function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_native_calcvals(vatlist, cfile, thresh,obj)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);
% mirroring is not supported yet

Eproj_mirror_enabled = 1;
if Eproj_mirror_enabled == 1
    numPatient = length(obj.allpatients)*2; 
else
    numPatient = length(obj.allpatients);  % no mirroring
end

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

    if side == 1
        side_tag = '_rh';
    else
        side_tag = '_lh';
    end

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient

        disp(['E-field projection ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        if pt <= length(obj.allpatients)
            E_proj_folder = [obj.allpatients{pt},filesep,'miscellaneous',filesep,obj.connectome,filesep,'gs_', obj.M.guid,side_tag];
            if isfile([E_proj_folder,filesep,'E_peak.mat'])
    
                E_proj_Peak = load([E_proj_folder,filesep,'E_peak.mat']);
                E_proj_5Peak = load([E_proj_folder,filesep,'E_5perc_peak.mat']);
    
                % no mirroring allowed atm
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: E-proj doesn''t exist!\n');
                continue;
            end
        else
            % mirrored patient, load from the other side
            if side == 1
                E_proj_folder = [obj.allpatients{pt-size(obj.allpatients,1)},filesep,'miscellaneous',filesep,obj.connectome,filesep,'gs_', obj.M.guid,'_lh'];
            else
                E_proj_folder = [obj.allpatients{pt-size(obj.allpatients,1)},filesep,'miscellaneous',filesep,obj.connectome,filesep,'gs_', obj.M.guid,'_rh'];
            end
                
            if isfile([E_proj_folder,filesep,'E_peak.mat'])
                E_proj_Peak_raw = load([E_proj_folder,filesep,'E_peak.mat']);
                E_proj_5Peak_raw = load([E_proj_folder,filesep,'E_5perc_peak.mat']);
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: E-proj doesn''t exist!\n');
                continue;
            end

            E_proj_Peak.E_peak = zeros(size(E_proj_Peak_raw.E_peak ,1),1);
            E_proj_5Peak.E_5perc_peak = zeros(size(E_proj_5Peak_raw.E_5perc_peak,1),1);


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
                    E_proj_Peak.E_peak(path_start:path_end) = E_proj_Peak_raw.E_peak(path_start_counter:end);
                    E_proj_5Peak.E_5perc_peak(path_start:path_end) = E_proj_5Peak_raw.E_5perc_peak(path_start_counter:end);
                elseif pathway_i == length(obj.map_list)
                    E_proj_Peak.E_peak(path_start:end) = E_proj_Peak_raw.E_peak(path_start_counter:path_end_counter);
                    E_proj_5Peak.E_5perc_peak(path_start:end) = E_proj_5Peak_raw.E_5perc_peak(path_start_counter:path_end_counter);
                else
                    E_proj_Peak.E_peak(path_start:path_end) = E_proj_Peak_raw.E_peak(path_start_counter:path_end_counter);
                    E_proj_5Peak.E_5perc_peak(path_start:path_end) = E_proj_5Peak_raw.E_5perc_peak(path_start_counter:path_end_counter);
                end
                %last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
            end
        end            

        % Find connected fibers, we don't use VATs here (but can create them from magnitude computed in native and warped)
        connected = E_proj_Peak.E_peak*1000.0 > thresh;

        % Generate binary fibsval for the T-test method
        fibsvalBin{side}(connected, pt)=1;

        % Generate fibsval for the Spearman's correlation method
        %fibsvalSum{side}(trimmedFiberInd(connected), pt) = cellfun(@sum, vals);
        %fibsvalMean{side}(trimmedFiberInd(connected), pt) = cellfun(@mean, vals);

        fibsvalPeak{side}(connected, pt) = E_proj_Peak.E_peak(connected)*1000.0;
        fibsval5Peak{side}(connected, pt) = E_proj_5Peak.E_5perc_peak(connected)*1000.0;
    end

    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected = any(fibsvalBin{side}, 2);
    fibsvalBin{side} = sparse(fibsvalBin{side}(fibIsConnected, :));
    fibsvalSum{side} = sparse(fibsvalSum{side}(fibIsConnected, :));
    fibsvalMean{side} = sparse(fibsvalMean{side}(fibIsConnected, :));
    fibsvalPeak{side} = sparse(fibsvalPeak{side}(fibIsConnected, :));
    fibsval5Peak{side} = sparse(fibsval5Peak{side}(fibIsConnected, :));

    % Extract connected fiber cell
    connFiberInd{side} = find(fibIsConnected);
    connFiber = fibers(ismember(fibers(:,4), connFiberInd{side}), 1:3);
    fibcell{side} = mat2cell(connFiber, idx(connFiberInd{side}));
end
