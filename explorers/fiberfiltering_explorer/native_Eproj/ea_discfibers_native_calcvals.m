function [fibsvalBin_proj, fibsvalSum_proj, fibsvalMean_proj, fibsvalPeak_proj, fibsval5Peak_proj, fibcell_proj, connFiberInd_proj,fibsvalBin_magn, fibsvalSum_magn, fibsvalMean_magn, fibsvalPeak_magn, fibsval5Peak_magn, fibcell_magn, connFiberInd_magn, totalFibers] = ea_discfibers_native_calcvals(vatlist, cfile, space, obj)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile);

if strcmp(space,'native')
    space = '';
else
    space = 'MNI/';
end

prefs = ea_prefs;
try
    thresh = obj.calcthreshold;
catch
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[~, numSide] = size(vatlist);

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

fibsvalBin_magn = cell(1, numSide);
fibsvalSum_magn = cell(1, numSide);
fibsvalMean_magn = cell(1, numSide);
fibsvalPeak_magn = cell(1, numSide);
fibsval5Peak_magn = cell(1, numSide);

fibcell_magn = cell(1, numSide);
connFiberInd_magn = cell(1, numSide);

fibsvalBin_proj = cell(1, numSide);
fibsvalSum_proj = cell(1, numSide);
fibsvalMean_proj = cell(1, numSide);
fibsvalPeak_proj = cell(1, numSide);
fibsval5Peak_proj = cell(1, numSide);

fibcell_proj = cell(1, numSide);
connFiberInd_proj = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

for side = 1:numSide
    fibsvalBin_magn{side} = zeros(length(idx), numPatient);
    fibsvalSum_magn{side} = zeros(length(idx), numPatient);
    fibsvalMean_magn{side} = zeros(length(idx), numPatient);
    fibsvalPeak_magn{side} = zeros(length(idx), numPatient);
    fibsval5Peak_magn{side} = zeros(length(idx), numPatient);

    fibsvalBin_proj{side} = zeros(length(idx), numPatient);
    fibsvalSum_proj{side} = zeros(length(idx), numPatient);
    fibsvalMean_proj{side} = zeros(length(idx), numPatient);
    fibsvalPeak_proj{side} = zeros(length(idx), numPatient);
    fibsval5Peak_proj{side} = zeros(length(idx), numPatient);

    if side == 1
        side_tag = '_rh';
    else
        side_tag = '_lh';
    end

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient

        disp(['E-field projection ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        if pt <= length(obj.allpatients)
            E_proj_folder = [obj.allpatients{pt},filesep,'connectomes',filesep,'dMRI',filesep,obj.connectome,filesep,'gs_', obj.M.guid,side_tag];
            if isfile([E_proj_folder,filesep,space,'E_metrics.mat'])
    
                load([E_proj_folder,filesep,space,'E_metrics.mat']);
    
                % no mirroring allowed atm
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: E_metrics file does not exist!\n');
                continue;
            end
        else
            % mirrored patient, load from the other side
            if side == 1
                E_proj_folder = [obj.allpatients{pt-size(obj.allpatients,1)},filesep,'connectomes',filesep,'dMRI',filesep,obj.connectome,filesep,'gs_', obj.M.guid,'_lh'];
            else
                E_proj_folder = [obj.allpatients{pt-size(obj.allpatients,1)},filesep,'connectomes',filesep,'dMRI',filesep,obj.connectome,filesep,'gs_', obj.M.guid,'_rh'];
            end
                
            if isfile([E_proj_folder,filesep,space,'E_metrics.mat'])
                load([E_proj_folder,filesep,space,'E_metrics.mat']);
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: E_metrics file does not exist!\n');
                continue;
            end
            
            % copy results from the other side
            proj_5perc_peak_raw = E_metrics.proj_5perc_peak;
            proj_peak_raw = E_metrics.proj_peak;
            proj_sum_raw = E_metrics.proj_sum;
            proj_mean_raw = E_metrics.proj_mean;

            magn_5perc_peak_raw = E_metrics.magn_5perc_peak;
            magn_peak_raw = E_metrics.magn_peak;
            magn_sum_raw = E_metrics.magn_sum;
            magn_mean_raw = E_metrics.magn_mean;

            % initialize actual mirrored metrics
            E_metrics.proj_peak = zeros(size(proj_peak_raw ,1),1);
            E_metrics.proj_5perc_peak = zeros(size(proj_5perc_peak_raw,1),1);
            E_metrics.proj_sum = zeros(size(proj_sum_raw ,1),1);
            E_metrics.proj_mean = zeros(size(proj_mean_raw,1),1);

            E_metrics.magn_peak = zeros(size(magn_peak_raw ,1),1);
            E_metrics.magn_5perc_peak = zeros(size(magn_5perc_peak_raw,1),1);
            E_metrics.magn_sum = zeros(size(magn_sum_raw ,1),1);
            E_metrics.magn_mean = zeros(size(magn_mean_raw,1),1);
            

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
                    E_metrics.proj_peak(path_start:path_end) = proj_peak_raw(path_start_counter:end);
                    E_metrics.proj_5perc_peak(path_start:path_end) = proj_5perc_peak_raw(path_start_counter:end);
                    E_metrics.proj_sum(path_start:path_end) = proj_sum_raw(path_start_counter:end);
                    E_metrics.proj_mean(path_start:path_end) = proj_mean_raw(path_start_counter:end);

                    E_metrics.magn_peak(path_start:path_end) = magn_peak_raw(path_start_counter:end);
                    E_metrics.magn_5perc_peak(path_start:path_end) = magn_5perc_peak_raw(path_start_counter:end);
                    E_metrics.magn_sum(path_start:path_end) = magn_sum_raw(path_start_counter:end);
                    E_metrics.magn_mean(path_start:path_end) = magn_mean_raw(path_start_counter:end);
                elseif pathway_i == length(obj.map_list)
                    E_metrics.proj_peak(path_start:end) = proj_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.proj_5perc_peak(path_start:end) = proj_5perc_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.proj_sum(path_start:end) = proj_sum_raw(path_start_counter:path_end_counter);
                    E_metrics.proj_mean(path_start:end) = proj_mean_raw(path_start_counter:path_end_counter);

                    E_metrics.magn_peak(path_start:end) = magn_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.magn_5perc_peak(path_start:end) = magn_5perc_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.magn_sum(path_start:end) = magn_sum_raw(path_start_counter:path_end_counter);
                    E_metrics.magn_mean(path_start:end) = magn_mean_raw(path_start_counter:path_end_counter);
                else
                    E_metrics.proj_peak(path_start:path_end) = proj_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.proj_5perc_peak(path_start:path_end) = proj_5perc_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.proj_sum(path_start:path_end) = proj_sum_raw(path_start_counter:path_end_counter);
                    E_metrics.proj_mean(path_start:path_end) = proj_mean_raw(path_start_counter:path_end_counter);

                    E_metrics.magn_peak(path_start:path_end) = magn_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.magn_5perc_peak(path_start:path_end) = magn_5perc_peak_raw(path_start_counter:path_end_counter);
                    E_metrics.magn_sum(path_start:path_end) = magn_sum_raw(path_start_counter:path_end_counter);
                    E_metrics.magn_mean(path_start:path_end) = magn_mean_raw(path_start_counter:path_end_counter);
                end
                %last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
            end
        end            

        % Find connected fibers, we don't use VATs here
        connected_magn = E_metrics.magn_peak*1000.0 > thresh;

        % Generate binary fibsval for the T-test method
        fibsvalBin_magn{side}(connected_magn, pt)=1;
        fibsvalPeak_magn{side}(connected_magn, pt) = E_metrics.magn_peak(connected_magn)*1000.0;
        fibsval5Peak_magn{side}(connected_magn, pt) = E_metrics.magn_5perc_peak(connected_magn)*1000.0;
        fibsvalSum_magn{side}(connected_magn, pt) = E_metrics.magn_sum(connected_magn)*1000.0;
        fibsvalMean_magn{side}(connected_magn, pt) = E_metrics.magn_mean(connected_magn)*1000.0;

        % Find connected fibers, we don't use VATs here
        connected_proj = E_metrics.proj_peak*1000.0 > thresh;

        % Generate binary fibsval for the T-test method
        fibsvalBin_proj{side}(connected_proj, pt)=1;
        fibsvalPeak_proj{side}(connected_proj, pt) = E_metrics.proj_peak(connected_proj)*1000.0;
        fibsval5Peak_proj{side}(connected_proj, pt) = E_metrics.proj_5perc_peak(connected_proj)*1000.0;
        fibsvalSum_proj{side}(connected_proj, pt) = E_metrics.proj_sum(connected_proj)*1000.0;
        fibsvalMean_proj{side}(connected_proj, pt) = E_metrics.proj_mean(connected_proj)*1000.0;
    end

    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected_magn = any(fibsvalBin_magn{side}, 2);
    fibsvalBin_magn{side} = sparse(fibsvalBin_magn{side}(fibIsConnected_magn, :));
    fibsvalSum_magn{side} = sparse(fibsvalSum_magn{side}(fibIsConnected_magn, :));
    fibsvalMean_magn{side} = sparse(fibsvalMean_magn{side}(fibIsConnected_magn, :));
    fibsvalPeak_magn{side} = sparse(fibsvalPeak_magn{side}(fibIsConnected_magn, :));
    fibsval5Peak_magn{side} = sparse(fibsval5Peak_magn{side}(fibIsConnected_magn, :));

    % Extract connected fiber cell
    connFiberInd_magn{side} = find(fibIsConnected_magn);
    connFiber_magn = fibers(ismember(fibers(:,4), connFiberInd_magn{side}), 1:3);
    fibcell_magn{side} = mat2cell(connFiber_magn, idx(connFiberInd_magn{side}));


    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected_proj = any(fibsvalBin_proj{side}, 2);
    fibsvalBin_proj{side} = sparse(fibsvalBin_proj{side}(fibIsConnected_proj, :));
    fibsvalSum_proj{side} = sparse(fibsvalSum_proj{side}(fibIsConnected_proj, :));
    fibsvalMean_proj{side} = sparse(fibsvalMean_proj{side}(fibIsConnected_proj, :));
    fibsvalPeak_proj{side} = sparse(fibsvalPeak_proj{side}(fibIsConnected_proj, :));
    fibsval5Peak_proj{side} = sparse(fibsval5Peak_proj{side}(fibIsConnected_proj, :));

    % Extract connected fiber cell
    connFiberInd_proj{side} = find(fibIsConnected_proj);
    connFiber_proj = fibers(ismember(fibers(:,4), connFiberInd_proj{side}), 1:3);
    fibcell_proj{side} = mat2cell(connFiber_proj, idx(connFiberInd_proj{side}));
end
