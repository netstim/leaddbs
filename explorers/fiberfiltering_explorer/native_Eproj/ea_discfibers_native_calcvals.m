function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals(vatlist, cfile, thresh,obj)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);
% mirroring is not supported yet
numPatient = size(obj.allpatients,1);

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
        disp(['E-field prjection ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        E_proj_folder = [obj.allpatients{pt},filesep,'miscellaneous',filesep,obj.connectome,filesep,'gs_', obj.M.guid,side_tag];
        if isfile([E_proj_folder,filesep,'E_peak.mat'])

            E_proj_Peak = load([E_proj_folder,filesep,'E_peak.mat']);
            E_proj_5Peak = load([E_proj_folder,filesep,'E_5perc_peak.mat']);

            % no mirroring allowed atm
        else
            ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: VTA doesn''t exist!\n');
            continue;
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
