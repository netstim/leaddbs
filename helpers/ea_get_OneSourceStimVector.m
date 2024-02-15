function [Phi_vector, current_control, Case_grounding] = ea_get_OneSourceStimVector(S, eleNum, conNum,source_i)
% Get stimulation vectors (for OSS-DBS) based on S.
% eleNum: number of electrodes
% conNum: number of contacts per electrode
% source_i = eleNum,1. Index of the source, NaN entry for no stimulation

% Initialize current control flag
current_control = nan(eleNum, 1);

% Initalize Phi vector
Phi_vector = nan(eleNum, conNum);

% Initialize grounding
Case_grounding = zeros(eleNum,1);

% Get the stimulation parameters from S in case stimSetMode is 0, otherwise
% they will be loaded directly from the Current_protocols_[0|1].csv files

% 0 - Voltage Control; 1 - Current Control
for i=1:eleNum
    switch i
        case 1
            sideCode = 'R';
        case 2
            sideCode = 'L';
    end

    if ~isnan(source_i(i))
        current_control(i) = uint8(S.([sideCode, 's', num2str(source_i(i))]).va==2);
    end
end

% Get stimulation amplitude
amp = nan(eleNum,1);
for i=1:eleNum
    if ~isnan(source_i(i))
        amp(i) = S.amplitude{i}(source_i(i));
    end
end

%  For VC, check all sources
if (current_control(1) == 0 && current_control(2) == 0) || (isnan(current_control(1)) && current_control(2) == 0) || (current_control(1) == 0 && isnan(current_control(2)))% both hemisphere MUST have the same mode
    %amp = nan(eleNum,1);   
    for i=1:eleNum
        switch i
            case 1
                sideCode = 'R';
            case 2
                sideCode = 'L';
        end

        %amp(i) = S.amplitude{i}(source_i(i));

        % check grounding for the active source
        if amp(i) ~= 0 && ~isnan(amp(i))
            stimSource = S.([sideCode, 's', num2str(source_i(i))]);
            if stimSource.case.perc == 100
                Case_grounding(i) = 1;
            end
        end  

    end
end


for i = 1:eleNum
    switch i
        case 1
            sideCode = 'R';
            cntlabel = {'k0','k1','k2','k3','k4','k5','k6','k7'};
        case 2
            sideCode = 'L';
            cntlabel = {'k8','k9','k10','k11','k12','k13','k14','k15'};
    end

    if ~isnan(source_i(i))
        stimSource = S.([sideCode, 's', num2str(source_i(i))]);

        % Split voltage in case contacts have both polarities
        if stimSource.va == 1

            v_source = S.([sideCode, 's', num2str(source_i(i))]);
            if v_source.case.pol == 0
                amp(i) = amp(i)/2;
            end

        end

        for cnt = 1:conNum
            if current_control(i) == 1
                if S.activecontacts{i}(cnt)
                    switch stimSource.(cntlabel{cnt}).pol
                        case 1 % Negative, cathode
                            Phi_vector(i, cnt) = -amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                        case 2 % Postive, anode
                            Phi_vector(i, cnt) = amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                    end
                end
            else

                v_source = S.([sideCode, 's', num2str(source_i(i))]);
                if v_source.(cntlabel{cnt}).perc    % the contact is active for this source
%                     if isnan(Phi_vector(i, cnt))
%                         Phi_vector(i, cnt) = 0.0;  % initialize
%                     end
                    switch v_source.(cntlabel{cnt}).pol
                        case 1
                            Phi_vector(i, cnt) = -1*amp(i);
                        case 2
                            Phi_vector(i, cnt) =  amp(i);
                    end
                end

            end
        end
        if stimSource.case.perc == 100
            Case_grounding(i) = 1;
        end
    end
end
