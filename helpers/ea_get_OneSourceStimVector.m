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

amp = nan(eleNum,1);
for i = 1:eleNum

    % should be subsituted to an automatic solution
    % but keep in mind the skip
    if conNum <= 8
        switch i
            case 1
                sideCode = 'R';
                cntlabel = {'k0','k1','k2','k3','k4','k5','k6','k7'};
            case 2
                sideCode = 'L';
                cntlabel = {'k8','k9','k10','k11','k12','k13','k14','k15'};
        end
    else
        switch i
            case 1
                sideCode = 'R';
                cntlabel = {'k0','k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12','k13','k14','k15'};
            case 2
                sideCode = 'L';
                cntlabel = {'k16','k17','k18','k19','k20','k21','k22','k23','k24','k25','k26','k27','k28','k29','k30','k31'};
        end
    end

    if ~isnan(source_i(i))

        % get stimulation mode (0 - Voltage Control; 1 - Current Control) and amplitude
        current_control(i) = uint8(S.([sideCode, 's', num2str(source_i(i))]).va==2);
        amp(i) = S.amplitude{i}(source_i(i));
        stimSource = S.([sideCode, 's', num2str(source_i(i))]);

        % VC case pre-set
        if stimSource.va == 1

            % Split voltage in case contacts have both polarities
            v_source = S.([sideCode, 's', num2str(source_i(i))]);
            if v_source.case.pol == 0
                amp(i) = amp(i)/2;
            end
        end

        for cnt = 1:conNum
            v_source = S.([sideCode, 's', num2str(source_i(i))]);
            if current_control(i) == 1
                % OSS-DBS will make sure that all excess of current is
                % grounded
                if S.activecontacts{i}(cnt)
                    switch stimSource.(cntlabel{cnt}).pol
                        case 1 % Negative, cathode
                            Phi_vector(i, cnt) = -amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                        case 2 % Positive, anode
                            Phi_vector(i, cnt) = amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                    end
                end
            else
                % IMPORTANT: OSS-DBS will shift everything to
                % negative volts and grounding
                % if you have bipolar stim 3V k0+ k1-, OSS-DBS will solve
                % for 0V vs -3V
                if v_source.(cntlabel{cnt}).perc    % the contact is active for this source
                    switch v_source.(cntlabel{cnt}).pol
                        case 1
                            Phi_vector(i, cnt) = -1*amp(i);
                        case 2
                            Phi_vector(i, cnt) =  amp(i);
                    end
                end
            end
        end
        % check grounding for the active source
        if amp(i) ~= 0 && ~isnan(amp(i))
            if v_source.case.perc == 100
                Case_grounding(i) = 1;
            end
        end  
    end
end
