function [Phi_vector, current_control, Case_grounding] = ea_getStimVector(S, eleNum, conNum)
% Get stimulation vectors (for OSS-DBS) based on S.
% eleNum: number of electrodes
% conNum: number of contacts per electrode

% Initialize current control flag
current_control = nan(eleNum, 1);

% Initalize Phi vector
Phi_vector = nan(eleNum, conNum);

% Initialize grounding
Case_grounding = zeros(eleNum,1);

% Get the stimulation parameters from S in case stimSetMode is 0, otherwise
% they will be loaded directly from the Current_protocols_[0|1].csv files
source = nan(eleNum,1);
for i=1:eleNum
    if any(S.amplitude{i})
        source(i) = find(S.amplitude{i},1);
    end
end

% 0 - Voltage Control; 1 - Current Control
for i=1:eleNum
    switch i
        case 1
            sideCode = 'R';
        case 2
            sideCode = 'L';
    end

    if ~isnan(source(i))
        current_control(i) = uint8(S.([sideCode, 's', num2str(source(i))]).va==2);
    end
end

% Get stimulation amplitude
amp = nan(eleNum,1);
for i=1:eleNum
    if ~isnan(source(i))
        amp(i) = S.amplitude{i}(source(i));
    end
end

%  For VC, check all sources
if (current_control(1) == 0 && current_control(2) == 0) || (isnan(current_control(1)) && current_control(2) == 0) || (current_control(1) == 0 && isnan(current_control(2)))% both hemisphere MUST have the same mode
    numSources = 4;
    amp = nan(eleNum,numSources);   % 4 - number of sources
    for i=1:eleNum
        for j=1:numSources
            amp(i,j) = S.amplitude{i}(j);
        end
    end
else
    amp = nan(eleNum,1);  % For CC, one source is used, add a check
    for i=1:eleNum
        if ~isnan(source(i))
            amp(i) = S.amplitude{i}(source(i));
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

    if ~isnan(source(i))
        stimSource = S.([sideCode, 's', num2str(source(i))]);

        % Split voltage in case contacts have both polarities
        if stimSource.va == 1
            for j=1:numSources
                v_source = S.([sideCode, 's', num2str(j)]);
                if v_source.case.pol == 0
                    amp(i,j) = amp(i,j)/2;
                end
            end
        end

        for cnt = 1:conNum
            if current_control(i) == 1  % only one source for CC
                if S.activecontacts{i}(cnt)
                    switch stimSource.(cntlabel{cnt}).pol
                        case 1 % Negative, cathode
                            Phi_vector(i, cnt) = -amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                        case 2 % Postive, anode
                            Phi_vector(i, cnt) = amp(i)*stimSource.(cntlabel{cnt}).perc/100;
                    end
                end
            else
                for j=1:numSources  % go over sources
                    v_source = S.([sideCode, 's', num2str(j)]);
                    if v_source.(cntlabel{cnt}).perc    % the contact is active for this source
                        if isnan(Phi_vector(i, cnt))
                            Phi_vector(i, cnt) = 0.0;  % initialize
                        end
                        switch v_source.(cntlabel{cnt}).pol  % sanity check needed: same polarity for a contact over all sources
                            case 1
                                Phi_vector(i, cnt) = Phi_vector(i, cnt) - amp(i,j);
                            case 2
                                Phi_vector(i, cnt) = Phi_vector(i, cnt) + amp(i,j);
                        end
                    end
                end
            end
        end
        if stimSource.case.perc == 100
            Case_grounding(i) = 1;
        end
    end
end
