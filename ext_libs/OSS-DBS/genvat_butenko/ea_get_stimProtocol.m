function settings = ea_get_stimProtocol(options, S, settings, activeSources, source_index)
% Get stimulation settings for particular source
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options         % Lead-DBS options for electrode reconstruction and stimulation
    S               % Lead-DBS stimulation settings
    settings        % parameters for OSS-DBS simulation
    activeSources   % 2(sides)x4 array of sources' indices. Non active are NaNs
    source_index    {mustBeNumeric} % index of the current source
end

conNum = options.elspec.numContacts;
settings.Activation_threshold_VTA = []; % initialize
nActiveSources = [nnz(~isnan(activeSources(1,:))), nnz(~isnan(activeSources(2,:)))];
if any(nActiveSources > 1)
    settings.multisource = 1;
else
    settings.multisource = 0;
end
settings.stim_center = nan(2, 3);

try
    settings.pulseWidth = [double(S.(['Rs', num2str(source_index)]).pulseWidth);double(S.(['Ls', num2str(source_index)]).pulseWidth)];
catch % Fallback to default pulseWith in case it's not updated in S
    settings.pulseWidth = [options.prefs.machine.vatsettings.butenko_pulseWidth; options.prefs.machine.vatsettings.butenko_pulseWidth];
end

if ~settings.stimSetMode
    [settings.Phi_vector, settings.current_control, settings.Case_grounding] = ea_get_OneSourceStimVector(S, 2, conNum, activeSources(:,source_index));

    for side = 1:2

        % estimate center of VAT grid
        if ~isnan(activeSources(side,source_index))

            % Threshold for Astrom VTA (V/m)
            if settings.adaptive_threshold
                settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;ea_get_adaptiveEthreshold(settings.pulseWidth(side))];
            else
                if options.prefs.machine.vatsettings.butenko_ethresh < 1.0
                    ea_warndlg("The E-threshold is too low (see OSS-DBS settings), likely wrong units are used, upscaling to V/m!")
                    settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;options.prefs.machine.vatsettings.butenko_ethresh*1000.0];
                else
                    settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;options.prefs.machine.vatsettings.butenko_ethresh];
                end
            end

            stimamp = sum(abs(settings.Phi_vector(side,:)),"all",'omitnan');
            settings.stim_center(side,:) = sum(settings.contactLocation{side}.*abs(settings.Phi_vector(side,:)')./stimamp,1,'omitnan');

             % estimate extent of the stimulation along the lead
            phi_temp = settings.Phi_vector(side,:);
            phi_temp(isnan(phi_temp)) = 0;
            first_active = find(phi_temp,1,'first');
            last_active = find(phi_temp,1,'last');
            length_active_span = norm(settings.contactLocation{side}(last_active,:)-settings.contactLocation{side}(first_active,:));
            if strcmp('Boston Scientific Vercise', options.elmodel)
                if length_active_span > 24.0
                    ea_warndlg("Large span of active contacts is detected. Consider extending VAT grid, see PointModel.Lattice.Shape in lead_settings.py")
                end
            elseif length_active_span > 18.0
                ea_warndlg("Large span of active contacts is detected. Consider extending VAT grid, see PointModel.Lattice.Shape in lead_settings.py")
            end
        elseif any(nActiveSources > 1)
            % enforce the settings value for all sources in multisource
            settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;options.prefs.machine.vatsettings.butenko_ethresh];
        else
            % non-active source in single source (no vta will be computed for this side)
            settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;-1.0];
        end

    end
else
    % StimSets
    settings.stim_center(1,:) = mean(settings.contactLocation{1});
    settings.stim_center(2,:) = mean(settings.contactLocation{2});
    settings.Phi_vector = 1000./conNum * ones(2, conNum);
    settings.Case_grounding = 1;
    if settings.exportVAT
        for side = 1:2
            % Threshold for Astrom VTA (V/m)
            if settings.adaptive_threshold
                settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;ea_get_adaptiveEthreshold(settings.pulseWidth(side))];
            else
                % potentially buggy for multisource
                if options.prefs.machine.vatsettings.butenko_ethresh < 1.0
                    ea_warndlg("The E-threshold is too low (see OSS-DBS settings), likely wrong units are used, upscaling to V/m!")
                    settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;options.prefs.machine.vatsettings.butenko_ethresh*1000.0];
                else
                    settings.Activation_threshold_VTA = [settings.Activation_threshold_VTA;options.prefs.machine.vatsettings.butenko_ethresh];
                end
            end
        end
    else
        % not needed for PAM
        settings.Activation_threshold_VTA = [-1;-1];
    end
end

if settings.calcAxonActivation
    settings.connectome = options.prefs.machine.vatsettings.butenko_connectome;
    settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
end

% for multisource PAM, filter fibers using settings.Phi_vector_max
% this simplifies merging results over sources
settings.Phi_vector_max = zeros(size(settings.Phi_vector));
if settings.multisource && settings.calcAxonActivation
    for inx = 1:4

        [Phi_vector, ~, ~] = ea_get_OneSourceStimVector(S, 2, conNum, activeSources(:,inx));
        Phi_vector(isnan(Phi_vector)) = 0;

        for side = 1:size(settings.Phi_vector_max,1)
            if ~isnan(activeSources(side,inx))
                for cnt = 1:size(settings.Phi_vector_max,2)
                    if abs(Phi_vector(side,cnt)) > settings.Phi_vector_max(side,cnt)
                        settings.Phi_vector_max(side,cnt) = abs(Phi_vector(side,cnt));
                    end
                end
            end
        end
    end
end