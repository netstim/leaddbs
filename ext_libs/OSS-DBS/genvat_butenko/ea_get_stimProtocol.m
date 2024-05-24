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

conNum = options.elspec.numel;
settings.Activation_threshold_VTA = []; % initialize
nActiveSources = [nnz(~isnan(activeSources(1,:))), nnz(~isnan(activeSources(2,:)))];
settings.stim_center = nan(2, 3);

if ~settings.stimSetMode
    [settings.Phi_vector, settings.current_control, settings.Case_grounding] = ea_get_OneSourceStimVector(S, 2, conNum, activeSources(:,source_index));
    settings.pulseWidth = [double(S.(['Rs', num2str(source_index)]).pulseWidth);double(S.(['Ls', num2str(source_index)]).pulseWidth)];
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
end
%else
%    settings.stim_center = [NaN;NaN];
%end

if settings.calcAxonActivation
    %settings.pulseWidth = [S.Rs1.pulseWidth;S.Ls1.pulseWidth];
    settings.connectome = options.prefs.machine.vatsettings.butenko_connectome;
    settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
end