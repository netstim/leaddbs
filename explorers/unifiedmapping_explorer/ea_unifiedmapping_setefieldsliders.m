function [statsettings, sliderState] = ...
    ea_unifiedmapping_setefieldsliders(statsettings, calcsettings)
%EA_UNIFIEDMAPPING_SETEFIELDSLIDERS
% Determine the limits, values, and labels for the unified-mapping
% threshold sliders.
%
% This function contains no App Designer dependencies. The .mlapp applies
% the returned sliderState values to the actual UI components.

%% Network mapping
% Network thresholds are stored as percentages:
% 30 corresponds to an absolute correlation threshold of |r| >= 0.30.

networkThreshold = getSetting( ...
    statsettings, 'efieldthreshold_network', 50);

networkThreshold = clampValue(networkThreshold, [0 100]);

statsettings.efieldthreshold_network = networkThreshold;

sliderState.network.used = true;
sliderState.network.limits = [0 100];
sliderState.network.value = networkThreshold;
sliderState.network.label = { ...
    'Voxels are retained when absolute connectivity', ...
    ['|r| is >= ', sprintf('%.2f', networkThreshold / 100)]};

%% Sweet spot and fiber filtering
switch statsettings.stimulationmodel

    case 'Electric Field'

        % Sweet spot
        spotLowerLimit = calcsettings.calcthreshold;
        spotUpperLimit = max(500, spotLowerLimit);
        spotLimits = [spotLowerLimit spotUpperLimit];

        spotThreshold = getSetting( ...
            statsettings, 'efieldthreshold_spot', 200);

        % Convert a previously stored sigmoid probability back to V/m.
        if spotThreshold > 0 && spotThreshold < 1
            spotThreshold = ea_EfieldFromSigmoid(spotThreshold);
        end

        spotThreshold = clampValue(spotThreshold, spotLimits);

        statsettings.efieldthreshold_spot = spotThreshold;

        sliderState.spot.used = true;
        sliderState.spot.limits = spotLimits;
        sliderState.spot.value = spotThreshold;
        sliderState.spot.label = { ...
            'Voxels are connected when E-field magnitude', ...
            ['is > ', sprintf('%.1f', spotThreshold), ' V/m']};

        %% Fiber filtering
        switch lower(statsettings.efieldmetric)
            case 'sum'
                tractLimits = [500 5000];

            case 'mean'
                tractLimits = [2.5 1000];

            case {'peak', 'peak 5%'}
                tractLimits = [50 5000];

            otherwise
                tractLimits = [50 5000];
        end

        tractThreshold = getSetting( ...
            statsettings, 'efieldthreshold_tract', 200);

        % Convert a previously stored sigmoid probability back to V/m.
        if tractThreshold > 0 && tractThreshold < 1
            tractThreshold = ea_EfieldFromSigmoid(tractThreshold);
        end

        tractThreshold = clampValue(tractThreshold, tractLimits);

        statsettings.efieldthreshold_tract = tractThreshold;

        sliderState.tract.used = true;
        sliderState.tract.limits = tractLimits;
        sliderState.tract.value = tractThreshold;

        metricName = lower(statsettings.efieldmetric);

        sliderState.tract.label = { ...
            ['Tracts are connected when ', metricName, ...
             ' E-field magnitude'], ...
            ['is > ', sprintf('%.1f', tractThreshold), ' V/m']};

    case 'Sigmoid Field'

        probabilityLowerLimit = round( ...
            ea_SigmoidFromEfield(calcsettings.calcthreshold), 2);

        % Ensure that the lower limit remains strictly below 1.
        probabilityLowerLimit = min( ...
            max(probabilityLowerLimit, 0), 0.99);

        probabilityLimits = [probabilityLowerLimit 1];

        %% Sweet spot
        spotThreshold = getSetting( ...
            statsettings, 'efieldthreshold_spot', 200);

        if spotThreshold > 1
            spotThreshold = ea_SigmoidFromEfield(spotThreshold);
        end

        spotThreshold = clampValue( ...
            spotThreshold, probabilityLimits);

        statsettings.efieldthreshold_spot = spotThreshold;

        sliderState.spot.used = true;
        sliderState.spot.limits = probabilityLimits;
        sliderState.spot.value = spotThreshold;
        sliderState.spot.label = { ...
            'Voxels are connected when activation probability', ...
            ['is > ', sprintf('%.2f', spotThreshold)]};

        %% Fiber filtering
        tractThreshold = getSetting( ...
            statsettings, 'efieldthreshold_tract', 200);

        if tractThreshold > 1
            tractThreshold = ea_SigmoidFromEfield(tractThreshold);
        end

        tractThreshold = clampValue( ...
            tractThreshold, probabilityLimits);

        statsettings.efieldthreshold_tract = tractThreshold;

        sliderState.tract.used = true;
        sliderState.tract.limits = probabilityLimits;
        sliderState.tract.value = tractThreshold;
        sliderState.tract.label = { ...
            'Tracts are connected when activation probability', ...
            ['is > ', sprintf('%.2f', tractThreshold)]};

    case 'VTA'

        % E-field thresholds are not used for VTA analyses. Preserve the
        % stored spot and tract values so that changing back to an E-field
        % model restores them.
        sliderState.spot.used = false;
        sliderState.spot.limits = [];
        sliderState.spot.value = [];
        sliderState.spot.label = ...
            'E-field threshold is not used for VTAs';

        sliderState.tract.used = false;
        sliderState.tract.limits = [];
        sliderState.tract.value = [];
        sliderState.tract.label = ...
            'E-field threshold is not used for VTAs';

    otherwise
        error( ...
            'Unknown stimulation model: %s', ...
            statsettings.stimulationmodel);
end

end


function value = getSetting(settings, fieldName, defaultValue)
% Return a valid scalar setting or its default.

if ~isfield(settings, fieldName) || ...
        isempty(settings.(fieldName)) || ...
        ~isnumeric(settings.(fieldName)) || ...
        ~isscalar(settings.(fieldName)) || ...
        ~isfinite(settings.(fieldName))

    value = defaultValue;
else
    value = settings.(fieldName);
end

end


function value = clampValue(value, limits)
% Restrict a scalar value to the supplied slider limits.

value = min(max(value, limits(1)), limits(2));

end