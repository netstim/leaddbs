function [statsettings, sliderState] = ea_unifiedmapping_setefieldsliders(explorer)
%EA_UNIFIEDMAPPING_SETEFIELDSLIDERS Determine the limits, values, scales,
% and labels for the unified-mapping threshold sliders.
%
% Thresholds in statsettings always remain in the raw data units. The
% returned slider values and limits are divided by sliderState.*.scale for
% compact display in App Designer.

statsettings = explorer.statsettings;
calcsettings = explorer.calcsettings;
isPseudoM = isprop(explorer, 'M') && ...
    isstruct(explorer.M) && isfield(explorer.M, 'pseudoM');

%% Network mapping
% Network thresholds are percentages: 30 corresponds to |r| >= 0.30.
networkThreshold = getSetting( ...
    statsettings, 'efieldthreshold_network', 50);
networkThreshold = clampValue(networkThreshold, [0 100]);
statsettings.efieldthreshold_network = networkThreshold;

sliderState.network.used = true;
sliderState.network.limits = [0 100];
sliderState.network.value = networkThreshold;
sliderState.network.scale = 1;
sliderState.network.label = { ...
    'Voxels are retained when absolute connectivity', ...
    ['|r| is >= ', sprintf('%.2f', networkThreshold / 100)]};

%% Sweet spot and fiber filtering
switch statsettings.stimulationmodel

    case 'Electric Field'
        [quantityName, unitName] = getFieldDescription(isPseudoM);

        %% Sweet spot
        spotLowerLimit = calcsettings.calcthreshold;
        spotUpperLimit = max(500, spotLowerLimit);
        hasSpotResults = false;

        if isfield(explorer.results, 'sweetspotmapping') && ...
                isfield(explorer.results.sweetspotmapping, 'efield') && ...
                ~isempty(explorer.results.sweetspotmapping.efield)
            [calculatedMaximum, hasSpotResults] = getCellMaximum( ...
                explorer.results.sweetspotmapping.efield);
            if hasSpotResults
                spotUpperLimit = max(calculatedMaximum, spotLowerLimit);
            end
        end

        if hasSpotResults
            [spotScale, spotExponent] = ...
                getEngineeringScale(spotUpperLimit);
            spotUpperLimit = makeNiceUpperLimit( ...
                spotUpperLimit, spotScale);
        else
            spotScale = 1;
            spotExponent = 0;
        end

        spotUpperLimit = ensureValidUpperLimit( ...
            spotLowerLimit, spotUpperLimit, spotScale);
        spotLimitsRaw = [spotLowerLimit, spotUpperLimit];
        spotLimitsDisplayed = spotLimitsRaw / spotScale;

        spotThreshold = getSetting( ...
            statsettings, 'efieldthreshold_spot', 200);

        % Convert a previously stored sigmoid probability back to raw units.
        if spotThreshold > 0 && spotThreshold < 1
            spotThreshold = ea_EfieldFromSigmoid(spotThreshold);
        end

        spotThreshold = clampValue(spotThreshold, spotLimitsRaw);
        statsettings.efieldthreshold_spot = spotThreshold;

        sliderState.spot.used = true;
        sliderState.spot.limits = spotLimitsDisplayed;
        sliderState.spot.value = spotThreshold / spotScale;
        sliderState.spot.scale = spotScale;
        sliderState.spot.exponent = spotExponent;
        sliderState.spot.quantityName = quantityName;
        sliderState.spot.unitLabel = ...
            makeUnitLabel(spotExponent, unitName);
        sliderState.spot.label = { ...
            ['Voxels are connected when ', quantityName], ...
            ['is > ', sprintf('%.1f', spotThreshold / spotScale), ...
             ' ', sliderState.spot.unitLabel]};

        %% Fiber filtering
        tractLimitsRaw = getFallbackTractLimits( ...
            statsettings.efieldmetric);
        hasTractResults = false;

        if isfield(explorer.results, 'fiberfiltering') && ...
                isfield(calcsettings, 'fibfilt_connectome') && ...
                ~isempty(calcsettings.fibfilt_connectome)
            try
                connid = ea_unifiedmapping_conn2connid( ...
                    calcsettings.fibfilt_connectome);
                methodID = ea_unifiedmapping_method2methodid(explorer);
                connectomeResults = ...
                    explorer.results.fiberfiltering.(connid);

                if isfield(connectomeResults, methodID) && ...
                        isfield(connectomeResults.(methodID), 'fibsval')
                    [calculatedMaximum, hasTractResults] = ...
                        getCellMaximum( ...
                        connectomeResults.(methodID).fibsval);
                    if hasTractResults
                        tractLimitsRaw = [0, calculatedMaximum];
                    end
                end
            catch
                % Retain the metric-specific fallback limits for older or
                % incomplete explorer results.
            end
        end

        if hasTractResults
            [tractScale, tractExponent] = ...
                getEngineeringScale(tractLimitsRaw(2));
            tractLimitsRaw(2) = makeNiceUpperLimit( ...
                tractLimitsRaw(2), tractScale);
        else
            tractScale = 1;
            tractExponent = 0;
        end

        tractLimitsRaw(2) = ensureValidUpperLimit( ...
            tractLimitsRaw(1), tractLimitsRaw(2), tractScale);
        tractLimitsDisplayed = tractLimitsRaw / tractScale;

        tractThreshold = getSetting( ...
            statsettings, 'efieldthreshold_tract', 200);

        % Convert a previously stored sigmoid probability back to raw units.
        if tractThreshold > 0 && tractThreshold < 1
            tractThreshold = ea_EfieldFromSigmoid(tractThreshold);
        end

        tractThreshold = clampValue(tractThreshold, tractLimitsRaw);
        statsettings.efieldthreshold_tract = tractThreshold;

        sliderState.tract.used = true;
        sliderState.tract.limits = tractLimitsDisplayed;
        sliderState.tract.value = tractThreshold / tractScale;
        sliderState.tract.scale = tractScale;
        sliderState.tract.exponent = tractExponent;
        sliderState.tract.quantityName = quantityName;
        sliderState.tract.unitLabel = ...
            makeUnitLabel(tractExponent, unitName);

        metricName = lower(statsettings.efieldmetric);
        sliderState.tract.label = { ...
            ['Tracts are connected when ', metricName, ...
             ' ', quantityName], ...
            ['is > ', sprintf('%.1f', tractThreshold / tractScale), ...
             ' ', sliderState.tract.unitLabel]};

    case 'Sigmoid Field'
        probabilityLowerLimit = round( ...
            ea_SigmoidFromEfield(calcsettings.calcthreshold), 2);
        probabilityLowerLimit = min( ...
            max(probabilityLowerLimit, 0), 0.99);
        probabilityLimits = [probabilityLowerLimit, 1];

        %% Sweet spot
        spotThreshold = getSetting( ...
            statsettings, 'efieldthreshold_spot', 200);
        if spotThreshold > 1
            spotThreshold = ea_SigmoidFromEfield(spotThreshold);
        end
        spotThreshold = clampValue(spotThreshold, probabilityLimits);
        statsettings.efieldthreshold_spot = spotThreshold;

        sliderState.spot.used = true;
        sliderState.spot.limits = probabilityLimits;
        sliderState.spot.value = spotThreshold;
        sliderState.spot.scale = 1;
        sliderState.spot.exponent = 0;
        sliderState.spot.unitLabel = '';
        sliderState.spot.label = { ...
            'Voxels are connected when activation probability', ...
            ['is > ', sprintf('%.2f', spotThreshold)]};

        %% Fiber filtering
        tractThreshold = getSetting( ...
            statsettings, 'efieldthreshold_tract', 200);
        if tractThreshold > 1
            tractThreshold = ea_SigmoidFromEfield(tractThreshold);
        end
        tractThreshold = clampValue(tractThreshold, probabilityLimits);
        statsettings.efieldthreshold_tract = tractThreshold;

        sliderState.tract.used = true;
        sliderState.tract.limits = probabilityLimits;
        sliderState.tract.value = tractThreshold;
        sliderState.tract.scale = 1;
        sliderState.tract.exponent = 0;
        sliderState.tract.unitLabel = '';
        sliderState.tract.label = { ...
            'Tracts are connected when activation probability', ...
            ['is > ', sprintf('%.2f', tractThreshold)]};

    case 'VTA'
        % Threshold sliders are unused for binary VTA analyses. Preserve the
        % stored raw thresholds so changing models can restore them.
        sliderState.spot.used = false;
        sliderState.spot.limits = [];
        sliderState.spot.value = [];
        sliderState.spot.scale = 1;
        sliderState.spot.exponent = 0;
        sliderState.spot.unitLabel = '';
        sliderState.spot.label = ...
            'E-field threshold is not used for VTAs';

        sliderState.tract.used = false;
        sliderState.tract.limits = [];
        sliderState.tract.value = [];
        sliderState.tract.scale = 1;
        sliderState.tract.exponent = 0;
        sliderState.tract.unitLabel = '';
        sliderState.tract.label = ...
            'E-field threshold is not used for VTAs';

    otherwise
        error('Unknown stimulation model: %s', ...
            statsettings.stimulationmodel);
end

end


function value = getSetting(settings, fieldName, defaultValue)
% Return a valid finite scalar setting or its default.

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
% Restrict a scalar value to the supplied raw limits.

value = min(max(value, limits(1)), limits(2));

end


function limits = getFallbackTractLimits(metric)
% Provide usable limits before matching calculated results are available.

switch lower(metric)
    case 'sum'
        limits = [500, 5000];
    case 'mean'
        limits = [2.5, 1000];
    case {'peak', 'peak 5%'}
        limits = [50, 5000];
    otherwise
        limits = [50, 5000];
end

end


function [maximumValue, hasValues] = getCellMaximum(values)
% Find the largest finite positive value in an array or cell array.

maximumValue = [];
hasValues = false;

if ~iscell(values)
    values = {values};
end

for valueNumber = 1:numel(values)
    thisValue = full(values{valueNumber});
    thisValue = thisValue(isfinite(thisValue) & thisValue > 0);

    if ~isempty(thisValue)
        thisMaximum = max(thisValue(:));
        if ~hasValues || thisMaximum > maximumValue
            maximumValue = thisMaximum;
        end
        hasValues = true;
    end
end

end


function [scaleFactor, exponent] = getEngineeringScale(maximumValue)
% Choose an engineering-notation scale in powers of 1000.

if isempty(maximumValue) || ...
        ~isfinite(maximumValue) || maximumValue <= 100000
    exponent = 0;
    scaleFactor = 1;
    return;
end

exponent = 3 * floor(log10(maximumValue) / 3);
exponent = max(exponent, 0);
scaleFactor = 10^exponent;

end


function upperLimit = makeNiceUpperLimit(rawMaximum, scaleFactor)
% Round the displayed maximum upward to a convenient slider endpoint.

displayedMaximum = rawMaximum / scaleFactor;

if displayedMaximum <= 10
    niceDisplayedMaximum = ceil(displayedMaximum);
elseif displayedMaximum <= 100
    niceDisplayedMaximum = 10 * ceil(displayedMaximum / 10);
else
    niceDisplayedMaximum = 50 * ceil(displayedMaximum / 50);
end

upperLimit = niceDisplayedMaximum * scaleFactor;

end


function upperLimit = ensureValidUpperLimit(lowerLimit, upperLimit, scale)
% MATLAB sliders require strictly increasing limits.

if upperLimit <= lowerLimit
    upperLimit = lowerLimit + max(scale, abs(lowerLimit) * 0.1);
end

end


function [quantityName, unitName] = getFieldDescription(isPseudoM)
% Use generic input-map terminology only for pseudoM analyses.

if isPseudoM
    quantityName = 'input-map magnitude';
    unitName = 'input-map units';
else
    quantityName = 'E-field magnitude';
    unitName = 'V/m';
end

end


function label = makeUnitLabel(exponent, unitName)
% Format the scale and physical/generic unit for display.

if exponent == 0
    label = unitName;
else
    label = sprintf('x10^%d %s', exponent, unitName);
end

end
