function explorer = ea_unifiedmapping_compat_statmetrics2statsettings(explorer)
% Migrate the old shared E-field threshold

legacyThreshold = [];

% Newer legacy files may contain the old value inside statsettings.
if isfield(explorer.statsettings, 'efieldthreshold') && ...
        ~isempty(explorer.statsettings.efieldthreshold)

    legacyThreshold = ...
        explorer.statsettings.efieldthreshold;

% Older objects may contain it as a top-level property.
elseif isobject(explorer) && ...
        isprop(explorer, 'efieldthreshold') && ...
        ~isempty(explorer.efieldthreshold)

    legacyThreshold = explorer.efieldthreshold;

% Support struct-form legacy data as well.
elseif isstruct(explorer) && ...
        isfield(explorer, 'efieldthreshold') && ...
        ~isempty(explorer.efieldthreshold)

    legacyThreshold = explorer.efieldthreshold;
end

if isempty(legacyThreshold)
    legacyThreshold = 200;
end

% The old E-field threshold applies to sweet spot and tract mapping.
% Never overwrite already-migrated tool-specific values.
if ~isfield(explorer.statsettings, 'efieldthreshold_spot') || ...
        isempty(explorer.statsettings.efieldthreshold_spot)

    explorer.statsettings.efieldthreshold_spot = ...
        legacyThreshold;
end

if ~isfield(explorer.statsettings, 'efieldthreshold_tract') || ...
        isempty(explorer.statsettings.efieldthreshold_tract)

    explorer.statsettings.efieldthreshold_tract = ...
        legacyThreshold;
end

% Network mapping has no equivalent legacy E-field threshold.
if ~isfield(explorer.statsettings, 'efieldthreshold_network') || ...
        isempty(explorer.statsettings.efieldthreshold_network)

    explorer.statsettings.efieldthreshold_network = 50;
end

%% Resolve legacy E-field metric

if isobject(explorer) && isprop(explorer, 'efieldmetric') && ...
        ~isempty(explorer.efieldmetric)

    legacyEfieldMetric = explorer.efieldmetric;

elseif isstruct(explorer) && ...
        isfield(explorer, 'efieldmetric') && ...
        ~isempty(explorer.efieldmetric)

    legacyEfieldMetric = explorer.efieldmetric;

elseif isfield(explorer.statsettings, 'efieldmetric') && ...
        ~isempty(explorer.statsettings.efieldmetric)

    legacyEfieldMetric = ...
        explorer.statsettings.efieldmetric;
else
    legacyEfieldMetric = 'Sum';
end

%% Resolve legacy connection threshold

if isobject(explorer) && isprop(explorer, 'connthreshold') && ...
        ~isempty(explorer.connthreshold)

    legacyConnThreshold = explorer.connthreshold;

elseif isstruct(explorer) && ...
        isfield(explorer, 'connthreshold') && ...
        ~isempty(explorer.connthreshold)

    legacyConnThreshold = explorer.connthreshold;

elseif isfield(explorer.statsettings, 'connthreshold') && ...
        ~isempty(explorer.statsettings.connthreshold)

    legacyConnThreshold = ...
        explorer.statsettings.connthreshold;
else
    legacyConnThreshold = 20;
end

if isnan(explorer.statmetric)
    return
end

switch explorer.statmetric
    case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)'}
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = '2-Sample Tests';
        explorer.statsettings.stattest = '2-Sample T-Test';
        explorer.statsettings.H0 = 'Zero';
    case {'One-Sample Tests / VTAs / PAM (OSS-DBS)'}
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = '1-Sample Tests';
        explorer.statsettings.stattest = '1-Sample T-Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Correlations / E-fields (Irmen 2020)'
        explorer.statsettings.stimulationmodel = 'Electric Field';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = 'Correlations';
        explorer.statsettings.stattest = explorer.corrtype;
        explorer.statsettings.H0 = 'Zero';
    case 'Weighted Linear Regression / EF-Sigmoid (Dembek 2023)'
        explorer.statsettings.stimulationmodel = 'Sigmoid Field';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = '1-Sample Tests';
        explorer.statsettings.stattest = '1-Sample Weighted Regression';
        explorer.statsettings.H0 = 'Zero';
    case 'Odds Ratios / EF-Sigmoid (Jergas 2023)'
        explorer.statsettings.stimulationmodel = 'Sigmoid Field';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = '2-Sample Tests';
        explorer.statsettings.stattest = 'Odds Ratios';
        explorer.statsettings.H0 = 'Zero';
    case 'Proportion Test (Chi-Square) / VTAs (binary vars)'
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = 'Binary-Outcome Tests';
        explorer.statsettings.stattest = 'Proportion Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Binomial Tests / VTAs (binary vars)'
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = 'Binary-Outcome Tests';
        explorer.statsettings.stattest = 'Binomial Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Reverse T-Tests / E-Fields (binary vars)'
        explorer.statsettings.stimulationmodel = 'Electric Field';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
        explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = 'Binary-Outcome Tests';
        explorer.statsettings.stattest = 'Reverse 2-Sample T-Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Plain Connections'
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = legacyEfieldMetric;
    explorer.statsettings.connthreshold = legacyConnThreshold;
        explorer.statsettings.statfamily = 'Descriptive';
        explorer.statsettings.stattest = 'N-Map';
        explorer.statsettings.H0 = 'Zero';
end

% Rename results entries from spearman_* to efield_*:
%these results are specific to fiberfiltering
connvals=fieldnames(explorer.results);
for connectome=1:length(connvals)
    if isfield(explorer.results.(connvals{connectome}),'spearman_sum')
        explorer.results.fiberfiltering.(connvals{connectome}).efield_sum=explorer.results.fiberfiltering.(connvals{connectome}).spearman_sum;
        explorer.results.fiberfiltering.(connvals{connectome})=rmfield(explorer.results.fiberfiltering.(connvals{connectome}),'spearman_sum');
        explorer.results.fiberfiltering.(connvals{connectome}).efield_peak=explorer.results.fiberfiltering.(connvals{connectome}).spearman_peak;
        explorer.results.fiberfiltering.(connvals{connectome})=rmfield(explorer.results.fiberfiltering.(connvals{connectome}),'spearman_peak');
        explorer.results.fiberfiltering.(connvals{connectome}).efield_mean=explorer.results.fiberfiltering.(connvals{connectome}).spearman_mean;
        explorer.results.fiberfiltering.(connvals{connectome})=rmfield(explorer.results.fiberfiltering.(connvals{connectome}),'spearman_mean');
        explorer.results.fiberfiltering.(connvals{connectome}).efield_5peak=explorer.results.fiberfiltering.(connvals{connectome}).spearman_5peak;
        explorer.results.fiberfiltering.(connvals{connectome})=rmfield(explorer.results.fiberfiltering.(connvals{connectome}),'spearman_5peak');
    end
end

