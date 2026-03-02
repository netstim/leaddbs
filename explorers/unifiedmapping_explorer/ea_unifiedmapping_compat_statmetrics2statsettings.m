function explorer = ea_unifiedmapping_compat_statmetrics2statsettings(explorer)

if isnan(explorer.statmetric)
    return
end

switch explorer.statmetric
    case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)'}
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = '2-Sample Tests';
        explorer.statsettings.stattest = '2-Sample T-Test';
        explorer.statsettings.H0 = 'Zero';
    case {'One-Sample Tests / VTAs / PAM (OSS-DBS)'}
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = '1-Sample Tests';
        explorer.statsettings.stattest = '1-Sample T-Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Correlations / E-fields (Irmen 2020)'
        explorer.statsettings.stimulationmodel = 'Electric Field';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = 'Correlations';
        explorer.statsettings.stattest = explorer.corrtype;
        explorer.statsettings.H0 = 'Zero';
    case 'Weighted Linear Regression / EF-Sigmoid (Dembek 2023)'
        explorer.statsettings.stimulationmodel = 'Sigmoid Field';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = '1-Sample Tests';
        explorer.statsettings.stattest = '1-Sample Weighted Regression';
        explorer.statsettings.H0 = 'Zero';
    case 'Odds Ratios / EF-Sigmoid (Jergas 2023)'
        explorer.statsettings.stimulationmodel = 'Sigmoid Field';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = '2-Sample Tests';
        explorer.statsettings.stattest = 'Odds Ratios';
        explorer.statsettings.H0 = 'Zero';
    case 'Proportion Test (Chi-Square) / VTAs (binary vars)'
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = 'Binary-Outcome Tests';
        explorer.statsettings.stattest = 'Proportion Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Binomial Tests / VTAs (binary vars)'
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = 'Binary-Outcome Tests';
        explorer.statsettings.stattest = 'Binomial Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Reverse T-Tests / E-Fields (binary vars)'
        explorer.statsettings.stimulationmodel = 'Electric Field';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
        explorer.statsettings.statfamily = 'Binary-Outcome Tests';
        explorer.statsettings.stattest = 'Reverse 2-Sample T-Test';
        explorer.statsettings.H0 = 'Zero';
    case 'Plain Connections'
        explorer.statsettings.stimulationmodel = 'VTA';
        explorer.statsettings.efieldmetric = explorer.efieldmetric;
        explorer.statsettings.efieldthreshold = explorer.efieldthreshold;
        explorer.statsettings.connthreshold = explorer.connthreshold;
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

