function tractset = ea_discfibers_compat_statmetrics2statsettings(tractset)

if isnan(tractset.statmetric)
    return
end

switch tractset.statmetric
    case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)'}
        tractset.statsettings.stimulationmodel = 'VTA';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = '2-Sample Tests';
        tractset.statsettings.stattest = '2-Sample T-Test';
        tractset.statsettings.H0 = 'Zero';
    case {'One-Sample Tests / VTAs / PAM (OSS-DBS)'}
        tractset.statsettings.stimulationmodel = 'VTA';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = '1-Sample Tests';
        tractset.statsettings.stattest = '1-Sample T-Test';
        tractset.statsettings.H0 = 'Zero';
    case 'Correlations / E-fields (Irmen 2020)'
        tractset.statsettings.stimulationmodel = 'Electric Field';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = 'Correlations';
        tractset.statsettings.stattest = tractset.corrtype;
        tractset.statsettings.H0 = 'Zero';
    case 'Weighted Linear Regression / EF-Sigmoid (Dembek 2023)'
        tractset.statsettings.stimulationmodel = 'Sigmoid Field';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = '1-Sample Tests';
        tractset.statsettings.stattest = '1-Sample Weighted Regression';
        tractset.statsettings.H0 = 'Zero';
    case 'Odds Ratios / EF-Sigmoid (Jergas 2023)'
        tractset.statsettings.stimulationmodel = 'Sigmoid Field';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = '2-Sample Tests';
        tractset.statsettings.stattest = 'Odds Ratios';
        tractset.statsettings.H0 = 'Zero';
    case 'Proportion Test (Chi-Square) / VTAs (binary vars)'
        tractset.statsettings.stimulationmodel = 'VTA';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = 'Binary-Outcome Tests';
        tractset.statsettings.stattest = 'Proportion Test';
        tractset.statsettings.H0 = 'Zero';
    case 'Binomial Tests / VTAs (binary vars)'
        tractset.statsettings.stimulationmodel = 'VTA';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = 'Binary-Outcome Tests';
        tractset.statsettings.stattest = 'Binomial Test';
        tractset.statsettings.H0 = 'Zero';
    case 'Reverse T-Tests / E-Fields (binary vars)'
        tractset.statsettings.stimulationmodel = 'Electric Field';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = 'Binary-Outcome Tests';
        tractset.statsettings.stattest = 'Reverse 2-Sample T-Test';
        tractset.statsettings.H0 = 'Zero';
    case 'Plain Connections'
        tractset.statsettings.stimulationmodel = 'VTA';
        tractset.statsettings.efieldmetric = tractset.efieldmetric;
        tractset.statsettings.efieldthreshold = tractset.efieldthreshold;
        tractset.statsettings.connthreshold = tractset.connthreshold;
        tractset.statsettings.statfamily = 'Descriptive';
        tractset.statsettings.stattest = 'N-Map';
        tractset.statsettings.H0 = 'Zero';
end

% Rename results entries from spearman_* to efield_*:
connvals=fieldnames(tractset.results);
for connectome=1:length(connvals)
    if isfield(tractset.results.(connvals{connectome}),'spearman_sum')
        tractset.results.(connvals{connectome}).efield_sum=tractset.results.(connvals{connectome}).spearman_sum;
        tractset.results.(connvals{connectome})=rmfield(tractset.results.(connvals{connectome}),'spearman_sum');
        tractset.results.(connvals{connectome}).efield_peak=tractset.results.(connvals{connectome}).spearman_peak;
        tractset.results.(connvals{connectome})=rmfield(tractset.results.(connvals{connectome}),'spearman_peak');
        tractset.results.(connvals{connectome}).efield_mean=tractset.results.(connvals{connectome}).spearman_mean;
        tractset.results.(connvals{connectome})=rmfield(tractset.results.(connvals{connectome}),'spearman_mean');
        tractset.results.(connvals{connectome}).efield_5peak=tractset.results.(connvals{connectome}).spearman_5peak;
        tractset.results.(connvals{connectome})=rmfield(tractset.results.(connvals{connectome}),'spearman_5peak');
    end
end

tractset.statmetric=nan; % switch off legacy statmetric.
tractset.efieldmetric=nan; % switch off legacy
tractset.efieldthreshold=nan;
tractset.connthreshold=nan;
tractset.fileformatversion=1.1;
