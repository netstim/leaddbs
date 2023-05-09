function id=ea_method2methodid(obj,efm)
if ~exist('efm','var')
    efm=obj.efieldmetric;
end

    obj.compat_statmetric; % compatibility for old statmetric notation (used to be stored as integers).

switch obj.connectivity_type
    case 2 % PAM
        switch obj.statmetric
            case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)','Proportion Test (Chi-Square) / VTAs (binary vars)','Binomial Tests / VTAs (binary vars)'}  % both refer to T-test atm
                id = 'PAM_Ttest'; 
            case 'Plain Connections'
                id = 'plainconn';
            otherwise
                disp('The metric is not supported by PAM')
        end
    otherwise
        switch obj.statmetric
            case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)','One-Sample Tests / VTAs / PAM (OSS-DBS)','Proportion Test (Chi-Square) / VTAs (binary vars)','Binomial Tests / VTAs (binary vars)'}  % VTAs
                id = 'VAT_Ttest';
            case {'Correlations / E-fields (Irmen 2020)','Reverse T-Tests / E-Fields (binary vars)'}  % E-fields
                id='spearman';
                switch efm
                    case 'Mean'
                        id=[id,'_mean'];
                    case 'Peak'
                        id=[id,'_peak'];
                    case 'Sum'
                        id=[id,'_sum'];
                    case 'Peak 5%'
                        id=[id,'_5peak'];
                end
            case 'Plain Connections' % do we even need an extra results entry for these?
                id = 'plainconn';
        end
end
