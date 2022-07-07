function id=ea_method2methodid(obj,efm)
if ~exist('efm','var')
    efm=obj.efieldmetric;
end
switch obj.connectivity_type
    case 2 % PAM
        switch obj.statmetric
            case {1,3,4}  % both refer to T-test atm
                id = 'PAM_Ttest'; 
            case 6
                id = 'plainconn';
            otherwise
                disp('The metric is not supported by PAM')
        end
    otherwise % E-field based
        switch obj.statmetric
            case {1,3,4}
                id = 'VAT_Ttest';
            case {2,5}
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
            case 6 % do we even need an extra results entry for these?
                id = 'plainconn';
        end
end
