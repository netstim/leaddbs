function id=ea_method2methodid(obj,efm)

obj.compat_statmetric; % old compatibility for old statmetric notation (used to be stored as integers).
obj=ea_discfibers_compat_statmetrics2statsettings(obj);

if ~exist('efm','var')
    efm=obj.statsettings.efieldmetric;
end

switch obj.connectivity_type
    case 2 % PAM
        switch obj.statsettings.stimulationmodel
            case 'VTA'
                % plainconn is reserved for e-fields
                id = 'PAM_Ttest';
            case {'Sigmoid Field'}
                id = 'PAM_probA';
            otherwise
                disp('The metric is not supported by PAM')
        end
    otherwise % Conventional / efield/VTA
        switch obj.statsettings.stimulationmodel
            case 'VTA'
                switch obj.statsettings.stattest
                    case 'N-Map' % do we even need an extra results entry for these?
                        if strcmp(obj.e_field_metric,'Projection')
                            id = 'plainconn_proj';
                        else
                            id = 'plainconn';
                        end
                    otherwise
                        if strcmp(obj.e_field_metric,'Projection')
                            id = 'VAT_Ttest_proj';
                        else
                            id = 'VAT_Ttest';
                        end
                end
            case {'Electric Field','Sigmoid Field'}  % E-fields
                id='efield';
                if strcmp(obj.e_field_metric,'Projection')
                    switch efm
                        case 'Mean'
                            id=[id,'_proj_mean'];
                        case 'Peak'
                            id=[id,'_proj_peak'];
                        case 'Sum'
                            id=[id,'_proj_sum'];
                        case 'Peak 5%'
                            id=[id,'_proj_5peak'];
                    end
                else
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
                end
        end
end
