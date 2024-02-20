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
                switch obj.statsettings.stattest
                    case 'N-Map'
                        id = 'plainconn';
                    otherwise
                        id = 'PAM_Ttest';
                end
            otherwise
                disp('The metric is not supported by PAM')
        end
    otherwise % Conventional / efield/VTA
        switch obj.statsettings.stimulationmodel
            case 'VTA'
                switch obj.statsettings.stattest
                    case 'N-Map' % do we even need an extra results entry for these?
                        id = 'plainconn';
                    otherwise
                        id = 'VAT_Ttest';
                end
            case {'Electric Field','Sigmoid Field'}  % E-fields
                id='efield';
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
