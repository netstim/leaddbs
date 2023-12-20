function id=ea_explorer_method2methodid(obj,efm)
if ~exist('efm','var')
    efm=obj.statsettings.efieldmetric;
end


switch obj.statsettings.stimulationmodel
    case 'VTA'
        id='peak';
    case {'Electric Field','Sigmoid Field'}
        switch efm
            case 'Mean'
                id=['mean'];
            case 'Peak'
                id=['peak'];
            case 'Sum'
                id=['sum'];
            case 'Peak 5%'
                id=['peak5'];
        end
end
end
