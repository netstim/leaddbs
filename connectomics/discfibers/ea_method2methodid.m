function id=ea_method2methodid(obj,efm)
if ~exist('efm','var')
    efm=obj.efieldmetric;
end
switch obj.statmetric
    case 1
        id='ttests';
    case 2
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
    case 3
        id = 'plainconn';
end
