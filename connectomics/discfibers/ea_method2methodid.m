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
            case 'mean'
                id=[id,'_mean'];
            case 'peak'
                id=[id,'_peak'];
            case 'sum'
                id=[id,'_sum'];
            case '5peak'
                id=[id,'_5peak'];
        end
end
